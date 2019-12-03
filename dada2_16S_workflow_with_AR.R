### Start script by running Rscript dada2_v1.10_workflow.R [Input folder with fastq files] [ID for run] [output folder]
### The input folder must contain forward and reverse read for the samples 
### if no RunID is given, it is generated from the standard folder structure
### if no output folder is given, standard folder defined in "outbase_std" is used
### TODO: create profiles (or automatic detection?) of V1/V2, V3/V4, Archaea, ITS, etc.

args <- commandArgs(trailingOnly = TRUE)

outbase_std="/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10"

path=gsub("/$","",args[1])
if(is.na(args[2])){runid=strsplit(rev(strsplit(path, split="/")[[1]])[2],split="_")[[1]][1]}else{runid=args[2]}
if(is.na(args[3])){outbase=outbase_std}else{outbase=gsub("/$","",args[3])}

numbases=1e+09
threads=8

.libPaths("/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/R_Librarires/")

library(dada2)
version=packageVersion("dada2")
set.seed(666)
cat(paste0("Running DADA2 version ",version,"\n"))
library(ShortRead)


fns <- list.files(path)
#fns

####################
### load your data
####################

fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
### make sure that R1 is for forward read and R2 for reverse

fnFs <- fastqs[grepl("R1_001.fastq.gz", fastqs)] ## Just the forward read files
fnRs <- fastqs[grepl("R2_001.fastq.gz", fastqs)] ## Just the reverse read files

## Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, "(-L1)*_((S[0-9]+)|([ACTG-]+))_L001_R1_001.fastq.gz"), `[`, 1)
cat(paste0("Starting processing for ",length(sample.names)," samples\n"))

## Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)


###############################################
### check if V3/V4 (351F - 806R) library
##############################################

cat("Checking if run is V3V4 sequencing run")

FWD="CCTACGGGAGGCAGCAG"
REV="GGACTACHVGGGTWTCTAAT"

allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn))[1:2000], fixed = FALSE)
    return(sum(nhits > 0))
}


i=1

nreads=length(sread(readFastq(fnFs[[i]])))
while(nreads<2500){i=i+1; nreads=length(sread(readFastq(fnFs[[i]])))}


v3v4=rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[i]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[i]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[i]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[i]]))/2000


is_v3v4=any(v3v4>.2)

if(is_v3v4){

if(outbase==outbase_std){outbase=paste0(outbase_std,"_V3V4")}

outdir=paste0(outbase,"/",runid)

cat(paste0("All output will be saved in: ",outdir,"\n"))

dir.create(outdir,recursive=T,showWarnings=F)
dir.create(paste0(outdir,"/plots"),recursive=T,showWarnings=F)
dir.create(paste0(outdir,"/errors"),recursive=T,showWarnings=F)


cutadapt <- "/home/sukmb276/Isilon/Microbiome/clean_data_from_dada2/DADA2_versions/cutadapt_env/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(outdir, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
			     "-j", threads, "--discard-untrimmed", 
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

fnFs=fnFs.cut
fnRs=fnRs.cut
}


}else{
outdir=paste0(outbase,"/",runid)

cat(paste0("All output will be saved in: ",outdir,"\n"))

dir.create(outdir,recursive=T,showWarnings=F)
dir.create(paste0(outdir,"/plots"),recursive=T,showWarnings=F)
dir.create(paste0(outdir,"/errors"),recursive=T,showWarnings=F)

}


###########################################
## Examine quality profiles of F & R reads
###########################################
pdf(paste0(outdir,"/plots/plotQualityProfile.pdf"), onefile=T)
plotQualityProfile(fnFs[1:2]) ## remove 20 plus 10 (primers and first odd bases)
plotQualityProfile(fnRs[1:2]) ## remove 20 plus 10 (primers and first odd bases)
dev.off()


##################################
## Perform filtering and trimming
##################################
filt_path <- file.path(outdir, "filtered") # Place filtered files in filtered/subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))


## Filter  the forward and reverse reads:
## Important to remove primers and low quality regions
if(is_v3v4){
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(265,245), ## Different settings tried,these are good for current primer constructs for MiSeq and HiSeq
                     maxN=0, maxEE=c(2,2), truncQ=5, rm.phix=TRUE,
                     compress=TRUE, multithread=threads) #
}else{
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,180), ## Different settings tried,these are good for current primer constructs for MiSeq and HiSeq
                     trimLeft=c(5, 5),
                     maxN=0, maxEE=c(2,2), truncQ=5, rm.phix=TRUE,
                     compress=TRUE, multithread=threads) #
}

## Examine quality profiles of filtered reads
pdf(paste0(outdir,"/plots/plotQualityProfile.filt.pdf"), onefile=T)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])
dev.off()


#########################
## Learn the Error Rates
#########################
## Learn forward error rates
errF <- learnErrors(filtFs, nbases=numbases, multithread=threads,randomize=T)
## Learn reverse error rates
errR <- learnErrors(filtRs, nbases=numbases, multithread=threads,randomize=T)

saveRDS(errF, paste0(outdir,"/errors/errF.Rds"))
saveRDS(errR, paste0(outdir,"/errors/errR.Rds"))


## Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

## Plot estimated error as sanity check
pdf(paste0(outdir,"/plots/plotErrors_F.pdf"), onefile=T)
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf(paste0(outdir,"/plots/plotErrors_R.pdf"), onefile=T)
plotErrors(errR, nominalQ=TRUE)
dev.off()


#########################
## Perform dereplication
#########################
## Dereplicate the filtered fastq files
derepRs <- derepFastq(filtRs, verbose=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


####################
## Sample Inference
####################
## Apply the core sequence-variant inference algorithm to the dereplicated data
## Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=threads)
dadaRs <- dada(derepRs, err=errR, multithread=threads)


## Inspect the dada-class object returned by dada
#dadaFs[[1]]

## Merge the denoised forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=10)
## Inspect the merger data.frame from the first sample
#head(mergers[[1]])


############################
## Construct sequence table
############################
seqtab <- makeSequenceTable(mergers)
## Get dimensions
#dim(seqtab)

## Inspect distribution of sequence lengths
#table(nchar(getSequences(seqtab)))


saveRDS(seqtab, paste0(outdir,"/seqtab.Rds"))

###################
## Remove chimeras
###################
## Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads, verbose=TRUE)
#dim(seqtab.nochim)
#sum(seqtab.nochim)/sum(seqtab)

saveRDS(seqtab.nochim, paste0(outdir,"/seqtab_nochim.Rds"))

####################################
## Track reads through the pipeline
####################################
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
#head(track)

write.table(track, paste0(outdir,"/track_reads.txt"),sep="\t",quote=FALSE)

###################
## Assign taxonomy
###################
taxHS <- assignTaxonomy(seqtab.nochim, "/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/rdp_train_set_16.fa.gz", multithread=threads) ## CHANGE to directory and pertinent database
colnames(taxHS) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
#unname(head(taxHS))
#unname(tail(taxHS))


write.table(taxHS, file = paste0(outdir,"/taxa_SV.tsv"), quote=FALSE)


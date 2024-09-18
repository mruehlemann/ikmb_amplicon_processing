### Start script by running Rscript dada2_v1.10_workflow.R [Input folder with fastq files] [ID for run] [output folder]
### The input folder must contain forward and reverse read for the samples 
### if no RunID is given, it is generated from the standard folder structure
### if no output folder is given, standard folder defined in "outbase_std" is used
### TODO: create profiles (or automatic detection?) of V1/V2, V3/V4, Archaea, ITS, etc.

.libPaths("/work_beegfs/sukmb276/Microbiome/clean_data_from_dada2/R_libraries_new/1.14")

numbases=1e+09
threads=24

args <- commandArgs(trailingOnly = TRUE)

profiles=c("V1V2","V3V4","PacBio","Archaea","Fungi")

profile=args[1]
if(! profile %in% profiles){
cat(paste0("Please specify one of these profiles:", paste(profiles, collapse=", "),"\n" ))
quit()
}

path=gsub("/$","",args[2])
if(is.na(args[3])){runid=strsplit(rev(strsplit(path, split="/")[[1]])[2],split="_")[[1]][1]}else{runid=args[3]}


if(profile=="V1V2"){
	outbase_std="/work_beegfs/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10"
	truncLens=c(230,180)
	trims=c(5,5)
	maxEEs=c(2,2)
	taxdatabase="/work_beegfs/sukmb276/Microbiome/clean_data_from_dada2/reference_data/GTDB_bac120_arc122_ssu_r202_Genus.fa.gz"
} else if (profile=="V3V4"){
	outbase_std="/work_beegfs/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10_V3V4"
	FWD="CCTACGGGAGGCAGCAG"
	REV="GGACTACHVGGGTWTCTAAT"
	truncLens=c(260,200)
	trims=c(5,5)
	maxEEs=c(2,2)
	taxdatabase="/work_beegfs/sukmb276/Microbiome/clean_data_from_dada2/reference_data/GTDB_bac120_arc122_ssu_r202_Genus.fa.gz"
} else if (profile=="PacBio"){
	outbase_std="/work_beegfs/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10_PacBio"
	FWD="AGRGTTYGATYMTGGCTCAG"
	REV="RGYTACCTTGTTACGACTT"
	taxdatabase="/work_beegfs/sukmb276/Microbiome/clean_data_from_dada2/reference_data/GTDB_bac120_arc122_ssu_r202_Genus.fa.gz"
} else if (profile=="Archaea"){
	outbase_std="/work_beegfs/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10_Archaea"
	FWD="CAGCMGCCGCGGTAA"
	REV="GGACTACVSGGGTATCTAAT"
	truncLens=c(250,200)
	trims=c(5,5)
	maxEEs=c(2,2)
	taxdatabase="/work_beegfs/sukmb276/Microbiome/clean_data_from_dada2/reference_data/GTDB_bac120_arc122_ssu_r202_Genus.fa.gz"
} else if (profile=="Fungi"){
	outbase_std="/work_beegfs/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10_Fungi"
	truncLens=c(230,150)
	trims=c(5,5)
	maxEEs=c(4,4)
	taxdatabase="/work_beegfs/sukmb276/Microbiome/clean_data_from_dada2/reference_data/sh_general_release_10.05.2021/sh_general_release_dynamic_10.05.2021.fasta"
}

if(is.na(args[4])){outbase=outbase_std}else{outbase=gsub("/$","",args[4])}
if(is.na(args[5])){assignTax=TRUE}else{assignTax=FALSE}
if(is.na(args[6])){doverbose=FALSE}else{doverbose=TRUE}

cat(paste0("Running DADA2 Pipeline in ", profile, " profile\n"))

cat("Loading packages...\n")
suppressPackageStartupMessages(library(dada2))
version=packageVersion("dada2")
set.seed(666)
cat(paste0("Running DADA2 version ",version,"\n"))
suppressPackageStartupMessages(library(ShortRead))
#suppressPackageStartupMessages(library(Biostrings))
#suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(reshape2))



### CREATE OUTPUT FOLDERS ETC.

outdir=paste0(outbase,"/",runid)

cat(paste0("All output will be saved in: ",outdir,"\n"))

dir.create(outdir,recursive=T,showWarnings=F)
dir.create(paste0(outdir,"/plots"),recursive=T,showWarnings=F)
dir.create(paste0(outdir,"/errors"),recursive=T,showWarnings=F)


### CREATE LIST FOR COUNTING

countlist=vector("list", 0)

####################
### load your data
####################

fns <- list.files(path)
#fns

fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
### make sure that R1 is for forward read and R2 for reverse

if(profile != "PacBio"){

fnFs <- fastqs[grepl("R1_001.fastq.gz", fastqs)] ## Just the forward read files
fnRs <- fastqs[grepl("R2_001.fastq.gz", fastqs)] ## Just the reverse read files

## Get sample names from the first part of the forward read filenames
#sample.names <- sapply(strsplit(fnFs, "(-L1)*_((S[0-9]+)|([ACTG-]+))_L001_R1_001.fastq.gz"), `[`, 1)
sample.names <- sapply(sapply(strsplit(fnFs, "(-L1)*_((S[0-9]+)|([ACTG-]+))_L001_R1_001.fastq.gz"), `[`, 1), function(x) strsplit(x, split="_")[[1]][2])
sample.names.all = sample.names 
cat(paste0("Starting processing for ",length(sample.names)," samples\n"))

## Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

} else {
fastqs = fastqs[grepl("i7", fastqs) & grepl("i5",fastqs)]

#sample.names = gsub("[.]hifi_reads[.]fastq[.]gz$", "", gsub("demultiplex[.]", "", basename(fastqs)))
sample.names = gsub("[.]fastq.gz$","", basename(fastqs))
sample.names.all = sample.names
fnFs=fastqs
fnFs <- file.path(path, fnFs)
}


countlist[[length(countlist)+1]] = data.frame(sample=sample.names, raw=ShortRead::countLines(fnFs)/4)


#####


if (profile == "PacBio"){

cat("Profile is", profile,"\nPrimer sequences will be removed using removePrimers() function\n")

path.cut <- file.path(outdir, "noprimers")
if(!dir.exists(path.cut)) dir.create(path.cut)

library(parallel)
fnFs.cut <- file.path(path.cut, paste0(sample.names, ".noprimer.fastq.gz"))
prime=mclapply(seq_along(fnFs), function(x) removePrimers(fnFs[x], fnFs.cut[x], primer.fwd=FWD, primer.rev=dada2:::rc(REV), orient=TRUE, verbose=doverbose), mc.cores=threads)

fnFs = fnFs.cut[file.exists(fnFs.cut)]

sample.names.remain=sample.names[file.exists(fnFs.cut)]
countlist[[length(countlist)+1]] = data.frame(sample=sample.names.remain, primercut=ShortRead::countLines(fnFs)/4)

filt_path <- file.path(outdir, "filtered") 
filtFs <- file.path(filt_path, paste0(sample.names.remain, ".filtered.fastq.gz"))

track <- filterAndTrim(fnFs, filtFs, minQ=3, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2, multithread = threads, qualityType="FastqQuality")

filtFs.remain = filtFs[file.exists(filtFs)]
sample.names.remain = sample.names.remain[file.exists(filtFs)]
countlist[[length(countlist)+1]] = data.frame(sample=sample.names.remain, qc=ShortRead::countLines(filtFs.remain)/4)
 
derepFs <- derepFastq(filtFs.remain, verbose=TRUE, qualityType="FastqQuality")              

errF <- learnErrors(derepFs, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=threads, qualityType="FastqQuality", nbases=numbases, randomize=T, verbose=doverbose)

saveRDS(errF, paste0(outdir,"/errors/errF.Rds"))

plotErrors(errF)

mergers <- dada(derepFs, err=errF, BAND_SIZE=32, multithread=threads)
countlist[[length(countlist)+1]] = data.frame(sample=sample.names.remain, denoised=sapply(mergers, function(x) sum(getUniques(x))))

saveRDS(mergers, paste0(outdir,"/errors/mergers.Rds"))              
                     
} else {


if(profile %in% c("V3V4","Archaea")){

cutadapt <- "/work_beegfs/sukmb276/Microbiome/clean_data_from_dada2/R_libraries_new/1.14/dada2_env/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
ca_ver=system2(cutadapt, args = "--version") # Run shell commands from R

cat("Profile is", profile,"\nPrimer sequences will be removed using cutadapt version",ca_ver,"\n")

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
  if(doverbose==F){cat("Processing Sample", i, "of", length(fnFs),"\r")} 
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-j", threads, "--discard-untrimmed", "--minimum-length", "100",
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i]), 
                             stdout=ifelse(doverbose==TRUE, "/dev/stdout","/dev/null"), 
                             stderr=ifelse(doverbose==TRUE, "/dev/stderr","/dev/null")) # input files
}
cat("\n")

#rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
#    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
#    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
#    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

fnFs=fnFs.cut
fnRs=fnRs.cut

countlist[[length(countlist)+1]] = data.frame(sample=sample.names, primercut=ShortRead::countLines(fnFs)/4)

} 


###########################################
## Examine quality profiles of F & R reads
###########################################
pdf(paste0(outdir,"/plots/plotQualityProfile.pdf"), onefile=T)
capture.output(plotQualityProfile(fnFs[1:2]))
capture.output(plotQualityProfile(fnRs[1:2]))
dev.off()


##################################
## Perform filtering and trimming
##################################
filt_path <- file.path(outdir, "filtered") # Place filtered files in filtered/subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))


## Filter  the forward and reverse reads:
## Important to remove primers and low quality regions
cat("Running quality control of sequencing data...\n")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncLens,
                     trimLeft=trims,
                     maxN=0, maxEE=maxEEs, truncQ=5, rm.phix=TRUE,
                     compress=TRUE, multithread=threads) #



exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]
sample.names=sample.names[exists]

## Examine quality profiles of filtered reads
pdf(paste0(outdir,"/plots/plotQualityProfile.filt.pdf"), onefile=T)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])
dev.off()


countlist[[length(countlist)+1]] = data.frame(sample=sample.names, filtered=ShortRead::countLines(filtFs)/4)

#########################
## Learn the Error Rates
#########################
cat("Learning error rates...\n")
## Learn forward error rates
errF <- learnErrors(filtFs, nbases=numbases, multithread=threads,randomize=T, verbose=doverbose)
## Learn reverse error rates
errR <- learnErrors(filtRs, nbases=numbases, multithread=threads,randomize=T, verbose=doverbose)

saveRDS(errF, paste0(outdir,"/errors/errF.Rds"))
saveRDS(errR, paste0(outdir,"/errors/errR.Rds"))


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
cat("Running dereplication...\n")

## Dereplicate the filtered fastq files
derepRs <- derepFastq(filtRs, verbose=doverbose)
derepFs <- derepFastq(filtFs, verbose=doverbose)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


####################
## Sample Inference
####################
## Apply the core sequence-variant inference algorithm to the dereplicated data
## Infer the sequence variants in each sample
cat("Running dada2 ASV inference...\n")

dadaFs <- dada(derepFs, err=errF, multithread=threads, verbose=doverbose)
dadaRs <- dada(derepRs, err=errR, multithread=threads, verbose=doverbose)

countlist[[length(countlist)+1]] = data.frame(sample=sample.names, denoised=sapply(dadaFs, function(x) sum(getUniques(x))))

if(profile=="Fungi"){
### Special case for fungi: if merging doesn't work, simply join FWD and REV read
merger <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, returnRejects=TRUE, verbose=doverbose)
concat <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE, verbose=doverbose)

mergers=merger
for(x in sample.names){m=merger[[x]]; c=concat[[x]]; m[!m$accept,] <- c[!m$accept,]; mergers[[x]] <- m}

}else{
## Merge the denoised forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=doverbose)
}

countlist[[length(countlist)+1]] = data.frame(sample=sample.names, merged=sapply(mergers, function(x) sum(getUniques(x))))

}


############################
## Construct sequence table
############################
cat("Constructing sequence table...\n")
seqtab <- makeSequenceTable(mergers)
countlist[[length(countlist)+1]] = data.frame(sample=sapply(rownames(seqtab), function(x) strsplit(x, split="[.]")[[1]][1]), tabled=rowSums(seqtab))

saveRDS(seqtab, paste0(outdir,"/seqtab.Rds"))
cat("Seqtab written to:", paste0(outdir,"/seqtab.Rds"), "\n")

###################
## Remove chimeras
###################
## Remove chimeric sequences:
cat("Removing chimeric sequences...\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads, verbose=doverbose)
countlist[[length(countlist)+1]] = data.frame(sample=sapply(rownames(seqtab.nochim), function(x) strsplit(x, split="[.]")[[1]][1]), nonchim=rowSums(seqtab.nochim))

saveRDS(seqtab.nochim, paste0(outdir,"/seqtab_nochim.Rds"))
cat("Chimera-free seqtab written to:", paste0(outdir,"/seqtab_nochim.Rds"), "\n")

asv.names = data.frame(name = paste0("ASV_",sprintf("%06d", seq_along(colnames(seqtab.nochim)))), seq = colnames(seqtab.nochim), stringsAsFactors = F)


if(assignTax == TRUE){
###################
## Assign taxonomy
###################
cat("Assigning taxonomy using:", taxdatabase, "\n")
taxHS <- assignTaxonomy(seqs=seqtab.nochim, refFasta=taxdatabase, multithread=threads, tryRC=ifelse(profile=="Fungi",TRUE,FALSE), verbose=doverbose) ## CHANGE to directory and pertinent database
colnames(taxHS) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")[1:length(colnames(taxHS))]
write.table(taxHS, file = paste0(outdir,"/taxa_SV.tsv"), quote=FALSE)
cat("Taxonomic annotations written to:",  paste0(outdir,"/taxa_SV.tsv"), "\n")
}

####################################
## Track reads through the pipeline
####################################

trackreads <- data.frame(do.call("cbind", lapply(countlist, function(x) data.frame(row.names=sample.names.all, x[match(sample.names.all, x$sample),2, drop=F]))))
write.table(trackreads, paste0(outdir,"/track_reads.txt"),sep="\t",quote=FALSE)

trackreads.rel <- t(apply(trackreads, 1, function(x) x/x[1]))
write.table(trackreads.rel, paste0(outdir,"/track_reads_rel.txt"),sep="\t",quote=FALSE)

system2("chmod", args=c("-R 777", outdir))

cat("Pipeline successfully completed!\n\n")
################################
### FOR APP?!
###########################

#st.final.renamed <- seqtab.nochim
#colnames(st.final.renamed) <- asv.names$name

#prefilt = sapply(fnFs, get_stats, simplify = F)
#postfilt = sapply(filtFs, get_stats, simplify = F)
#names(prefilt) <- names(postfilt) <- sample.names



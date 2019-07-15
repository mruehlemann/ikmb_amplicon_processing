### Start script by running Rscript dada2_v1.10_workflow.R [Input folder with fastq files] [ID for run] [output folder]
### The input folder must contain forward and reverse read for the samples 
### if no RunID is given, it is generated from the standard folder structure
### if no output folder is given, standard folder defined in "outbase_std" is used
### TODO: create profiles (or automatic detection?) of V1/V2, V3/V4, Archaea, ITS, etc.

args <- commandArgs(trailingOnly = TRUE)

outbase_std="/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10_Fungi"

path=gsub("/$","",args[1])
if(is.na(args[2])){runid=strsplit(rev(strsplit(path, split="/")[[1]])[2],split="_")[[1]][1]}else{runid=args[2]}
if(is.na(args[3])){outbase=outbase_std}else{outbase=gsub("/$","",args[3])}

numbases=1e+09
threads=16

outdir=paste0(outbase,"/",runid)

cat(paste0("All output will be saved in: ",outdir,"\n"))

dir.create(outdir,recursive=T,showWarnings=F)
dir.create(paste0(outdir,"/plots"),recursive=T,showWarnings=F)
dir.create(paste0(outdir,"/errors"),recursive=T,showWarnings=F)


.libPaths("/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/R_Librarires/")

library(dada2)
version=packageVersion("dada2")
set.seed(666)
cat(paste0("Running DADA2 version ",version,"\n"))

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
sample.names <- sapply(strsplit(fnFs, "-L1_S[0-9]+_L001_R1_001.fastq.gz"), `[`, 1) ## check if it is 1 or 2
cat(paste0("Starting processing for ",length(sample.names)," samples\n"))

## Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)


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
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,150), ## Different settings tried,these are good for current primer constructs for MiSeq and HiSeq
                     trimLeft=c(5, 5),
                     maxN=0, maxEE=c(4,4), truncQ=5, rm.phix=TRUE,
                     compress=TRUE, multithread=threads) #


## Examine quality profiles of filtered reads
pdf(paste0(outdir,"/plots/plotQualityProfile.filt.pdf"), onefile=T)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])
dev.off()

exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

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

###
sample.names=sample.names[exists]

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

##### Also table only fwd and rev due to bad quality
seqtab.Fs <- makeSequenceTable(dadaFs)
saveRDS(seqtab.Fs, paste0(outdir,"/seqtab_Fs.Rds"))
seqtab.Fs.nochim <- removeBimeraDenovo(seqtab.Fs, method="consensus", multithread=threads, verbose=TRUE)
saveRDS(seqtab.Fs.nochim, paste0(outdir,"/seqtab_Fs_nochim.Rds"))

seqtab.Rs <- makeSequenceTable(dadaRs)
saveRDS(seqtab.Rs, paste0(outdir,"/seqtab_Rs.Rds"))
seqtab.Rs.nochim <- removeBimeraDenovo(seqtab.Rs, method="consensus", multithread=threads, verbose=TRUE)
saveRDS(seqtab.Rs.nochim, paste0(outdir,"/seqtab_Rs_nochim.Rds"))




####################################
## Track reads through the pipeline
####################################
getN <- function(x) sum(getUniques(x))
track <- cbind(out[exists,], sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
head(track)
dim(track)
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
#head(track)

write.table(track, paste0(outdir,"/track_reads.txt"),sep="\t",quote=FALSE)

###################
## Assign taxonomy
###################
taxHS <- assignTaxonomy(seqtab.nochim, "/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/sh_general_release_dynamic_10.10.2017_dev.fasta", multithread=threads,tryRC=T) ## CHANGE to directory and pertinent database
taxHS.Fs <- assignTaxonomy(seqtab.Fs.nochim, "/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/sh_general_release_dynamic_10.10.2017_dev.fasta", multithread=threads,tryRC=T) ## CHANGE to directory and pertinent database
taxHS.Rs <- assignTaxonomy(seqtab.Rs.nochim, "/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/sh_general_release_dynamic_10.10.2017_dev.fasta", multithread=threads,tryRC=T) ## CHANGE to directory and pertinent database


#taxHS <- assignTaxonomy(seqtab.nochim, "/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/rdp_train_set_16.fa.gz", multithread=threads) ## CHANGE to directory and pertinent database
colnames(taxHS) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")
colnames(taxHS.Fs) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")
colnames(taxHS.Rs) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")


#unname(head(taxHS))
#unname(tail(taxHS))


write.table(taxHS, file = paste0(outdir,"/taxa_SV.tsv"), quote=FALSE)
write.table(taxHS.Fs, file = paste0(outdir,"/taxa_Fs_SV.tsv"), quote=FALSE)
write.table(taxHS.Rs, file = paste0(outdir,"/taxa_Rs_SV.tsv"), quote=FALSE)


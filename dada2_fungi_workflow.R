
### arguments taht need to be provided: folder of FastQ file locations, RunID and output folder location
args <- commandArgs(trailingOnly = TRUE)

#outbase_std="/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10_Fungi"

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


library(dada2)
version=packageVersion("dada2")
set.seed(666)
cat(paste0("Running DADA2 version ",version,"\n"))

fns <- list.files(path)

####################
### load your data
####################

fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) 

fnFs <- fastqs[grepl("R1_001.fastq.gz", fastqs)] 
fnRs <- fastqs[grepl("R2_001.fastq.gz", fastqs)] 

## Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, "-L1_S[0-9]+_L001_R1_001.fastq.gz"), `[`, 1) ## IKMB filenaming scheme
cat(paste0("Starting processing for ",length(sample.names)," samples\n"))

## Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

filt_path <- file.path(outdir, "filtered") # Place filtered files in filtered/subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

## Filter  the forward and reverse reads:
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,150), 
                     trimLeft=c(5, 5),
                     maxN=0, maxEE=c(4,4), truncQ=5, rm.phix=TRUE,
                     compress=TRUE, multithread=threads) 

### remove samples with 0 reads after filtering
exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]
sample.names=sample.names[exists]

#########################
## Learn Errors
#########################
## Learn forward error rates
errF <- learnErrors(filtFs, nbases=numbases, multithread=threads,randomize=T)
## Learn reverse error rates
errR <- learnErrors(filtRs, nbases=numbases, multithread=threads,randomize=T)

saveRDS(errF, paste0(outdir,"/errors/errF.Rds"))
saveRDS(errR, paste0(outdir,"/errors/errR.Rds"))

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


#### this is the crucial step for fungal data. because of the large amplicon variation, some sequences will not overlap
#### to still keep fwd and rev read information, those sequences not overlapping are concatenated
merger <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, returnRejects=TRUE)
concat <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE)
m2=merger
for(x in sample.names){m=merger[[x]]; c=concat[[x]]; m[!m$accept,] <- c[!m$accept,]; m2[[x]] <- m}


############################
## Construct sequence table
############################
seqtab <- makeSequenceTable(m2)

saveRDS(seqtab, paste0(outdir,"/seqtab.Rds"))

###################
## Remove chimeras
###################
## Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads, verbose=TRUE)

saveRDS(seqtab.nochim, paste0(outdir,"/seqtab_nochim.Rds"))

#######################
### Taxonomy assignment 
#######################

tax<-assignTaxonomy(seqtab.nochim, "sh_general_release_dynamic_10.10.2017_dev.fasta", tryRC=T, multithread=threads)
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")
write.table(tax, file = paste0(outdir,"/taxa_SV.tsv"), quote=FALSE)


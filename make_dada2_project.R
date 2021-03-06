#### Start by invoking Rscript make_dada2_project.R [Project name] [sequence ID file] [taxonomic assigment processes] [outputfolder] [basefolder where the original seqtabs are stored]
#### a folder with the name of the project will be created in the specified outputfolder
### tax assignment is at the moment coded in numbers
### 	1: RDP16 + Bayesian; 2: SILVA132 + Bayesian; 3: GTDB86 + Bayesian
### 	4: RDP16 + IDTAX; 5: SILVA132 + IDTAX; 6: GTDB86 + IDTAX
#### the sequence ID file should carry the sample name in the first, the sequencing run id in the second, and the sequencing ID in the third colum, no header should be used!
#### the basefolder has all the sequencing runs that are read from with folders named like the sequencing run IDS in the sequence ID file

args <- commandArgs(trailingOnly = TRUE)

.libPaths("/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/R_Librarires/")
#args=c("AllChimpsTest","DADA2_Projects/AllChimpsTest/seqs2.txt")
threads=16
add_species_bayesian=T
proceed_when_missing=T

project=args[1]

if(is.na(args[3])){use_taxassign=c(1,3,4,6)}else{use_taxassign=as.numeric(strsplit(args[3],split=",")[[1]])}
if(is.na(args[4])){outdirbase="/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/DADA2_Projects/"}else{outdirbase=args[4]}
if(is.na(args[5])){runsbasefolder="/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10"}else{runsbasefolder=args[5]}

outdir=paste0(outdirbase,"/",project)
dir.create(outdir,recursive=T,showWarnings=F)

cat(paste0("\nLoading dependencies...\n"))

suppressMessages(suppressWarnings(require(dada2)))
suppressMessages(suppressWarnings(require(DECIPHER)))

tax_assign_all=expand.grid(c("RDP16","SILVA132","GTDB86"),c("Bayesian","IDTAX"))
colnames(tax_assign_all)=c("Database","Algorithm")
tax_assign=tax_assign_all[use_taxassign,]

fetchfile=read.table(args[2],stringsAsFactors=F,head=F,sep="\t")
colnames(fetchfile)<-c("SampleName","RunID","NewID")

allrunids=unique(fetchfile$RunID)

allrundada=vector("list",length(allrunids))
names(allrundada)<-allrunids


cat(paste0("Fetching a total of ",nrow(fetchfile)," samples from ",length(allrunids), " sequencing runs\n"))
notfound=vector()

for(run in allrunids){
inthisrun<-fetchfile[fetchfile$RunID==run,"NewID"]
cat(paste0("Fetching ",length(inthisrun)," samples from ",run,"..."))
dadathis<-readRDS(paste0(runsbasefolder,"/",run,"/seqtab.Rds"))
dadasub<-dadathis[rownames(dadathis) %in% inthisrun,,drop=F]
dadasub<-dadasub[,colSums(dadathis)>0,drop=F]
if(nrow(dadasub)!=length(inthisrun)){cat(paste0("No data found in this run for ",length(inthisrun) - nrow(dadasub), " samples\n"));notfound<-c(notfound,inthisrun[!inthisrun %in% rownames(dadasub)])}else{cat(paste0("SUCCESSFUL\n"))}
allrundada[[run]]<-dadasub
}


cat(paste0("\nNo data found for ",length(notfound)," samples\n"))
print(fetchfile[fetchfile$NewID %in% notfound,])

if(proceed_when_missing==F){
  cat(paste0("\nABORTING\n"))
  break
}

cat(paste0("\nMerging sequence tables...\n"))


st.all<-allrundada[[1]]
for(mergerun in allrunids[-1]){
st.all <- mergeSequenceTables(st.all,allrundada[[mergerun]])
}

seqtab1_outfile=paste0(outdir,"/seqtab.Rds")

saveRDS(st.all, seqtab1_outfile)

cat(paste0("\nRemoving chimeras...\n"))


seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
#dim(seqtab.nochim)

fracbim=round(100*sum(seqtab.nochim)/sum(st.all),2)

cat(paste0("Of all data ",fracbim,"% were kept as non-chimeric\n"))
seqtab_outfile=paste0(outdir,"/seqtab_nochim.Rds")

cat(paste0("Writing Sequence Table to ",seqtab_outfile,"\n"))

saveRDS(seqtab.nochim, seqtab_outfile)


###################
## Assign taxonomy
###################

cat(paste0("\nRunning sequence classification using ",nrow(tax_assign)," difference approaches\n"))
dir.create(paste0(outdir,"/tax_assign"),recursive=T,showWarnings=F)

bayesian_dbs=matrix(c("RDP16","SILVA132","GTDB86","/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/rdp_train_set_16.fa.gz","/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/silva_nr_v132_train_set.fa.gz","/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/GTDB_bac_ssu_r86_dada2.fasta.gz"),ncol=2)
bayesian_species_dbs=matrix(c("RDP16","SILVA132","/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/rdp_species_assignment_16.fa.gz","/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/silva_species_assignment_v132.fa.gz"),ncol=2)
idtax_dbs=matrix(c("RDP16","SILVA132","GTDB86","/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/RDP_v16-mod_March2018.RData","/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/SILVA_SSU_r132_March2018.RData","/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/GTDB_r86-mod_September2018.RData"),ncol=2)

if("Bayesian" %in% tax_assign$Algorithm){
for(db in tax_assign[tax_assign$Algorithm=="Bayesian","Database"]){
cat(paste0("Assigning taxonomy with Bayesian classifier and ",db," database...\n"))

dbfile=bayesian_dbs[match(db,bayesian_dbs[,1]),2]

taxHS <- assignTaxonomy(seqtab.nochim, dbfile, multithread=threads) ## CHANGE to directory and pertinent database
colnames(taxHS) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

if(add_species_bayesian==T){
 if(db %in% c("RDP16","SILVA132")){
  cat(paste0("Assigning species with Bayesian classifier and ",db," database...\n"))
  species_dbfile=bayesian_species_dbs[match(db,bayesian_species_dbs[,1]),2]
  taxHS <- addSpecies(taxHS, species_dbfile)
}else{cat(paste0("No database for species level assignment available for ",db,"\n"))
}}

write.table(taxHS, file = paste0(outdir,"/tax_assign/taxa_Bayesian_",db,"_SV.tsv"), quote=FALSE)

genclass=round(100*mean(rowSums(seqtab.nochim[,!is.na(taxHS[,6])])/rowSums(seqtab.nochim)),2)
famclass=round(100*mean(rowSums(seqtab.nochim[,!is.na(taxHS[,5])])/rowSums(seqtab.nochim)),2)

cat(paste0("On average per sample, ",famclass,"% of the abundance could be assigned on family and ",genclass,"% on genus level\n"))

}}

######
if("IDTAX" %in% tax_assign$Algorithm){

for(db in tax_assign[tax_assign$Algorithm=="IDTAX","Database"]){
cat(paste0("Assigning taxonomy with IDTAX classifier and ",db," database...\n"))

dbfile=idtax_dbs[match(db,idtax_dbs[,1]),2]

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load(dbfile) # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="both", processors=threads, verbose=FALSE) # use all processors

#### "ranks" was hard coded before, which did not work for all databases
#ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # ranks of interest
#### now the script checks for levels, first occurence and then extracts the name; keeping no "sub" ranks => suborder etc.

allranks <- unique(trainingSet$ranks[sapply(sort(unique(trainingSet$levels)), function(x) min(which(trainingSet$levels==x)))])
ranks<-allranks[as.vector(unlist(sapply(c("^kingdom$|^domain$","^phylum$","^class$","^order$","^family$","^genus$","^species$"), function(x) grep(x, allranks, perl=T, ignore.case=T))))]

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        gsub(" ","_",gsub("\"","",taxa))
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

write.table(taxid, file = paste0(outdir,"/tax_assign/taxa_IDTAX_",db,"_SV.tsv"), quote=FALSE)

genclass=round(100*mean(rowSums(seqtab.nochim[,!is.na(taxid[,6])])/rowSums(seqtab.nochim)),2)
famclass=round(100*mean(rowSums(seqtab.nochim[,!is.na(taxid[,5])])/rowSums(seqtab.nochim)),2)

cat(paste0("On average per sample, ",famclass,"% of the abundance could be assigned on family and ",genclass,"% on genus level\n"))

}}

cat(paste0("\nAll steps finished. All output data can be found in ",outdir,"\n"))


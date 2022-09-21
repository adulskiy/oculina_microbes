# 16S Analysis of Ana's Oculina Data
# Author: Ana Dulskiy
# Date: 9 Aug 2022
### Based on DADA2 Pipeline 1.16 Walkthrough 
### & Nicola Kriefall's Moorea holobiont analysis

#~########################~#
##### PRE-PROCESSING #######
#~########################~#

# Redoing analysis with all data together (both original run + july 22)
# Skipping cutadapt steps as I've already done it with all the data (skip to dada2)
#fastq files should have R1 & R2 designations for PE reads
#Also - some pre-trimming. Retain only PE reads that match amplicon primer. 
# Remove reads containing Illumina sequencing adapters
# I ran cutadapt part (pre-dada2) all on longleaf

#in Terminal home directory:
#following instructions of installing BBtools from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/
#1. download BBMap package, sftp to installation directory
## I followed the instructions at this link: https://astrobiomike.github.io/unix/installing_tools#first-we-need-a-happy-bin
#2. untar: 
#tar -xvzf BBMap_(version).tar.gz
#3. test package:
#cd bbmap
#~/happy_bin/bbmap/stats.sh in=~/bin/bbmap/resources/phix174_ill.ref.fa.gz


# my adaptors for 16S, which I saved as "adaptors.fasta"
# >forward
# AATGATACGGCGACCAC
# >forwardrc
# GTGGTCGCCGTATCATT
# >reverse
# CAAGCAGAAGACGGCATAC
# >reverserc
# GTATGCCGTCTTCTGCTTG


#primers for 16S, which I saved as "primers.fasta": 
# >forward
# GTGYCAGCMGCCGCGGTA
# >reverse
# GGACTACHVGGGTWTCTAAT

### primer check below showed some rev complement primers still in one of my samples, so redoing part of the preprocessing

## Still in terminal - making a sample list based on the first phrase before 
## the underscore in the .fastq name
#ls *R1.fastq | cut -d '_' -f 1 > samples.list

## Get rid of reads that still have the adaptor sequence, shouldn't be there, I didn't have any
# $ for file in $(cat samples.list); do ~/happy_bin/bbmap/bbduk.sh in1=${file}_16S_R1.fastq in2=${file}_16S_R2.fastq ref=adaptors.fasta out1=${file}_R1_NoIll.fastq out2=${file}_R2_NoIll.fastq; done &>bbduk_NoIll.log

## Get rid of first 4 bases (degenerate primers created them)
# $ for file in $(cat samples.list); do ~/happy_bin/bbmap/bbduk.sh in1=${file}_R1_NoIll.fastq in2=${file}_R2_NoIll.fastq ftl=4 out1=${file}_R1_NoIll_No4N.fastq out2=${file}_R2_NoIll_No4N.fastq; done &>bbduk_No4N.log

## Only keep reads that start with the 16S primer
# $ for file in $(cat samples.list); do ~/happy_bin/bbmap/bbduk.sh in1=${file}_R1_NoIll_No4N.fastq in2=${file}_R2_NoIll_No4N.fastq restrictleft=20 k=10 literal=GTGYCAGCMGCCGCGGTAA,GGACTACNVGGGTWTCTAAT copyundefined=t outm1=${file}_R1_NoIll_No4N_16S.fastq outu1=${file}_R1_check.fastq outm2=${file}_R2_NoIll_No4N_16S.fastq outu2=${file}_R2_check.fastq; done &>bbduk_16S.log
## higher k = more reads removed, but can't surpass k=20 or 21

## Install cutadapt if not installed (need to have python installed in order to run this command):
# $ python3 -m pip install --user --upgrade cutadapt

## Use cutadapt to remove primer
# $ for file in $(cat samples.list); do cutadapt -g GTGYCAGCMGCCGCGGTAA -a ATTAGAWACCCBNGTAGTCC -G GGACTACHVGGGTWTCTAAT -A TTACCGCGGCKGCTGRCAC -n 2 --discard-untrimmed -o ${file}_R1.fastq -p ${file}_R2.fastq ${file}_R1_NoIll_No4N_16S.fastq ${file}_R2_NoIll_No4N_16S.fastq; done &>clip.log
# $ for file in $(cat samples.list); do cutadapt -g GTGYCAGCMGCCGCGGTAA -a ATTAGAWACCCBNGTAGTCC -G GGACTACHVGGGTWTCTAAT -A TTACCGCGGCKGCTGRCAC -n 2 --discard-untrimmed -o ${file}_R1.fastq -p ${file}_R2.fastq ${file}_16S_R1.fastq ${file}_16S_R2.fastq; done &>clip.log

##-g regular 5' forward primer 
##-G regular 5' reverse primer
##-o forward out
##-p reverse out
##-max-n 0 means 0 Ns allowed
##this overwrote my original renamed files 

# did sftp of *_R1.fastq & *_R2.fastq files to the folder to be used in dada2


#~########################~#
##### DADA2 BEGINS #########
#~########################~#

#installing/loading packages:
library(dada2); packageVersion("dada2")
#Version 1.22.0

library(ShortRead); packageVersion("ShortRead")
#1.52.0

library(Biostrings); packageVersion("Biostrings")
#2.62.0


path <- "~/oculina_old/16S_preprocessed" # CHANGE ME to the directory containing the fastq files after unzipping.

fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))

fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))

#### check for primers ####
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME to your reverse primer sequence

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
REV.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress=TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[4]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[4]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[4]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[4]]))

#no primers -- skip cutadapt part 2

### cutadapt part 2 ###
## Install cutadapt to local computer and set path here (python3 -m pip install --user --upgrade cutadapt)
#to find where cutadapt is installed do $which cutadapt and then set path to that:
#$export PATH=/nas/longleaf/home/adulskiy/.local/bin/cutadapt/$PATH
cutadapt <- "/nas/longleaf/home/adulskiy/.local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

## look for cutadapt path and create one if it doesn't exist
path_cut <- file.path(path, "cutadapt")
if(!dir.exists(path_cut)) dir.create(path_cut)

## create list of forward/reverse samples to put the cutadapt samples
fnFs_cut <- file.path(path_cut, basename(fnFs)) # 16S path with the subdirectory 'cutadapt' then the sample names
fnRs_cut <- file.path(path_cut, basename(fnRs))

# forward/reverse complements of primers
FWD_RC <- dada2:::rc(FWD)
REV_RC <- dada2:::rc(REV)

R1_flags <- paste("-g", FWD, "-a", REV_RC) # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R2_flags <- paste("-G", REV, "-A", FWD_RC) # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1_flags, R2_flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs_cut[i], "-p", fnRs_cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
#Overview of removed sequences
#length	count	expect	max.err	error counts
#3	1	5.4	0	1
#59	1	0.0	1	0 1
#216	2	0.0	1	2
#219	1	0.0	1	0 1
#221	1	0.0	1	1
#227	3	0.0	1	0 3

#towards the beginning
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs_cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs_cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs_cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs_cut[[1]]))

#at the end
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs_cut[[4]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs_cut[[4]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs_cut[[4]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs_cut[[4]]))




#### Visualizing raw data ####


#First, lets look at quality profile of R1 reads
plotQualityProfile(fnFs.filtN[c(1,2,3,4)])

#Then look at quality profile of R2 reads
plotQualityProfile(fnRs.filtN[c(1,2,3,4)])

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#changing a bit from default settings - maxEE=1 (1 max expected error, more conservative), truncating length at 200 bp for both forward & reverse [leaves ~50bp overlap], added "trimleft" to cut off primers [18 for forward, 20 for reverse]
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(175,175), #leaves ~50bp overlap
                     maxN=0, #DADA does not allow Ns
                     maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     #trimLeft=c(18,20), #N nucleotides to remove from the start of each read
                     rm.phix=TRUE, #remove reads matching phiX genome
                     matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)
tail(out)


#~############################~#
##### Learn Error Rates ########
#~############################~#

#setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE) 
plotErrors(errR, nominalQ=TRUE) 



#~############################~#
##### Dereplicate reads ########
#~############################~#

#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#~###############################~#
##### Infer Sequence Variants #####
#~###############################~#

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 

dadaFs[[1]]
dadaRs[[1]]


#~############################~#
##### Merge paired reads #######
#~############################~#

#To further cull spurious sequence variants
#Merge the denoised forward and reverse reads
#Paired reads that do not exactly overlap are removed

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

summary((mergers[[1]]))

#We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs.

#~##################################~#
##### Construct sequence table #######
#~##################################~#
#a higher-resolution version of the “OTU table” produced by classical methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

plot(table(nchar(getSequences(seqtab)))) #real variants appear to be right in that 244-264 window


#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(240,260)] #again, being fairly conservative wrt length


#~############################~#
##### Remove chimeras ##########
#~############################~#
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#Identified 798 bimeras out of 8697 input sequences. 

sum(seqtab.nochim)/sum(seqtab2)
#0.9749089
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. 


#~############################~#
##### Track Read Stats #########
#~############################~#

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

#setwd("~/oculina_old/data")
write.csv(track,file="oculina_rev_readstats.csv",row.names=TRUE,quote=FALSE)


#~############################~#
##### Assign Taxonomy ##########
#~############################~#

# #Using package DECIPHER as an alternative to 'assignTaxonomy'
# install.packages("BiocManager")
# BiocManager::install("DECIPHER")
# library(DECIPHER); packageVersion("DECIPHER")
# version 2.22.0
# #citation("DECIPHER")
# 
# #http://DECIPHER.codes/Downloads.html. Download the SILVA SSU r132 (modified) file to follow along.
# 

### not doing this method (skip to the Assign Taxonomy section)
# dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
# load("~/Downloads/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
#ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE, threshold=50) # use all processors
# ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
# taxid <- t(sapply(ids, function(x) {
#   m <- match(ranks, x$rank)
#   taxa <- x$taxon[m]
#   taxa[startsWith(taxa, "unclassified_")] <- NA
#   taxa
# }))
# colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)


#also doing other taxonomy method:
#Assign Taxonomy


taxa <- assignTaxonomy(seqtab.nochim, "~/oculina_old/setup/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

unname(head(taxa))
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


taxa.plus <- addSpecies(taxa, "~/oculina_old/setup/tax/silva_species_assignment_v132.fa.gz",tryRC=TRUE,verbose=TRUE)

saveRDS(taxa.plus, file="oculina16s_rev_taxaplus.rds") 
#86 out of 7899 were assigned to the species level.
#Of which 77 had genera consistent with the input table

saveRDS(taxa, file="oculina16s_rev_taxa.rds")
#write.csv(taxa.plus, file="oculina16s_rev_taxaplus.csv")
#write.csv(taxa, file="oculina16s_rev_taxa.csv")

saveRDS(seqtab.nochim, file="oculina16s_rev_seqtab.nochim.rds")
write.csv(seqtab.nochim, file="oculina16s_rev_seqtab.nochim.csv")

#### Read in previously saved datafiles ####
#setwd(~/oculina_old/data)
seqtab.nochim <- readRDS("oculina16s_rev_seqtab.nochim.rds")
taxa <- readRDS("oculina16s_rev_taxa.rds")
taxa.plus <- readRDS("oculina16s_rev_taxaplus.rds")


#~############################~#
##### handoff 2 phyloseq #######
#~############################~#

#BiocManager::install("phyloseq")
library('phyloseq')
library('ggplot2')
library('Rmisc')
library(cowplot)
library(ShortRead)

#import dataframe holding sample information
samdf<-read.csv("oculina16s_sampledata_symdens.csv")
head(samdf)
rownames(samdf) <- samdf$id

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa.plus))

ps # 7899 taxa, 77 samples

#### first look at data ####
ps_glom <- tax_glom(ps, "Family")
plot_bar(ps_glom, x="site", fill="Family")+
  theme(legend.position="none")

#phyloseq object with shorter names - doing this one instead of one above
ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))
#making output fasta file 
path='~/oculina_old/data/oculina16s_rev.fasta'
uniquesToFasta(seqtab.nochim, path, ids = ids, mode = "w", width = 20000)


colnames(seqtab.nochim)<-ids
taxa2 <- cbind(taxa.plus, rownames(taxa.plus)) #retaining raw sequence info before renaming
rownames(taxa2)<-ids

#phyloseq object with new taxa ids
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa2))
ps #7899 taxa, 77 samples


#### remove mitochondria, chloroplasts, non-bacteria #### 
ps.mito <- subset_taxa(ps, (Family=="Mitochondria"))
ps.mito #97 taxa to remove
ps.chlor <- subset_taxa(ps, (Order=="Chloroplast"))
ps.chlor #414 taxa to remove
ps.notbact <- subset_taxa(ps, (Kingdom!="Bacteria") | is.na(Kingdom))
ps.notbact #242 taxa to remove


ps.nomito <- subset_taxa(ps, (Family!="Mitochondria") | is.na(Family))
ps.nomito #7802 taxa
ps.nochlor <- subset_taxa(ps.nomito, (Order!="Chloroplast") | is.na(Order))
ps.nochlor #7388 taxa
ps.clean <- subset_taxa(ps.nochlor, (Kingdom=="Bacteria"))
ps.clean #7146 taxa

#just archaea
ps.arch <- subset_taxa(ps.nomito, (Kingdom=="Archaea"))
ps.arch #164 taxa

#### identifying contamination ####
#BiocManager::install("decontam")
library(decontam)

df <- as.data.frame(sample_data(ps.clean)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps.clean)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=site)) + 
  geom_point()

ggplot(data=df, aes(x=Index, y=LibrarySize, color=season, group = site)) + 
  geom_point(aes(shape = site))


sample_data(ps.clean)$is.neg <- sample_data(ps.clean)$site == "neg"
contamdf.prev <- isContaminant(ps.clean, neg="is.neg",threshold=0.5)
table(contamdf.prev$contaminant)

# FALSE  TRUE 
# 7140    6

# Make phyloseq object of presence-absence in negative controls and true samples

ps.pa <- transform_sample_counts(ps.clean, function(abund) 1*(abund>0))
#ps.pa.neg <- prune_samples(sample_data(ps.pa)$site == "neg", ps.pa) #none
ps.pa.pos <- prune_samples(sample_data(ps.pa)$site != "neg", ps.pa)
# Make data.frame of prevalence in positive and negative samples
#df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa),
#                    contaminant=contamdf.prev$contaminant)
#ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
#  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


#remove from ps.clean:
ps.clean1 <- prune_taxa(!contamdf.prev$contaminant,ps.clean)
#also remove negative controls, don't need them anymore I think
ps.cleaner <- subset_samples(ps.clean1,(site!="neg"))
#7140 taxa, 74 samples
#saveRDS(ps.cleaner, file ="ps.cleaner_rev.RDS")


#### blast asvs to NCBI to see if any eukaryotes got through ####
##Running blast on longleaf to make organism match files for my 16s data
##used 'oculina16s_rev.fasta' made way above

## On longleaf:
##submitted the following job:
#nano blast_taxid.sh
## 
#!/bin/bash
#SBATCH -p general
#SBATCH -t 2-
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8gb
#SBATCH --output <oculina16s_rev_taxids.out>

#module load blast
#blastn -query oculina16s.fasta -db nt -outfmt "6 std staxids sskingdoms" -evalue 1e-5 -max_target_seqs 5 -out oculina16s_taxids.out -remote
#[exit]
#qsub -pe omp 28 blast_taxid.sh

##takes a very long time (I had ~7300 ASVs for my larger dataset, took like 32 hours)


##now getting taxonomy info:

##download/install taxonkit things, more instructions here:
##https://bioinf.shenwei.me/taxonkit/usage/

# cd /oculina/tax
# conda install -c bioconda taxonkit
# #command taxonkit should work now

##extracting taxa ids from blast output for taxonkit (in terminal):
# awk -F " " '{print $13}' oculina16s_rev_taxids.out > ids_rev
# taxonkit lineage ids_rev > ids_rev.tax
# cut -f1 oculina16s_rev_taxids.out > ids_rev.seq; paste ids_rev.seq ids_rev.tax > ids_rev.seq.tax
# grep "Eukaryota" ids_rev.seq.tax | cut -f1 | sort | uniq > euk.contam.asvs_rev

# convert euk.contam.asvs_rev to .csv file
eukasvs.rev <- read_table("euk.contam.asvs_rev", col_names = FALSE)
write_csv(eukasvs.rev, file = "euk.contam.asvs_rev.csv", col_names = FALSE)

##transferring euk.contam.asvs_rev to back here
##remove from ps.cleaner
##should be 151 to remove
euks <- read.csv("euk.contam.asvs_rev",header=FALSE)
euks_names <- euks$V1
alltaxa <- taxa_names(ps.cleaner) #should be 6582
`%!in%` <- Negate("%in%")
keepers <- alltaxa[(alltaxa %!in% euks_names)] #keepers = 6573, so that means only 9 euks got through above
ps.cleanest <- prune_taxa(keepers, ps.cleaner) 
str(ps.cleanest)
#7005 in ps.cleanest

seqtab.cleanest <- data.frame(otu_table(ps.cleanest))
#write.csv(seqtab.cleanest,file="oculina16s_rev_seqtab.cleanest.csv")
seqtab.cleanest <- read.csv("oculina16s_rev_seqtab.cleanest.csv",row.names=1)

##re-read in cleaned phyloseq object
saveRDS(ps.cleanest,file="ps.cleanest_rev.rds")


#### rarefy #####
library(vegan)

seqtab.cleanest <- data.frame(ps.cleanest@otu_table)
samdf.cleanest <- data.frame(ps.cleanest@sam_data)

rarecurve(seqtab.cleanest,step=100,label=FALSE) #after removing contaminants

total <- rowSums(seqtab.cleanest)
#make some plots to figure out where to rarefy
total.df <- data.frame(total)
ggplot(total.df, aes(x = total)) +
  geom_histogram()

ggplot(total.df, aes(x=1, y=total))+
  geom_boxplot()

ggplot(total.df, aes(x=1, y=total))+
  geom_jitter()+
  scale_y_log10()

total.df %>%
  arrange(total) %>%
  ggplot(aes(x=1:nrow(.),y=total))+
  geom_line()  #shows # of samples on x and # of seqs in each sample on y

total.df %>%
  arrange(total) %>%
  print(20) 
subset(total, total <600) #3 samples need to be removed (D6, G9,)

#originally got rid of all samples with less than 2000 seqs, changed to 600 on 6/21:
#subset(total, total <2000)
#6 samples
# identified by MCMC.OTU below as being too low 

row.names.remove <- c("D6","H2","G9")
seqtab.less <- seqtab.cleanest[!(row.names(seqtab.cleanest) %in% row.names.remove),]
samdf.less <- samdf.cleanest[!(row.names(samdf.cleanest) %in% row.names.remove), ]
str(samdf.less) #71 samples left
ps.less <- phyloseq(otu_table(seqtab.less, taxa_are_rows=FALSE), 
                    sample_data(samdf.less), 
                    tax_table(taxa2))
###saving
#setwd("~/oculina/data")
#write.csv(seqtab.less, file="oculina16s_rev_seqtab.less")
#save(ps.less,file="ps.less_rev.Rdata")


#for rarefied data, remove everything < 2000
subset(total, total <2000) # 7 samples
row.names.remove <- c("D6","F3","H2","F7","H4","H7", "G9")
seqtab.less <- seqtab.cleanest[!(row.names(seqtab.cleanest) %in% row.names.remove),]
samdf.rare <- samdf.cleanest[!(row.names(samdf.cleanest) %in% row.names.remove), ]
seqtab.rare <- rrarefy(seqtab.less,sample=2000)
rarecurve(seqtab.rare,step=100,label=FALSE)
#phyloseq object but rarefied
ps.rare <- phyloseq(otu_table(seqtab.rare, taxa_are_rows=FALSE), 
                    sample_data(samdf.rare), 
                    tax_table(taxa2))
ps.rare #7005 taxa, 67 samples

#removing missing taxa - lost after rarefying
ps.rare <- prune_taxa(taxa_sums(ps.rare) > 0, ps.rare)
ps.rare #5237 taxa, 67 samples

seqtab.rare <- data.frame(otu_table(ps.rare))

#saving
#### data files - rarefied, decontaminated ####
#saving
#setwd("~/oculina/data")
#write.csv(seqtab.rare, file="oculina16s_rev_seqtab.rare.2k")
#save(taxa2,file="taxa2_rev.Rdata")

seqtab.rare <- read.csv("oculina16s_rev_seqtab.rare.2k",row.names=1)
load("taxa2_rev.Rdata")

ps.rare <- phyloseq(otu_table(seqtab.rare, taxa_are_rows=FALSE), 
                    sample_data(samdf.rare), 
                    tax_table(taxa2)) 
ps.rare #5237 taxa; 67 samples

#### trim underrepresented otus ####
#install.packages("MCMC.OTU")
library(MCMC.OTU)

#formatting the table for mcmc.otu - requires one first column that's 1 through whatever
#& has "X" as column name
nums <- 1:nrow(seqtab.less)  #originally used seqtab.cleanest but I think seqtab.less better
samples <- rownames(seqtab.less)

int <- cbind(sample = 0, seqtab.less)
seq.formcmc <- cbind(X = 0, int)

seq.formcmc$X <- nums
seq.formcmc$sample <- samples

#change second columns value to equal total variables
str(seq.formcmc)
seq.trim.allinfo <- purgeOutliers(seq.formcmc,count.columns=3:7007,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.02)
#no samples with counts below z-score
#1170 ASVs pass filters using seqtab.less

#remove sample info
str(seq.trim.allinfo)
seq.trim <- seq.trim.allinfo[,3:825] # number variables of seq.trim.allinfo, i think

#write.csv(seq.trim,file="oculina16s_rev_seqtab.trim.csv")
seq.trim <- read.csv("oculina16s_rev_seqtab.trim.csv",row.names=1)

#remake phyloseq objects
ps.trim <- phyloseq(otu_table(seq.trim, taxa_are_rows=FALSE), 
                    sample_data(samdf.less), 
                    tax_table(taxa2))
ps.trim #823 taxa, 67 samples

#save(ps.trim,file="ps.trim_rev.Rdata")


#### rarefy - trimmed #####
library(vegan)

seqtab.trim <- data.frame(ps.trim@otu_table)
samdf.trim <- data.frame(ps.trim@sam_data)

rarecurve(seqtab.trim,step=100,label=FALSE) 

total <- rowSums(seqtab.trim)
subset(total, total <2000)
# 3 samples

row.names.remove <- c("E8","E9","H3")
seqtab.less <- seqtab.trim[!(row.names(seqtab.trim) %in% row.names.remove),]

seqtab.trim.rare <- rrarefy(seqtab.less,sample=2000)
rarecurve(seqtab.trim.rare,step=100,label=FALSE)

#phyloseq object but rarefied & trimmed
ps.trim.rare <- phyloseq(otu_table(seqtab.trim.rare, taxa_are_rows=FALSE), 
                         sample_data(samdf.rare), 
                         tax_table(taxa2))
ps.trim.rare #823 taxa, 64 samples

#saving
#### data files - rarefied, decontaminated, trimmed ####
#saving
#write.csv(seqtab.trim.rare, file="oculina16s_rev_seqtab.trim.rare.2k.csv")
#saveRDS(ps.trim.rare,file="ps.trim.rare_rev.Rdata")

#### making fasta file for picrust2 - trimmed not rarefied ####
library(phyloseq)
library(dada2)

#if needed:
#setwd("~/oculina/data")
seqtab.trim <- read.csv(file="oculina16s_rev_seqtab.trim.rare.2k.csv",row.names=1)
load("taxa2.Rdata")
samdf <- read.csv(file="oculina16s_sampledata_symdens.csv")
row.names(samdf) <- samdf$id

ps.trim <- phyloseq(otu_table(seqtab.trim, taxa_are_rows=FALSE), 
                    sample_data(samdf), 
                    tax_table(taxa2))
ps.trim #796 taxa during revision

trim.otu <- as.matrix(ps.trim@otu_table)
trim.taxa <- data.frame(ps.trim@tax_table)
rownames(trim.taxa)==colnames(trim.otu)

colnames(trim.otu) <- trim.taxa$V8
ids <- rownames(trim.taxa)

path="~/oculina/data/oculina16s_rev_trimmed.fasta"
uniquesToFasta(trim.otu, path, ids = ids, mode = "w", width = 20000)

#setwd("~/oculina_old/setup/tax")
#re-formatting seq table so picrust likes it:
#a tab-delimited table with ASV ids as the first column and sample abundances as all subsequent columns
seqtab.trim.t <- t(seqtab.trim)
write.table(seqtab.trim.t,file="oculina16s_rev_seqtab.trim.t.txt")
#manually removed the quotation marks that appeared in the file, and converted to tab delimited file from Excel

#### moving on to oculina16s_revised_analysis.R script in other folder ####



# check read depths
load("~/oculina/data/ps.less_rev.Rdata")
mean(sample_sums(ps.less))
#15488.52
sd(sample_sums(ps.less))
#16011.4
sum(sample_sums(ps.less))
#1,099,685
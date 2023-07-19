

# getting ready
library(dada2)
packageVersion("dada2")
# 1.22.0
library(ShortRead)
packageVersion("ShortRead")
# 1.52.0
library(Biostrings)
packageVersion("Biostrings")
# 2.62.0

args = commandArgs(trailingOnly=TRUE)
input_table = args[1]
data_dir = args[2]
odir = args[3]
indata <- read.table(input_table,sep='\t',header=1,)
#indata <- read.table("/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/results/METcnkiM/input.tab",sep='\t',header=1,)
#odir = '/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/results/METcnkiM/Pdada2_output'
filtFs <- file.path(data_dir,paste(indata$sample.ID,'_R1.clean.fq.gz',sep=''))
filtRs <- file.path(data_dir,paste(indata$sample.ID,'_R2.clean.fq.gz',sep=''))

sample.names <- as.array(indata$sample.ID)
# learn errors
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# for (i in 1:length(derepFs)){
# m<-as.matrix(derepFs[[i]])
# n <- names(derepFs)[i]
# dir.create(file.path(odir,'derep'))
# write.csv(m,file.path(odir,'derep',paste0(n,'.csv')), row.names = FALSE)
# }

# name derep-class objects by the sample names
names(derepFs) <- indata$sample.ID
names(derepRs) <- indata$sample.ID

# sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# for (i in 1:length(mergers)){
# mergers1<-as.matrix(mergers[[i]])
# n <- names(mergers)[i]
# #print(n)
# dir.create(file.path(odir,'merged'))
# write.csv(mergers1,file.path(odir,'merged',paste0(n,'.csv')), row.names = FALSE)
# }

# construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# [1]   48 2068

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# Identified 964 bimeras out of 2068 input sequences.

# seq length distribution
table(nchar(getSequences(seqtab.nochim)))
# As expected, quite a bit of length variability in the amplified nirS region.
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
# track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(
    sapply(derepFs, getN),
    sapply(dadaFs, getN),  
    sapply(mergers,getN), 
    rowSums(seqtab.nochim))


# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("dereplicateF", "denoisedF", "merged","nonchim")
rownames(track) <- sample.names
# kept the majority of our raw reads, run successfully

### ### ### ### ### ###
###  Fasta output   ###
### ### ### ### ### ###
uniquesToFasta(seqtab.nochim, file.path(odir,"rep.fasta"))
### ### ### ### ### ###
### Export table   ####
### ### ### ### ### ###

samples.out <- rownames(seqtab.nochim)
asv.1=as.data.frame(t(seqtab.nochim))
asv.1$sum=rowSums(asv.1)
asv.final=asv.1[order(-asv.1$sum),]
asv.final$Sequence=rownames(asv.final)
rownames(asv.final) = sprintf("preASV%06d", 1:nrow(asv.final))
write.table(asv.final, file.path(odir,"ASV_table.tsv"),  quote = FALSE,sep='\t')
write.table(track, file.path(odir,"ASV_stats.txt"), quote = FALSE,sep='\t')
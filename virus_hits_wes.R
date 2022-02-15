###################################################
## sequencing statistics
## for viruses
###################################################

library(reshape2)
library(ggplot2)


# system('mkdir plots');
dirElem  = strsplit(getwd(),'/')[[1]]
expname  = paste('_',dirElem[length(dirElem)],"_wes",sep='');

## virus statistics
files = list.files(path="bam_wes",pattern=".u.idxstats$", full.names=T)
samples = gsub("^bam_wes/","",files)
samples = gsub(".u.idxstats","",samples)
labs <- gsub("bam_wes/", "", gsub(".u.idxstats$", "", files, perl=TRUE))

i = 1
tmp = read.csv(files[i], sep="\t",as.is=T,header=F)[,-4]
tmp = tmp[c(grep("^[NAF]",tmp$V1),which(tmp$V1 == "chrEBV")),] # get viruses only
stats = tmp
names(stats) = c("virus", "genome_length", samples[i])
if (length(labs)>1) {
  for (i in 2:length(labs)) {
	tmp = read.csv(files[i], sep="\t",as.is=T,header=F)[,-4]
	tmp = tmp[c(grep("^[NAF]",tmp$V1),which(tmp$V1 == "chrEBV")),c(1,3)]
	names(tmp) = c("virus", samples[i])
	stats = merge(stats, tmp, by="virus")
  }
}

# EBV         HBV         HCV          HHV-8       HIV-1       HIV-2       HPV         HTLV-1      HTLV-2      MLV        SMRV  XMRV
# NC_007605.1 NC_003977.2 NC_004102.1  NC_009333.1 NC_001802.1 NC_001722.1 NC_004500.1 NC_001436.1 NC_001488.1 AF221065.1 NC_001514.1 FN692043.2
# stats$virus
#  [1] "AF221065.1"  "chrEBV"      "FN692043.2"  "NC_001436.1" "NC_001488.1"
#  [6] "NC_001514.1" "NC_001722.1" "NC_001802.1" "NC_003977.2" "NC_004102.1"
# [11] "NC_004500.1" "NC_009333.1"
genome_length = stats$genome_length
stats = stats[c(2,9,10,12,8,7, 11,4,5,1,6,3),-c(1:2)]
rownames(stats) = c("EBV", "HBV", "HCV", "HHV-8", "HIV-1", "HIV-2", "HPV", "HTLV-1", "HTLV-2", "MLV", "SMRV", "XMRV")



## get total reads of merged sequences
files = list.files(path="star_wes",pattern="join.Log.final.out$", full.names=T)
cols = c("TotalReads", "UniqueHits")
i = 1
stats2 = read.csv(files[i], sep="\t",as.is=T)
if (length(files)>1) {
  for (i in 2:length(files)) {
	  stats2 = merge(stats2, read.csv(files[i], sep="\t",as.is=T), by="Started.job.on..")
  }
}
colnames(stats2) = c("category", labs)
stats2 = stats2[c(grep("input reads", stats2$category), grep("Uniquely mapped reads %", stats2$category)),]
rownames(stats2) = cols
## get total reads of unmerged, separated seqs and sum up with merged sequences
files = list.files(path="star_wes",pattern="un.Log.final.out$", full.names=T)
for (i in 1:length(files)) {
	stats2[1,1+i] = as.numeric(stats2[1,1+i]) + as.numeric(read.csv(files[i], sep="\t",as.is=T)[7,2])
}


stats = rbind(stats, stats2[1,-1])
write.table(cbind(names(stats),t(stats)), paste('tables/virus_hits',expname,'.csv',sep=''),row.names=F,quote=F,sep='\t',na="");



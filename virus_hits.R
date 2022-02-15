###################################################
## sequencing statistics
## for viruses
###################################################

library(reshape2)
library(ggplot2)


# system('mkdir plots');
dirElem  = strsplit(getwd(),'/')[[1]]
expname  = paste('_',dirElem[length(dirElem)],sep='');

## virus statistics
files = list.files(path="bam",pattern=".u.idxstats$", full.names=T)
samples = gsub("^bam/","",files)
samples = gsub(".u.idxstats","",samples)
labs <- gsub("bam/", "", gsub(".u.idxstats$", "", files, perl=TRUE))

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

tab = cbind(data.frame(t(stats)),Cell=colnames(stats))
names(tab) = c(rownames(stats),"Cell")
tab = melt(tab,id.vars="Cell")

pdf(paste("plots/sum_all_reads_samples", expname, ".pdf",sep=""))
ggplot(tab,aes(variable,value,fill=Cell)) + 
	theme_bw(base_size=16, base_family = "") + ylab("Sum of all reads and samples for one virus") + xlab("") +
	geom_bar(stat="identity") + scale_colour_gradientn(colours=rainbow(4)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none")
dev.off()




## get total reads of merged sequences
files = list.files(path="star",pattern="join.Log.final.out$", full.names=T)
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
files = list.files(path="star",pattern="un.Log.final.out$", full.names=T)
for (i in 1:length(files)) {
	stats2[1,1+i] = as.numeric(stats2[1,1+i]) + as.numeric(read.csv(files[i], sep="\t",as.is=T)[7,2])
}

stats3 = rbind(stats, as.numeric(stats2[-2,-1]))
for (i in 1:ncol(stats3)) {
	stats3[,i] = stats3[,i]/stats3[13,i]*100
}
tab = cbind(data.frame(t(stats3[-nrow(stats3),])),Cell=colnames(stats3))
names(tab) = c(rownames(stats),"Cell")
tab = melt(tab,id.vars="Cell")

pdf(paste("plots/percent_reads_samples", expname, ".pdf",sep=""))
ggplot(tab,aes(variable,value,fill=Cell)) + 
	theme_bw(base_size=16, base_family = "") + ylab("Sum % reads per virus") + xlab("") +
	geom_bar(stat="identity") + scale_colour_gradientn(colours=rainbow(4)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none")
dev.off()
W
# ebv only
tab = cbind(data.frame(t(stats3[1,stats3[1,]>0.01])),Cell=names(stats3[1,stats3[1,]>0.01]))
pdf(paste("plots/ebv_percent", expname, ".pdf",sep=""))
ggplot(tab,aes(Cell,EBV)) + 
	theme_bw(base_size = 16, base_family = "") + ylab("% EBV reads to total reads") + xlab("") +
	geom_bar(stat="identity", position="dodge") +  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

# mlv only
tab = cbind(data.frame(t(stats3[10,stats3[10,]>0.01])),Cell=names(stats3[10,stats3[10,]>0.01]))
pdf(paste("plots/mlv_percent", expname, ".pdf",sep=""))
ggplot(tab,aes(Cell,MLV)) + 
	theme_bw(base_size = 16, base_family = "") + ylab("% MLV reads to total reads") + xlab("") +
	geom_bar(stat="identity", position="dodge") +  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

# mlv+xmrv only
tab = cbind(data.frame(t(stats3[c(10,12),stats3[10,]>0.01])),Cell=names(stats3[c(10,12),stats3[10,]>0.01]))
tab = melt(tab,id.vars="Cell")
pdf(paste("plots/mlv_xmrv_percent", expname, ".pdf",sep=""))
ggplot(tab,aes(Cell,value,fill=variable)) + 
	theme_bw(base_size = 16, base_family = "") + ylab("% MLV reads to total reads") + xlab("") +
	geom_bar(stat="identity", position="dodge") + scale_fill_brewer(palette="Set1") +  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()




stats = rbind(stats, stats2[1,-1])
write.table(cbind(names(stats),t(stats)), paste('tables/virus_hits',expname,'.csv',sep=''),row.names=F,quote=F,sep='\t',na="");
write.table(cbind(names(stats3[-13,]),t(stats3[-13,])), paste('tables/virus_hits_percent',expname,'.csv',sep=''),row.names=F,quote=F,sep='\t',na="");


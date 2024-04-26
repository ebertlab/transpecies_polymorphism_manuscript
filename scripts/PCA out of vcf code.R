
setwd("~/Desktop/new_vcf/minGQ30/")
library(gdsfmt)
library(SNPRelate)

#Create GDS file
snpgdsVCF2GDS("all_186_start_pca_minGQ30_1k.recode.vcf", "1000BP.gds")

# Summary
snpgdsSummary("1000BP.gds")

# Open the GDS file
genofile1000BP <- snpgdsOpen("1000BP.gds")

#PCA
pca <- snpgdsPCA(genofile1000BP, snp.id = NULL, num.thread=2, autosome.only = FALSE)

pc.percent<-pca$varprop*100
head(round(pc.percent,2))

# make a data.frame
tab1000BP <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  EV6 = pca$eigenvect[,6],
                  stringsAsFactors = FALSE)
head(tab1000BP)

write.table(tab1000BP, "tab_1000bp_all.txt",sep = "\t",quote = F,row.names = F,col.names = F)
#PCA Plots
plot(tab1000BP$EV1, tab1000BP$EV2, xlab=paste("PC 1 (",round(pc.percent[1],1),"%)",sep = ""), ylab=paste("PC 2 (",round(pc.percent[2],1),"%)",sep = ""),cex=1.5)

b<-read.table("tab_1000bp_all.txt",header = T)
plot(b$EV1, b$EV2, xlab=paste("PC 1 (",round(pc.percent[1],1),"%)",sep = ""), ylab=paste("PC 2 (",round(pc.percent[2],1),"%)",sep = ""),cex=2, pch= c(21,24,23)[as.numeric(b$species)])

showfile.gds(closeall=TRUE)




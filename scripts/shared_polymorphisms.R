
setwd("~/Desktop/new_TSP/shared_SNPs_aug/polymorphisms/")

## magna

WE<-read.table("magna_WE.tab",header = T,comment.char = "")
EA<-read.table("magna_EA.tab",header = T,comment.char = "")
US<-read.table("magna_NA.tab",header = T,comment.char = "")

colnames(WE)[1]<-c("CHROM")
colnames(EA)[1]<-c("CHROM")
colnames(US)[1]<-c("CHROM")
WE1<-WE[,c("CHROM","POS")]
EA1<-EA[,c("CHROM","POS")]
US1<-US[,c("CHROM","POS")]

## similis


sim_all<-read.table("similis_all.tab",header = T,comment.char = "")

colnames(sim_all)[1]<-c("CHROM")
sim_all1<-sim_all[,c("CHROM","POS")]


## sinensis

sin<-read.table("sinensis.tab",header = T,comment.char = "")
colnames(sin)[1]<-c("CHROM")
sin1<-sin[,c("CHROM","POS")]

## overlap NOT considering the two similis groups

eurasia_magna<-rbind(WE1,EA1)
shared_magna_eurasia<-eurasia_magna[duplicated(eurasia_magna), , drop=F]

magna_all<-rbind(shared_magna_eurasia,US1)
shared_magna_all<-magna_all[duplicated(magna_all), , drop=F]  ## shared SNPs in the three magna groups (30126)

mag_sim<-rbind(shared_magna_all,sim_all1)
shared_mag_sim<-mag_sim[duplicated(mag_sim), , drop=F]  ## shared SNPs between magna and similis (522)

mag_sim_sin<-rbind(shared_mag_sim,sin1)
shared_mag_sim_sin<-mag_sim_sin[duplicated(mag_sim_sin), , drop=F]  ## shared SNPs in between magna and similis and sinensis (131)

write.table(shared_mag_sim_sin,file = "shared_polymorphisms_sim_all_md05_dp10_gq30_maf005.txt",quote = F,row.names = F,col.names = F,sep = "\t")

## overlap NOT considering the two similis groups and the NA magna clones

eurasia_magna<-rbind(WE1,EA1)
shared_magna_eurasia<-eurasia_magna[duplicated(eurasia_magna), , drop=F]

mag_sim<-rbind(shared_magna_eurasia,sim_all1)
shared_mag_sim<-mag_sim[duplicated(mag_sim), , drop=F]  ## shared SNPs between magna (no NA) and similis (2020)

mag_sim_sin<-rbind(shared_mag_sim,sin1)
shared_mag_sim_sin<-mag_sim_sin[duplicated(mag_sim_sin), , drop=F]  ## shared SNPs in between magna and similis and sinensis (307)

write.table(shared_mag_sim_sin,file = "shared_polymorphisms_all_no_sim_no_NA_md05_dp10_gq30_maf01.txt",quote = F,row.names = F,col.names = F,sep = "\t")


## overlap magna all - sinensis

eurasia_magna<-rbind(WE1,EA1)
shared_magna_eurasia<-eurasia_magna[duplicated(eurasia_magna), , drop=F]

magna_all<-rbind(shared_magna_eurasia,US1)
shared_magna_all<-magna_all[duplicated(magna_all), , drop=F]  ## shared SNPs in the three magna groups (30126)

mag_sin<-rbind(shared_magna_all,sin1)
shared_mag_sin<-mag_sin[duplicated(mag_sin), , drop=F]  ## shared SNPs in between magna and similis (140)

write.table(shared_mag_sin,file = "shared_polymorphisms_magna_all_sin_md05_dp10_gq30_maf01.txt",quote = F,row.names = F,col.names = F,sep = "\t")





## overlap similis - sinensis


all_similis<-rbind(sim_rus1,sim_isr1)
shared_all_similis<-all_similis[duplicated(all_similis), , drop=F]  ## shared SNPs in the two similis groups (74458)

sim_sin<-rbind(shared_all_similis,sin1)
shared_sim_sin<-sim_sin[duplicated(sim_sin), , drop=F]  ## shared SNPs in between sinensis and similis (6241)

write.table(shared_sim_sin,file = "shared_polymorphisms_sim_all_sin_md05_dp10_gq30_maf01.txt",quote = F,row.names = F,col.names = F,sep = "\t")


#vcftools --gzvcf only_SNPs_DP6.recode.vcf.gz --positions shared_polymorphisms_md08.txt --recode --out shared_all_md08

setwd("~/Desktop/phased_big_new_contigs")

a<-shared_mag_sim_sin
a<-read.table("shared_polymorphisms_all_md05_dp10_gq30_maf01.txt",header = F,col.names = c("CHROM","POS"))

a1<-cbind(as.character(a$CHROM), as.numeric(a$POS-25),as.numeric(a$POS+25)) ### to get a XXXbp sequence including the share polymorphism in the middle
write.table(a1,file = "sequence_to_extract_50bp.bed",quote = F,sep = "\t",row.names = F,col.names = F)


bed<-read.table("sequence_to_extract_1500bp.bed",header = F)

setwd("~/Desktop/phased_big_new_contigs")
a<-read.table("1k_random_snps.txt",header = F,col.names = c("CHROM","POS"))

setwd("~/Desktop/phased_big_new_contigs/control_phased_vcfs/")

setwd("~/Desktop/phased_big_new_contigs")
bed<-read.table("sequence_to_extract_100bp.bed",header = F)
setwd("~/Desktop/phased_big_new_contigs/phased_vcfs/1500bp/")
for (i in 1:nrow(bed)) {
  cmd<-paste("/usr/local/bin/vcftools --vcf ../../merged_phased_no_indels.recode.vcf --remove-indels --max-alleles 2 --chr \"", bed[i,1], "\" --from-bp ", bed[i,2], " --to-bp ", bed[i,3], " --recode --out \"", bed[i,1], "_", a[i,2],"\"", sep = "")
  system(cmd)
}

## test_11F imputation
a<-a[a$CHROM=="000011F|quiver",]
bed<-bed[bed$V1=="000011F|quiver",]

setwd("~/Desktop/phased_big_new_contigs/test_11F_imputation/")
for (i in 1:nrow(bed)) {
  cmd<-paste("/usr/local/bin/vcftools --vcf /Users/luca/Desktop/phased_big_new_contigs/phased_dif_missingness/000011F.25perc.vcf --remove-indels --max-alleles 2 --chr \"", bed[i,1], "\" --from-bp ", bed[i,2], " --to-bp ", bed[i,3], " --recode --out \"", bed[i,1], "_", a[i,2],"\"", sep = "")
  system(cmd)
}

setwd("~/Desktop/phased_big_new_contigs/test_11F_imputation/")
for (i in 1:nrow(bed)) {
  cmd<-paste("/usr/local/bin/vcftools --vcf /Users/luca/Desktop/phased_big_new_contigs/phased_dif_missingness/000011F.50perc.vcf --remove-indels --max-alleles 2 --chr \"", bed[i,1], "\" --from-bp ", bed[i,2], " --to-bp ", bed[i,3], " --recode --out \"", bed[i,1], "_", a[i,2],"\"", sep = "")
  system(cmd)
}

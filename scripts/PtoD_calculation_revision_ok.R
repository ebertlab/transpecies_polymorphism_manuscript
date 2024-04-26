
setwd("~/Desktop/new_vcf/PtoD/all_genes/")
dir<-list.files(pattern="mRNA")
su<-NULL
lung<-read.table("../gene_length_correct.txt",header = T)
for(i in 1:length(dir)){
  setwd(paste("~/Desktop/new_vcf/PtoD/all_genes/",dir[i],sep = ""))
  sim_mag<-read.table(paste("fst_sim_mag_",dir[i],".weir.fst",sep = ""),header = T)
  sin_mag<-read.table(paste("fst_sin_mag_",dir[i],".weir.fst",sep = ""),header = T)
  sin<-read.table(paste("pi_sin_",dir[i],".sites.pi",sep=""),header = T)
  sim<-read.table(paste("pi_sim_",dir[i],".sites.pi",sep=""),header = T)
  mag<-read.table(paste("pi_mag_",dir[i],".sites.pi",sep=""),header = T)
  PtoD_sin<-(length(which(sin$PI!="NaN" & sin$PI>0))/lung[i,2])/((length(which(sin_mag$WEIR_AND_COCKERHAM_FST=="1"))/lung[i,2]) + 1)
  PtoD_sim<-(length(which(sim$PI!="NaN" & sim$PI>0))/lung[i,2])/((length(which(sim_mag$WEIR_AND_COCKERHAM_FST=="1"))/lung[i,2]) + 1)
  PtoD_mag<-(length(which(mag$PI!="NaN" & mag$PI>0))/lung[i,2])/((length(which(sin_mag$WEIR_AND_COCKERHAM_FST=="1"))/lung[i,2]) + 1)
  PtoD_mag1<-(length(which(mag$PI!="NaN" & mag$PI>0))/lung[i,2])/((length(which(sim_mag$WEIR_AND_COCKERHAM_FST=="1"))/lung[i,2]) + 1)
  su<-rbind(su,cbind(dir[i],PtoD_mag,PtoD_mag1,PtoD_sim,PtoD_sin))
  print(i)
}
write.table(su,"../../summary_PtoD_new_corrected_gene_length_revision.txt",quote = F,sep = "\t",row.names = F,col.names = c("gene","magna","magna1","similis","sinensis"))



### for revision Feb 2024

setwd("~/Desktop/Nature_submission/revision/")
a<-read.table("ptod_last.txt",header = T)
bad<-a[a$magna_corr==0 & a$similis_corr==0 & a$sinensis_corr==0 & a$magna1_corr==0,]
`%notin%` <- Negate(`%in%`)
a1<-a[a$gene %notin% bad$gene,]
hist(a1$magna)
ecdf_fun <- function(x,perc) ecdf(x)(perc)
p_mag<-1-ecdf_fun(a1$magna_corr,a1$magna_corr)
p_mag1<-1-ecdf_fun(a1$magna1_corr,a1$magna1_corr)
p_sim<-1-ecdf_fun(a1$similis_corr,a1$similis_corr)
p_sin<-1-ecdf_fun(a1$sinensis_corr,a1$sinensis_corr)
a2<-cbind(a1,p_mag,p_mag1,p_sim,p_sin)

write.table(a2,"PtoD_new_corrected_gene_length_revision.txt",quote = F,sep = "\t",row.names = F,col.names = c("gene","magna","magna1","similis","sinensis","old_percentile_magna","old_percentile_magna1","old_percentile_similis","old_percentile_sinensis","length","magna_corr","magna1_corr","similis_corr","sinensis_corr","percentile_magna","percentile_magna1","percentile_similis","percentile_sinensis"))

high<-a2[a2$p_mag<0.05 & a2$p_sim<0.05 & a2$p_sin<0.05,]
write.table(high,"PtoD_all_genes_top_perc_last.txt",quote = F,sep = "\t",row.names = F,col.names = c("gene","magna","similis","sinensis","old_percentile_magna","old_percentile_similis","old_percentile_sinensis","length","magna_corr","similis_corr","sinensis_corr","percentile_magna","percentile_similis","percentile_sinensis"))



library(ape)
library(distances)
## for background tree generation


setwd("~/Desktop/last_phasing/1500bp/with_sub/")

n<-read.table("000000F|quiver_4431683.res.names",header = F)
b<-c()
c<-c()
for (i in 1:nrow(n)) {
  b<-rbind(paste(as.character(n[i,]),"a",sep = "_"),paste(as.character(n[i,]),"b",sep = "_"))
  c<-rbind(c,b)
}

c

file.names <- dir(path="./", pattern =".seqs$")
pos.names<-dir(path = "./",pattern = ".pos$")
p<-c()
s<-c()
nn<-c()



for(i in 1:length(pos.names)){
  f<-c()
  p<-read.table(pos.names[i],header = F,skip = 1)
  s<-strsplit(pos.names[i],split = "[._]")
  seq<-read.table(file.names[i],header = F)
  seq1<-seq$V1
  nn<-paste(c,substr(seq1,which(p==as.numeric(s[[1]][2])),which(p==as.numeric(s[[1]][2]))),sep = "_")
  nn<-gsub(".bam_","",nn)
  seq2<-vector("list",length(seq1))
  for(x in 1:length(seq1)){
    d<-strsplit(as.character(seq1[x]),split = "")
    seq2[[x]] <-d[[1]][-which(p==s[[1]][2])]
  }
  for (h in 1:length(seq2) ) {
    f<-c(f,paste(seq2[[h]],collapse = ""))
  }
  write.table(cbind(nn,as.data.frame(f)),file = paste(file.names[i],"_red.new.fa",sep = ""),quote = F,row.names = F,col.names = F,sep = "\n")
  print(i)
}


# for i in *.seqs.new.fa; do sed -e 's/RU_AST/D_sim_RU_AST/g' -e 's/RU_B3/D_sim_RU_B3/g' -e 's/D_sinRU/D_sin_RU/g' $i > ${i%.res.seqs.new.fa*}.renamed.fa; done

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

trees <- dir(path="./", pattern ="*seqs_red.new.fa.renamed.fa$")
res<-NULL
for(w in 1:length(trees)){
  dna<-read.dna(trees[w],format = "fasta")
  temp <- as.matrix(dist.dna(dna,model = "N"))
  lab<-colnames(temp)
  lab_sin<-lab[grep("D_sin_*",lab)]
  lab_sim<-lab[grep("D_sim_*",lab)]
  lab_mag<-setdiff(lab,c(lab_sin,lab_sim))
  lab_sin_1<-lab_sin[grep(unique(substrRight(lab,1))[1],substrRight(lab_sin,1))]
  lab_sin_2<-lab_sin[grep(unique(substrRight(lab,1))[2],substrRight(lab_sin,1))]
  lab_sim_1<-lab_sim[grep(unique(substrRight(lab,1))[1],substrRight(lab_sim,1))]
  lab_sim_2<-lab_sim[grep(unique(substrRight(lab,1))[2],substrRight(lab_sim,1))]
  lab_mag_1<-lab_mag[grep(unique(substrRight(lab,1))[1],substrRight(lab_mag,1))]
  lab_mag_2<-lab_mag[grep(unique(substrRight(lab,1))[2],substrRight(lab_mag,1))]
  n_sim<-1000
  esito<-NULL
  for(x in 1:n_sim){
    rand<-c(sample(lab_mag_1,1),sample(lab_mag_2,1),sample(lab_sim_1,1),sample(lab_sim_2,1),sample(lab_sin_1,1),sample(lab_sin_2,1))
    temp2<-temp[rownames(temp)%in%rand,colnames(temp)%in%rand] 
    my_dist<-distances(temp2)
    d_d<-nearest_neighbor_search(distances = my_dist,k = 3,query_indices = NULL)
    vec<-NULL
    for(j in 1:length(my_dist)){
      somma<-sum(my_dist[d_d[,j],j])
      vec<-c(vec,somma)
    }
    row.names(temp2)[d_d[,which.min(vec)]]
    f<-substrRight(row.names(temp2)[d_d[,which.min(vec)]],1)
    esito<-c(esito,mode(sapply(list(f[1], f[2]), all.equal, f[3])) == "logical")
    #plot(nj(temp2))
    print(x)
  }
  res<-c(res,(length(which(esito=="TRUE"))/n_sim)*100)
}  
cbind(trees,res)
write.table(cbind(trees,res),"ratio_allelic_trees_1500_bp.txt",quote = F,sep = "\t",row.names = F,col.names = c("tree","ratio_AT_100bp"))

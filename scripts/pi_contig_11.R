library(ggplot2)
library(zoo)
library(ggpubr)
require(gridExtra)


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

## background 20 longest contigs

setwd("~/Desktop/new_vcf/pop_gen_bg/")
pi_sin_bg<-read.table("pi_sin.windowed.pi",header = T)
pi_sin_bg_1<-as.data.frame(cbind("pi_sin_bg",pi_sin_bg$PI),stringsAsFactors = FALSE)
pi_sin_bg_1$V2<-as.numeric(pi_sin_bg_1$V2)

pi_sim_bg<-read.table("pi_sim.windowed.pi",header = T)
pi_sim_bg_1<-as.data.frame(cbind("pi_sim_bg",pi_sim_bg$PI),stringsAsFactors = FALSE)
pi_sim_bg_1$V2<-as.numeric(pi_sim_bg_1$V2)

pi_mag_bg<-read.table("pi_mag.windowed.pi",header = T)
pi_mag_bg_1<-as.data.frame(cbind("pi_mag_bg",pi_mag_bg$PI),stringsAsFactors = FALSE)
pi_mag_bg_1$V2<-as.numeric(pi_mag_bg_1$V2)

# contig 11

setwd("~/Desktop/new_vcf/pop_gen_contig_11/pi/")

# sinensis

pi_sin<-read.table("pi_sin.windowed.pi",header = T)

`%notin%` <- Negate(`%in%`)

pr_sin1<-pi_sin[pi_sin$BIN_START>2141056-50000 & pi_sin$BIN_START<2141056+50000,]
pr_sin2<-pi_sin[pi_sin$BIN_START>2343857-50000 & pi_sin$BIN_START<2343857+50000,]
pr<-pi_sin[pi_sin$BIN_START>2189740 & pi_sin$BIN_START<2332489,]
bg_1<-pi_sin[pi_sin$BIN_START %notin% pr$BIN_START & pi_sin$BIN_START %notin% pr_sin2$BIN_START ,]
bg_2<-pi_sin[pi_sin$BIN_START %notin% pr$BIN_START & pi_sin$BIN_START %notin% pr_sin1$BIN_START ,]

pr_11_sin1<-as.data.frame(cbind("TSP_11F_sin1",pr_sin1$PI),stringsAsFactors = FALSE)
pr_11_sin1$V2<-as.numeric(pr_11_sin1$V2)

pr_11_sin2<-as.data.frame(cbind("TSP_11F_sin2",pr_sin2$PI),stringsAsFactors = FALSE)
pr_11_sin2$V2<-as.numeric(pr_11_sin2$V2)

bg_11_sin1<-as.data.frame(cbind("bg_11F_sin1",bg_1$PI),stringsAsFactors = FALSE)
bg_11_sin1$V2<-as.numeric(bg_11_sin1$V2)

bg_11_sin2<-as.data.frame(cbind("bg_11F_sin2",bg_2$PI),stringsAsFactors = FALSE)
bg_11_sin2$V2<-as.numeric(bg_11_sin2$V2)

sin_p<-ggplot(pi_sin, aes(BIN_START, PI)) +
        geom_point(pch=21) +
        geom_point(aes(x = 2141056, y = 0.004),bg="green",pch=24, cex=2) +
        geom_point(aes(x = 2343857, y = 0.004),bg="green",pch=24, cex=2) +
        geom_line(aes(y=rollmean(PI,10,fill = NA)),col="red") +
        geom_line(aes(y=mean(PI)),lty="dotted") +
        geom_vline(xintercept = 2189740, lty="dashed", colour="blue") +
        geom_vline(xintercept = 2332489, lty="dashed", colour="blue") +
        geom_vline(xintercept = 2329948, lty="dashed", colour="brown") +
        geom_vline(xintercept = 2361178, lty="dashed", colour="brown") +
        theme_bw() 
sin_p

c11_sin<-as.data.frame(rbind(pr_11_sin1, pr_11_sin2,bg_11_sin1, bg_11_sin2,pi_sin_bg_1))
colnames(c11_sin)<-c("cat","PI")
c11_sin$cat<-factor(c11_sin$cat, levels = c("TSP_11F_sin1","bg_11F_sin1","TSP_11F_sin2","bg_11F_sin2","pi_sin_bg"))

my_comparisons_sinensis <- list(c("TSP_11F_sin1", "bg_11F_sin1"),c("TSP_11F_sin2", "bg_11F_sin2"), c("TSP_11F_sin1","pi_sin_bg"), c("TSP_11F_sin2","pi_sin_bg"))
sinensis_c_p<-ggplot(c11_sin, aes(x=cat,y=PI,fill=cat)) + ylab("PI - D. sinensis") + xlab("all") + geom_violin() + stat_summary(fun.data=data_summary, color="black") + theme_bw() + theme(legend.position = "none") + stat_compare_means(comparisons = my_comparisons_sinensis) 
sinensis_c_p

# similis

pi_sim<-read.table("pi_sim.windowed.pi",header = T)

`%notin%` <- Negate(`%in%`)

pr_sim1<-pi_sim[pi_sim$BIN_START>2141056-50000 & pi_sim$BIN_START<2141056+50000,]
pr_sim2<-pi_sim[pi_sim$BIN_START>2343857-50000 & pi_sim$BIN_START<2343857+50000,]
pr<-pi_sim[pi_sim$BIN_START>2189740 & pi_sim$BIN_START<2332489,]
bg_1<-pi_sim[pi_sim$BIN_START %notin% pr$BIN_START & pi_sim$BIN_START %notin% pr_sim2$BIN_START ,]
bg_2<-pi_sim[pi_sim$BIN_START %notin% pr$BIN_START & pi_sim$BIN_START %notin% pr_sim1$BIN_START ,]

pr_11_sim1<-as.data.frame(cbind("TSP_11F_sim1",pr_sim1$PI),stringsAsFactors = FALSE)
pr_11_sim1$V2<-as.numeric(pr_11_sim1$V2)

pr_11_sim2<-as.data.frame(cbind("TSP_11F_sim2",pr_sim2$PI),stringsAsFactors = FALSE)
pr_11_sim2$V2<-as.numeric(pr_11_sim2$V2)

bg_11_sim1<-as.data.frame(cbind("bg_11F_sim1",bg_1$PI),stringsAsFactors = FALSE)
bg_11_sim1$V2<-as.numeric(bg_11_sim1$V2)

bg_11_sim2<-as.data.frame(cbind("bg_11F_sim2",bg_2$PI),stringsAsFactors = FALSE)
bg_11_sim2$V2<-as.numeric(bg_11_sim2$V2)

sim_p<-ggplot(pi_sim, aes(BIN_START, PI)) + ylim(0, 0.0045) +
  geom_point(pch=21) +
  geom_point(aes(x = 2141056, y = 0.004),bg="green",pch=24, cex=2) +
  geom_point(aes(x = 2343857, y = 0.004),bg="green",pch=24, cex=2) +
  geom_line(aes(y=rollmean(PI,10,fill = NA)),col="red") +
  geom_line(aes(y=mean(PI)),lty="dotted") +
  geom_vline(xintercept = 2189740, lty="dashed", colour="blue") +
  geom_vline(xintercept = 2332489, lty="dashed", colour="blue") +
  geom_vline(xintercept = 2329948, lty="dashed", colour="brown") +
  geom_vline(xintercept = 2361178, lty="dashed", colour="brown") +
  theme_bw() 
sim_p

c11_sim<-as.data.frame(rbind(pr_11_sim1, pr_11_sim2,bg_11_sim1, bg_11_sim2,pi_sim_bg_1))
colnames(c11_sim)<-c("cat","PI")
c11_sim$cat<-factor(c11_sim$cat, levels = c("TSP_11F_sim1","bg_11F_sim1","TSP_11F_sim2","bg_11F_sim2","pi_sim_bg"))

my_comparisons_similis <- list(c("TSP_11F_sim1", "bg_11F_sim1"),c("TSP_11F_sim2", "bg_11F_sim2"), c("TSP_11F_sim1","pi_sim_bg"), c("TSP_11F_sim2","pi_sim_bg"))
similis_c_p<-ggplot(c11_sim, aes(x=cat,y=PI,fill=cat)) + ylab("PI - D. similis") + xlab("all") + geom_violin() + stat_summary(fun.data=data_summary, color="black") + theme_bw() + theme(legend.position = "none") + stat_compare_means(comparisons = my_comparisons_similis) 
similis_c_p


# magna

pi_mag<-read.table("pi_mag.windowed.pi",header = T)

`%notin%` <- Negate(`%in%`)

pr_mag1<-pi_mag[pi_mag$BIN_START>2141056-50000 & pi_mag$BIN_START<2141056+50000,]
pr_mag2<-pi_mag[pi_mag$BIN_START>2343857-50000 & pi_mag$BIN_START<2343857+50000,]
pr<-pi_mag[pi_mag$BIN_START>2189740 & pi_mag$BIN_START<2332489,]
bg_1<-pi_mag[pi_mag$BIN_START %notin% pr$BIN_START & pi_mag$BIN_START %notin% pr_mag2$BIN_START ,]
bg_2<-pi_mag[pi_mag$BIN_START %notin% pr$BIN_START & pi_mag$BIN_START %notin% pr_mag1$BIN_START ,]

pr_11_mag1<-as.data.frame(cbind("TSP_11F_mag1",pr_mag1$PI),stringsAsFactors = FALSE)
pr_11_mag1$V2<-as.numeric(pr_11_mag1$V2)

pr_11_mag2<-as.data.frame(cbind("TSP_11F_mag2",pr_mag2$PI),stringsAsFactors = FALSE)
pr_11_mag2$V2<-as.numeric(pr_11_mag2$V2)

bg_11_mag1<-as.data.frame(cbind("bg_11F_mag1",bg_1$PI),stringsAsFactors = FALSE)
bg_11_mag1$V2<-as.numeric(bg_11_mag1$V2)

bg_11_mag2<-as.data.frame(cbind("bg_11F_mag2",bg_2$PI),stringsAsFactors = FALSE)
bg_11_mag2$V2<-as.numeric(bg_11_mag2$V2)

mag_p<-ggplot(pi_mag, aes(BIN_START, PI)) +
  geom_point(pch=21) +
  geom_point(aes(x = 2141056, y = 0.012),bg="green",pch=24, cex=2) +
  geom_point(aes(x = 2343857, y = 0.012),bg="green",pch=24, cex=2) +
  geom_line(aes(y=rollmean(PI,10,fill = NA)),col="red") +
  geom_line(aes(y=mean(PI)),lty="dotted") +
  geom_vline(xintercept = 2189740, lty="dashed", colour="blue") +
  geom_vline(xintercept = 2332489, lty="dashed", colour="blue") +
  geom_vline(xintercept = 2329948, lty="dashed", colour="brown") +
  geom_vline(xintercept = 2361178, lty="dashed", colour="brown") +
  theme_bw() 
mag_p

c11_mag<-as.data.frame(rbind(pr_11_mag1, pr_11_mag2,bg_11_mag1, bg_11_mag2,pi_mag_bg_1))
colnames(c11_mag)<-c("cat","PI")
c11_mag$cat<-factor(c11_mag$cat, levels = c("TSP_11F_mag1","bg_11F_mag1","TSP_11F_mag2","bg_11F_mag2","pi_mag_bg"))

my_comparisons_magensis <- list(c("TSP_11F_mag1", "bg_11F_mag1"),c("TSP_11F_mag2", "bg_11F_mag2"), c("TSP_11F_mag1","pi_mag_bg"), c("TSP_11F_mag2","pi_mag_bg"))
magna_c_p<-ggplot(c11_mag, aes(x=cat,y=PI,fill=cat)) + ylab("PI - D. magna") + xlab("all") + geom_violin() + stat_summary(fun.data=data_summary, color="black") + theme_bw() + theme(legend.position = "none") + stat_compare_means(comparisons = my_comparisons_magensis) 
magna_c_p

grid.arrange(mag_p, magna_c_p, sim_p, similis_c_p, sin_p, sinensis_c_p, ncol=2,top="Nucleotide diversity - contig 11")


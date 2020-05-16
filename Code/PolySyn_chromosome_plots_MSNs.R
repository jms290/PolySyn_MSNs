# chromosome rank plots for Polysyndrome project

library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

# SCAx

tmp <- data.frame(chromosome_an_ranks[order(chromosome_an_ranks[,1]),1],chromosome_an_ranks_se[order(chromosome_an_ranks[,1]),1])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,3] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","SE","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","red","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAx_chromosome_ranks_PLS.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  geom_pointrange(aes(x=reorder(Chromosome,-Rank), ymin=(Rank-7620.5)-SE, ymax=(Rank-7620.5)+SE))+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp2 <- data.frame(Ranks=SCA_MSNs_nodes_effects_an_AIBS_PLS_ranks[,1]-7620.5, Chr=chromosome_labels)
tmp2 <- tmp2[-which(tmp2$V1=="" | tmp2$V1=="Un"),]
ggplot(tmp2, aes(tmp2$Ranks,reorder(tmp2$V1,-tmp2$Ranks)))+
  geom_density_ridges(jittered_points = TRUE, quantile_lines = TRUE, quantiles=2,scale = 0.9, alpha = 0.7,
                      vline_size = 1, vline_color = "red", point_size = 0.4,
                      position = position_raincloud(),point_shape = '|')
  
pal <- c("red","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAx_chromosome_ranks_PLS_flip.png",width = 1200,height = 800)
ggplot(tmp, aes(reorder(Chromosome,Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+
  geom_pointrange(aes(x=reorder(Chromosome,Rank), ymin=(Rank-7620.5)-SE, ymax=(Rank-7620.5)+SE))+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_an_ranks[order(chromosome_an_ranks[,1]),1],chromosome_an_ranks_se[order(chromosome_an_ranks[,1]),1])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)==""),])
tmp[,3] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","SE","Chromosome")
pal <- c("blue","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAx_chromosome_ranks_PLS_flip_noreorder_pointrange.png",width = 1200,height = 800)
ggplot(tmp, aes(reorder(Chromosome,as.numeric(as.character(Chromosome))),(Rank-7620.5)))+geom_hline(yintercept = 0,linetype="dashed",size=1,colour="black")+#geom_bar(stat="identity",fill=pal)+
  geom_pointrange(aes(x=reorder(Chromosome,as.numeric(as.character(Chromosome))), ymin=(Rank-7620.5)-SE, ymax=(Rank-7620.5)+SE),colour=pal,size=3)+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48),axis.text.x=element_text(angle=90,vjust = 0.5))
dev.off()

tmp <- data.frame(chromosome_an_ranks_SCT[order(chromosome_an_ranks_SCT[,1]),1])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","red","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCTx_chromosome_ranks_PLS.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

# SCAy

tmp <- data.frame(chromosome_an_ranks[order(chromosome_an_ranks[,2]),2],chromosome_an_ranks_se[order(chromosome_an_ranks[,2]),2])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,3] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","SE","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","red","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAy_chromosome_ranks_PLS.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+
  geom_pointrange(aes(x=reorder(Chromosome,-Rank), ymin=(Rank-7620.5)-SE, ymax=(Rank-7620.5)+SE))+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

pal <- c("red","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAy_chromosome_ranks_PLS_flip.png",width = 1200,height = 800)
ggplot(tmp, aes(reorder(Chromosome,Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+
  geom_pointrange(aes(x=reorder(Chromosome,Rank), ymin=(Rank-7620.5)-SE, ymax=(Rank-7620.5)+SE))+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_an_ranks[order(chromosome_an_ranks[,2]),2],chromosome_an_ranks_se[order(chromosome_an_ranks[,2]),2])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)==""),])
tmp[,3] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","SE","Chromosome")
pal <- c("blue","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAy_chromosome_ranks_PLS_flip_noreorder_pointrange.png",width = 1200,height = 800)
ggplot(tmp, aes(reorder(Chromosome,as.numeric(as.character(Chromosome))),(Rank-7620.5)))+geom_hline(yintercept = 0,linetype="dashed",size=1,colour="black")+#geom_bar(stat="identity",fill=pal)+
  geom_pointrange(aes(x=reorder(Chromosome,as.numeric(as.character(Chromosome))), ymin=(Rank-7620.5)-SE, ymax=(Rank-7620.5)+SE),colour=pal,size=3)+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48),axis.text.x=element_text(angle=90,vjust = 0.5))
dev.off()

tmp <- data.frame(chromosome_an_ranks_SCT[order(chromosome_an_ranks_SCT[,2]),2])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","red","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCTy_chromosome_ranks_PLS.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

# Turner

tmp <- data.frame(chromosome_Turner_ranks[order(chromosome_Turner_ranks[,1]),])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,3] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","SE","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","blue","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Turner_chromosome_ranks_PLS.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+
  geom_pointrange(aes(x=reorder(Chromosome,-Rank), ymin=(Rank-7620.5)-SE, ymax=(Rank-7620.5)+SE))+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

pal <- c(3)
png("/data/NIMH_CHP/PolySyn/Turner_chromosome_ranks_PLS_flip.png",width = 1200,height = 800)
ggplot(tmp, aes(reorder(Chromosome,Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+
  geom_pointrange(aes(x=reorder(Chromosome,Rank), ymin=(Rank-7620.5)-SE, ymax=(Rank-7620.5)+SE))+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_Turner_ranks[order(chromosome_Turner_ranks[,1]),])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)==""),])
tmp[,3] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","SE","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","red")
png("/data/NIMH_CHP/PolySyn/Turner_chromosome_ranks_PLS_flip_noreorder_pointrange.png",width = 1200,height = 800)
ggplot(tmp, aes(reorder(Chromosome,as.numeric(as.character(Chromosome))),(Rank-7620.5)))+geom_hline(yintercept = 0,linetype="dashed",size=1,colour="black")+#geom_bar(stat="identity",fill=pal)+
  geom_pointrange(aes(x=reorder(Chromosome,as.numeric(as.character(Chromosome))), ymin=(Rank-7620.5)-SE, ymax=(Rank-7620.5)+SE),colour=pal,size=3)+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48),axis.text.x=element_text(angle=90,vjust = 0.5))
dev.off()

# Downs

tmp <- data.frame(chromosome_Downs_ranks[order(chromosome_Downs_ranks[,1]),])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,3] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","SE","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","red","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Downs_chromosome_ranks_PLS.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+
  geom_pointrange(aes(x=reorder(Chromosome,-Rank), ymin=(Rank-7620.5)-SE, ymax=(Rank-7620.5)+SE))+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","red","gray")
png("/data/NIMH_CHP/PolySyn/Downs_chromosome_ranks_PLS_flip.png",width = 1200,height = 800)
ggplot(tmp, aes(reorder(Chromosome,Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+
  geom_pointrange(aes(x=reorder(Chromosome,Rank), ymin=(Rank-7620.5)-SE, ymax=(Rank-7620.5)+SE))+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_Downs_ranks[order(chromosome_Downs_ranks[,1]),])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)==""),])
tmp[,3] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","SE","Chromosome")
pal <- c("blue","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Downs_chromosome_ranks_PLS_flip_noreorder_pointrange.png",width = 1200,height = 800)
ggplot(tmp, aes(reorder(Chromosome,as.numeric(as.character(Chromosome))),(Rank-7620.5)))+geom_hline(yintercept = 0,linetype="dashed",size=1,colour="black")+#geom_bar(stat="identity",fill=pal)+
  geom_pointrange(aes(x=reorder(Chromosome,as.numeric(as.character(Chromosome))), ymin=(Rank-7620.5)-SE, ymax=(Rank-7620.5)+SE),colour=pal,size=3)+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48),axis.text.x=element_text(angle=90,vjust = 0.5))
dev.off()

# VCFS 

VELOC_MSNs_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(VELOC_PLS_genes2$loadings[,1]))
tmp <- data.frame(median(VELOC_MSNs_nodes_effects_AIBS_PLS_ranks[match(VELOC_genes_critical_region$V1,colnames(gene_expression)),],na.rm=T))
tmp1 <- sd(VELOC_MSNs_nodes_effects_AIBS_PLS_ranks[match(VELOC_genes_critical_region$V1,colnames(gene_expression)),],na.rm=T)/sqrt(length(VELOC_MSNs_nodes_effects_AIBS_PLS_ranks[match(read.table("Guna_2015_VCFS_ADbreakpoint_Genes.txt")[,1],colnames(gene_expression)),]))
#tmp <- data.frame(median(VELOC_MSNs_nodes_effects_AIBS_PLS_ranks[match(read.table("Guna_2015_VCFS_ADbreakpoint_Genes.txt")[,1],colnames(gene_expression)),],na.rm=T))
#tmp1 <- sd(VELOC_MSNs_nodes_effects_AIBS_PLS_ranks[match(read.table("Guna_2015_VCFS_ADbreakpoint_Genes.txt")[,1],colnames(gene_expression)),],na.rm=T)/sqrt(length(VELOC_MSNs_nodes_effects_AIBS_PLS_ranks[match(read.table("Guna_2015_VCFS_ADbreakpoint_Genes.txt")[,1],colnames(gene_expression)),]))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(VELOC_MSNs_nodes_effects_AIBS_PLS_ranks[,1])[1:90])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(VELOC_MSNs_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels_old=="22"),1])[1:90])
}
png("/data/NIMH_CHP/PolySyn/VELOC_chromosome_ranks_PLS.png",width = 1200,height = 1000)
emp <- -(tmp-15241)-7620.5
up <- (-(tmp-15241)-7620.5)-tmp1
down <- (-(tmp-15241)-7620.5)+tmp1
rand_ranks <- -(rand_ranks-15241)-7620.5
rand_ranks2 <- -(rand_ranks2-15241)-7620.5
ggplot()+theme_classic()+geom_density(aes(rand_ranks),size=2)+geom_density(aes(rand_ranks2),colour="gray",size=2)+
  geom_vline(aes(xintercept=emp[,1]), colour="red",size=3)+geom_vline(aes(xintercept=up[,1]), colour="red",size=1,linetype="dashed")+geom_vline(aes(xintercept=down[,1]), colour="red",size=1,linetype="dashed")+
  ylab("Density")+xlab("PLS Rank")+theme(text=element_text(size=48))+
  scale_y_continuous(expand = c(0,0))
dev.off()

# WAGR

tmp <- data.frame(median(WAGR_MSNs_nodes_effects_AIBS_PLS_ranks[match(WAGR_genes$V1,colnames(gene_expression_v3)),],na.rm=T))
tmp1 <- sd(WAGR_MSNs_nodes_effects_AIBS_PLS_ranks[match(WAGR_genes$V1,colnames(gene_expression_v3)),],na.rm=T)/sqrt(length(WAGR_MSNs_nodes_effects_AIBS_PLS_ranks[match(WAGR_genes$V1,colnames(gene_expression_v3)),]))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(WAGR_MSNs_nodes_effects_AIBS_PLS_ranks[,1])[1:59])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(WAGR_MSNs_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels_v2[,6]=="11"),1])[1:59])
}
png("/data/NIMH_CHP/PolySyn/WAGR_chromosome_ranks_PLS.png",width = 850,height = 600)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),size=2)+geom_density(aes(rand_ranks2-7620.5),colour="gray",size=2)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+geom_vline(aes(xintercept=(tmp-7620.5)-tmp1), colour="blue",size=1,linetype="dashed")+geom_vline(aes(xintercept=(tmp-7620.5)+tmp1), colour="blue",size=1,linetype="dashed")+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

WAGR_breakpoints <- read.csv("WAGR_sample_breakpoints.csv")
WAGR_breakpoints <- WAGR_breakpoints[-c(1:6),]
library(GenomicRanges)
library(Homo.sapiens)
geneRanges <- 
  function(db, column="ENTREZID")
  {
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
  }
splitColumnByOverlap <-
  function(query, subject, column="ENTREZID", ...)
  {
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
  }
splitByOverlap <- function(x, f, column="ENTREZID", ...)
{
  olaps <- findOverlaps(query, subject, ...)
  f1 <- factor(subjectHits(olaps), levels=seq_len(subjectLength(olaps)))
  splitAsList(mcols(x)[[column]][queryHits(olaps)], f1)
}
cnv_min <- median(WAGR_breakpoints$Del_Start)
cnv_max <- median(WAGR_breakpoints$Del_Stop)
cnv = GRanges(paste0("chr11:",cnv_min, "-", cnv_max))
gns = geneRanges(Homo.sapiens, column="SYMBOL")
overlap = as.vector(splitColumnByOverlap(gns, cnv, "SYMBOL"))
tmp <- data.frame(median(WAGR_MSNs_nodes_effects_AIBS_PLS_ranks[match(overlap$`1`,colnames(gene_expression)),],na.rm=T))
tmp1 <- sd(WAGR_MSNs_nodes_effects_AIBS_PLS_ranks[match(overlap$`1`,colnames(gene_expression)),],na.rm=T)/sqrt(length(VELOC_MSNs_nodes_effects_AIBS_PLS_ranks[match(read.table("Guna_2015_VCFS_ADbreakpoint_Genes.txt")[,1],colnames(gene_expression)),]))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(WAGR_MSNs_nodes_effects_AIBS_PLS_ranks[,1])[1:length(overlap$`1`)])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(WAGR_MSNs_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="11"),1])[1:length(overlap$`1`)])
}
png("/data/NIMH_CHP/PolySyn/WAGR_chromosome_ranks_PLS_breakpoint_median.png",width = 1200,height = 1000)
emp <- -(tmp-15043)-7521.5
up <- (-(tmp-15043)-7521.5)-tmp1
down <- (-(tmp-15043)-7521.5)+tmp1
rand_ranks <- -(rand_ranks-15043)-7521.5
rand_ranks2 <- -(rand_ranks2-15043)-7521.5
ggplot()+theme_classic()+geom_density(aes(rand_ranks),size=2)+geom_density(aes(rand_ranks2),colour="gray",size=2)+
  geom_vline(aes(xintercept=emp[,1]), colour="red",size=3)+geom_vline(aes(xintercept=up[,1]), colour="red",size=1,linetype="dashed")+geom_vline(aes(xintercept=down[,1]), colour="red",size=1,linetype="dashed")+
  ylab("Density")+xlab("PLS Rank")+theme(text=element_text(size=48))+
  scale_y_continuous(expand = c(0,0))
dev.off()

cnv_min <- max(WAGR_breakpoints$Del_Start)
cnv_max <- min(WAGR_breakpoints$Del_Stop)
cnv = GRanges(paste0("chr11:",cnv_min, "-", cnv_max))
gns = geneRanges(Homo.sapiens, column="SYMBOL")
overlap = as.vector(splitColumnByOverlap(gns, cnv, "SYMBOL"))
tmp <- data.frame(median(WAGR_MSNs_nodes_effects_AIBS_PLS_ranks[match(overlap$`1`,colnames(gene_expression)),],na.rm=T))
tmp1 <- sd(VELOC_MSNs_nodes_effects_AIBS_PLS_ranks[match(overlap$`1`,colnames(gene_expression)),],na.rm=T)/sqrt(length(VELOC_MSNs_nodes_effects_AIBS_PLS_ranks[match(read.table("Guna_2015_VCFS_ADbreakpoint_Genes.txt")[,1],colnames(gene_expression)),]))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(WAGR_MSNs_nodes_effects_AIBS_PLS_ranks[,1])[1:length(overlap$`1`)])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(WAGR_MSNs_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="11"),1])[1:length(overlap$`1`)])
}
png("/data/NIMH_CHP/PolySyn/WAGR_chromosome_ranks_PLS_breakpoint_max.png",width = 1200,height = 1000)
emp <- -(tmp-15241)-7620.5
up <- (-(tmp-15241)-7620.5)-tmp1
down <- (-(tmp-15241)-7620.5)+tmp1
rand_ranks <- -(rand_ranks-15241)-7620.45
rand_ranks2 <- -(rand_ranks2-15241)-7620.45
ggplot()+theme_classic()+geom_density(aes(rand_ranks),size=2)+geom_density(aes(rand_ranks2),colour="gray",size=2)+
  geom_vline(aes(xintercept=emp[,1]), colour="orange",size=3)+geom_vline(aes(xintercept=up[,1]), colour="orange",size=1,linetype="dashed")+geom_vline(aes(xintercept=down[,1]), colour="orange",size=1,linetype="dashed")+
  ylab("Density")+xlab("PLS Rank")+theme(text=element_text(size=48))+
  scale_y_continuous(expand = c(0,0))
dev.off()






# PAX6

tmp <- data.frame(PAX6_MSNs_nodes_effects_AIBS_ranks[which(colnames(gene_expression)=="PAX6"),])
# tmp1 <- c(SCA_MSNs_nodes_effects_an_AIBS_PLS_ranks[which(colnames(gene_expression)=="PAX6"),1:2],
#                    Turner_MSNs_nodes_effects_AIBS_PLS_ranks[which(colnames(gene_expression)=="PAX6"),],
#                    Downs_MSNs_nodes_effects_AIBS_PLS_ranks[which(colnames(gene_expression)=="PAX6"),],
#                    VELOC_MSNs_nodes_effects_AIBS_ranks[which(colnames(gene_expression)=="PAX6"),],
#                    WAGR_MSNs_nodes_effects_AIBS_ranks[which(colnames(gene_expression)=="PAX6"),])
png("/data/NIMH_CHP/PolySyn/PAX6_chromosome_ranks_PLS.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(PAX6_MSNs_nodes_effects_AIBS_ranks-7620.5),size=2)+
  geom_density(aes(PAX6_MSNs_nodes_effects_AIBS_ranks[which(chromosome_labels=="11"),]-7620.5),size=2,colour="gray")+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))+scale_x_continuous(breaks=c(-5000,0,5000))
dev.off()


# ASD and Down differential expression

tmp2 <- data.frame(rank(-abideI_PLS_genes$loadings[,1]))
tmp <- data.frame(median(tmp2[match(Parikshak_Nature16_S2a_genes$HGNC.Symbol[which(Parikshak_Nature16_S2a_genes$log2.FC..ASD.vs.CTL > 0 & Parikshak_Nature16_S2a_genes$FDR.adjusted.P.value..ASD.vs.CTL <= 0.05)],colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(tmp2[,1])[1:454])
}
tmp1 <- data.frame(median(tmp2[match(Parikshak_Nature16_S2a_genes$HGNC.Symbol[which(Parikshak_Nature16_S2a_genes$log2.FC..ASD.vs.CTL < 0 & Parikshak_Nature16_S2a_genes$FDR.adjusted.P.value..ASD.vs.CTL <= 0.05)],colnames(gene_expression)),],na.rm=T))
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(tmp2[,1])[1:427])
}
png("/data/NIMH_CHP/PolySyn/abideI_DGE_ranks_PLS.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),colour="white",fill=alpha("blue",0.3))+geom_density(aes(rand_ranks2-7620.5),colour="white",fill=alpha("red",0.3))+
  geom_vline(aes(xintercept=tmp1-7620.5), colour="blue",size=3)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="red",size=3)+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp2 <- data.frame(rank(-abideI_PLS_genes$loadings[,1]))
tmp <- data.frame(median(tmp2[match(Parikshak_Nature16_S2a_genes$HGNC.Symbol[which(Parikshak_Nature16_S2a_genes$log2.FC..ASD.vs.CTL > 0 & Parikshak_Nature16_S2a_genes$P.value..ASD.vs.CTL <= 0.05)],colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(tmp2[,1])[1:1493])
}
tmp1 <- data.frame(median(tmp2[match(Parikshak_Nature16_S2a_genes$HGNC.Symbol[which(Parikshak_Nature16_S2a_genes$log2.FC..ASD.vs.CTL < 0 & Parikshak_Nature16_S2a_genes$P.value..ASD.vs.CTL <= 0.05)],colnames(gene_expression)),],na.rm=T))
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(tmp2[,1])[1:1827])
}
png("/data/NIMH_CHP/PolySyn/abideI_DGE_noFDR_p_ranks_PLS.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),colour="white",fill=alpha("blue",0.3))+geom_density(aes(rand_ranks2-7620.5),colour="white",fill=alpha("red",0.3))+
  geom_vline(aes(xintercept=tmp1-7620.5), colour="blue",size=3)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="red",size=3)+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

abideI_MSNs_nodes_effects_UM_AIBS_PLS_ranks <- data.frame(rank(-abideI_PLS_genes_UM$loadings[,1]))
tmp <- data.frame(median(abideI_MSNs_nodes_effects_UM_AIBS_PLS_ranks[match(Parikshak_Nature16_S2a_genes$HGNC.Symbol[which(Parikshak_Nature16_S2a_genes$log2.FC..ASD.vs.CTL > 0 & Parikshak_Nature16_S2a_genes$FDR.adjusted.P.value..ASD.vs.CTL <= 0.05)],colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(abideI_MSNs_nodes_effects_UM_AIBS_PLS_ranks[,1])[1:454])
}
tmp1 <- data.frame(median(abideI_MSNs_nodes_effects_UM_AIBS_PLS_ranks[match(Parikshak_Nature16_S2a_genes$HGNC.Symbol[which(Parikshak_Nature16_S2a_genes$log2.FC..ASD.vs.CTL < 0 & Parikshak_Nature16_S2a_genes$FDR.adjusted.P.value..ASD.vs.CTL <= 0.05)],colnames(gene_expression)),],na.rm=T))
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(abideI_MSNs_nodes_effects_UM_AIBS_PLS_ranks[,1])[1:427])
}
png("/data/NIMH_CHP/PolySyn/abideI_UM_DGE_ranks_PLS.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),colour="white",fill=alpha("blue",0.3))+geom_density(aes(rand_ranks2-7620.5),colour="white",fill=alpha("red",0.3))+
  geom_vline(aes(xintercept=tmp1-7620.5), colour="blue",size=3)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="red",size=3)+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp2 <- plsr(DownsUK_MSNs_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=1)
tmp2 <- data.frame(rank(tmp2$loadings[,1]))
tmp <- data.frame(median(tmp2[match(Sestan_Neuron16_DownsDE$gene_symbol[which(Sestan_Neuron16_DownsDE$fold_difference..log2. > 0)],colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(tmp2[,1])[1:315])
}
tmp1 <- data.frame(median(tmp2[match(Sestan_Neuron16_DownsDE$gene_symbol[which(Sestan_Neuron16_DownsDE$fold_difference..log2. < 0)],colnames(gene_expression)),],na.rm=T))
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(tmp2[,1])[1:350])
}
png("/data/NIMH_CHP/PolySyn/Downs_DGE_noFDR_ranks_PLS.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),colour="white",fill=alpha("blue",0.3))+geom_density(aes(rand_ranks2-7620.5),colour="white",fill=alpha("red",0.3))+
  geom_vline(aes(xintercept=tmp1-7620.5), colour="blue",size=3)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="red",size=3)+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp2 <- plsr(DownsUK_MSNs_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=1)
tmp2 <- data.frame(rank(tmp2$loadings[,1]))
tmp <- data.frame(median(tmp2[match(Sestan_Neuron16_DownsDE$gene_symbol[which(Sestan_Neuron16_DownsDE$fold_difference..log2. > 0 & Sestan_Neuron16_DownsDE$qval <= 0.05)],colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(tmp2[,1])[1:92])
}
tmp1 <- data.frame(median(tmp2[match(Sestan_Neuron16_DownsDE$gene_symbol[which(Sestan_Neuron16_DownsDE$fold_difference..log2. < 0 & Sestan_Neuron16_DownsDE$qval <= 0.05)],colnames(gene_expression)),],na.rm=T))
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(tmp2[,1])[1:87])
}
png("/data/NIMH_CHP/PolySyn/Downs_DGE_ranks_PLS.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),colour="white",fill=alpha("blue",0.3))+geom_density(aes(rand_ranks2-7620.5),colour="white",fill=alpha("red",0.3))+
  geom_vline(aes(xintercept=tmp1-7620.5), colour="blue",size=3)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="red",size=3)+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp <- median(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. < 0)]), colnames(gene_expression)),],na.rm=T)
tmp <- data.frame(tmp,sd(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. < 0)]), colnames(gene_expression)),],na.rm=T)/sqrt(length(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. < 0)]), colnames(gene_expression)),])))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[,1])[1:784])
}
tmp1 <- median(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0)]), colnames(gene_expression)),],na.rm=T)
tmp1 <- data.frame(tmp1,sd(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0)]), colnames(gene_expression)),],na.rm=T)/sqrt(length(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0)]), colnames(gene_expression)),])))
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[,1])[1:729])
}
tmp2 <- sapply(unique(Sestan_Neuron16_DownsDE_s4$perm), function(x) median(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0 & Sestan_Neuron16_DownsDE_s4$perm==x)]), colnames(gene_expression)),],na.rm=T))
tmp3 <- sapply(unique(Sestan_Neuron16_DownsDE_s4$perm), function(x) median(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. < 0 & Sestan_Neuron16_DownsDE_s4$perm==x)]), colnames(gene_expression)),],na.rm=T))
png("/data/NIMH_CHP/PolySyn/Downs_DGE_ranks_PLS_v2.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),colour="white",fill=alpha("blue",0.3))+geom_density(aes(rand_ranks2-7620.5),colour="white",fill=alpha("red",0.3))+
  geom_vline(aes(xintercept=tmp1[,1]-7620.5), colour="red",size=3)+geom_vline(aes(xintercept=(tmp1[,1]-7620.5)-tmp1[,2]), colour="red",size=1,linetype="dashed")+geom_vline(aes(xintercept=(tmp1[,1]-7620.5)+tmp1[,2]), colour="red",size=1,linetype="dashed")+
  geom_vline(aes(xintercept=tmp[,1]-7620.5), colour="blue",size=3)+geom_vline(aes(xintercept=(tmp[,1]-7620.5)-tmp[,2]), colour="blue",size=1,linetype="dashed")+geom_vline(aes(xintercept=(tmp[,1]-7620.5)+tmp[,2]), colour="blue",size=1,linetype="dashed")+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp4 <- median(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0)]), colnames(gene_expression))[which(chromosome_labels_v2[,6]!="21")],],na.rm=T)
tmp4 <- data.frame(tmp4,sd(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0)]), colnames(gene_expression))[which(chromosome_labels!="21")],],na.rm=T)/sqrt(length(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0 & Sestan_Neuron16_DownsDE_s4$qval<=0.05)]), colnames(gene_expression))[which(chromosome_labels_v2[,6]!="21")],])))

png("/data/NIMH_CHP/PolySyn/Downs_DGE_ranks_PLS_v2_perm.png",width = 1200,height = 800)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),colour="white",fill=alpha("blue",0.3))+geom_density(aes(rand_ranks2-7620.5),colour="white",fill=alpha("red",0.3))+
  geom_rect(aes(xmin=min(tmp2-7620.5),xmax=max(tmp2-7620.5), ymin=0, ymax=0.0015),fill="red",alpha=0.25)+
  geom_rect(aes(xmin=min(tmp3-7620.5),xmax=max(tmp3-7620.5), ymin=0, ymax=0.0015),fill="blue",alpha=0.25)+
  #geom_vline(aes(xintercept=tmp2-7620.5),size=1, colour="coral")+
  #geom_vline(aes(xintercept=tmp3-7620.5),size=1, colour="lightblue")+
  #geom_vline(aes(xintercept=tmp1[,1]-7620.5), colour="red",size=3, alpha=0.5)+geom_vline(aes(xintercept=(tmp1[,1]-7620.5)-tmp1[,2]), colour="red",size=1,linetype="dashed")+geom_vline(aes(xintercept=(tmp1[,1]-7620.5)+tmp1[,2]), colour="red",size=1,linetype="dashed")+
  geom_vline(aes(xintercept=tmp[,1]-7620.5), colour="blue",size=3)+geom_vline(aes(xintercept=(tmp[,1]-7620.5)-tmp[,2]), colour="blue",size=1,linetype="dashed")+geom_vline(aes(xintercept=(tmp[,1]-7620.5)+tmp[,2]), colour="blue",size=1,linetype="dashed")+
  geom_vline(aes(xintercept=tmp4[,1]-7620.5), colour="red",size=3)+geom_vline(aes(xintercept=(tmp4[,1]-7620.5)-tmp4[,2]), colour="red",size=1,linetype="dashed")+geom_vline(aes(xintercept=(tmp4[,1]-7620.5)+tmp4[,2]), colour="red",size=1,linetype="dashed")+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))+scale_x_continuous(breaks = c(-2500,0,2500))
dev.off()


tmp1 <- median(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0 & Sestan_Neuron16_DownsDE_s4$qval<=0.05)]), colnames(gene_expression))[which(chromosome_labels_v2[,6]!="21")],],na.rm=T)
tmp1 <- data.frame(tmp1,sd(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0 & Sestan_Neuron16_DownsDE_s4$qval<=0.05)]), colnames(gene_expression))[which(chromosome_labels!="21")],],na.rm=T)/sqrt(length(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0 & Sestan_Neuron16_DownsDE_s4$qval<=0.05)]), colnames(gene_expression))[which(chromosome_labels_v2[,6]!="21")],])))
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[,1])[1:690])
}
png("/data/NIMH_CHP/PolySyn/Downs_DGE_ranks_PLS_no21_v2.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),colour="white",fill=alpha("blue",0.3))+geom_density(aes(rand_ranks2-7620.5),colour="white",fill=alpha("red",0.3))+
  geom_vline(aes(xintercept=tmp1[,1]-7620.5), colour="red",size=3)+geom_vline(aes(xintercept=(tmp1[,1]-7620.5)-tmp1[,2]), colour="red",size=1,linetype="dashed")+geom_vline(aes(xintercept=(tmp1[,1]-7620.5)+tmp1[,2]), colour="red",size=1,linetype="dashed")+
  geom_vline(aes(xintercept=tmp[,1]-7620.5), colour="blue",size=3)+geom_vline(aes(xintercept=(tmp[,1]-7620.5)-tmp[,2]), colour="blue",size=1,linetype="dashed")+geom_vline(aes(xintercept=(tmp[,1]-7620.5)+tmp[,2]), colour="blue",size=1,linetype="dashed")+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

# tmp <- median(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. < 0)]), colnames(gene_expression)),],na.rm=T)
# tmp <- data.frame(tmp,sd(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. < 0)]), colnames(gene_expression)),],na.rm=T)/sqrt(length(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. < 0)]), colnames(gene_expression)),])))
# rand_ranks <- matrix(0, ncol=1, nrow=10000)
# for (i in 1:10000){
#   rand_ranks[i,] <- median(sample(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[,1])[1:785])
# }
# tmp1 <- median(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0)]), colnames(gene_expression)),],na.rm=T)
# tmp1 <- data.frame(tmp1,sd(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0)]), colnames(gene_expression)),],na.rm=T)/sqrt(length(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0)]), colnames(gene_expression)),])))
# rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
# for (i in 1:10000){
#   rand_ranks2[i,] <- median(sample(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[,1])[1:733])
# }
# png("/data/NIMH_CHP/PolySyn/Downs_DGE_ranks_PLS_v2.png",width = 800,height = 400)
# ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),colour="white",fill=alpha("blue",0.3))+geom_density(aes(rand_ranks2-7620.5),colour="white",fill=alpha("red",0.3))+
#   geom_vline(aes(xintercept=tmp1[,1]-7620.5), colour="red",size=3)+geom_vline(aes(xintercept=(tmp1[,1]-7620.5)-tmp1[,2]), colour="red",size=1,linetype="dashed")+geom_vline(aes(xintercept=(tmp1[,1]-7620.5)+tmp1[,2]), colour="red",size=1,linetype="dashed")+
#   geom_vline(aes(xintercept=tmp[,1]-7620.5), colour="blue",size=3)+geom_vline(aes(xintercept=(tmp[,1]-7620.5)-tmp[,2]), colour="blue",size=1,linetype="dashed")+geom_vline(aes(xintercept=(tmp[,1]-7620.5)+tmp[,2]), colour="blue",size=1,linetype="dashed")+
#   ylab("")+xlab("")+theme(text=element_text(size=36))+
#   scale_y_continuous(expand = c(0,0))
# dev.off()
# 
# 
# tmp1 <- median(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0)]), colnames(gene_expression))[which(chromosome_labels!="21")],],na.rm=T)
# tmp1 <- data.frame(tmp1,sd(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0)]), colnames(gene_expression))[which(chromosome_labels!="21")],],na.rm=T)/sqrt(length(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[match(unique(Sestan_Neuron16_DownsDE_s4$gene_symbol[which(Sestan_Neuron16_DownsDE_s4$fold_difference..log2. > 0)]), colnames(gene_expression))[which(chromosome_labels!="21")],])))
# rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
# for (i in 1:10000){
#   rand_ranks2[i,] <- median(sample(Downs_MSNs_nodes_effects_AIBS_PLS_ranks[,1])[1:702])
# }
# png("/data/NIMH_CHP/PolySyn/Downs_DGE_ranks_PLS_no21_v2.png",width = 800,height = 400)
# ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),colour="white",fill=alpha("blue",0.3))+geom_density(aes(rand_ranks2-7620.5),colour="white",fill=alpha("red",0.3))+
#   geom_vline(aes(xintercept=tmp1[,1]-7620.5), colour="red",size=3)+geom_vline(aes(xintercept=(tmp1[,1]-7620.5)-tmp1[,2]), colour="red",size=1,linetype="dashed")+geom_vline(aes(xintercept=(tmp1[,1]-7620.5)+tmp1[,2]), colour="red",size=1,linetype="dashed")+
#   geom_vline(aes(xintercept=tmp[,1]-7620.5), colour="blue",size=3)+geom_vline(aes(xintercept=(tmp[,1]-7620.5)-tmp[,2]), colour="blue",size=1,linetype="dashed")+geom_vline(aes(xintercept=(tmp[,1]-7620.5)+tmp[,2]), colour="blue",size=1,linetype="dashed")+
#   ylab("")+xlab("")+theme(text=element_text(size=36))+
#   scale_y_continuous(expand = c(0,0))
# dev.off()


chromosome_ranks_list <- list()
for (i in 1:24){
  chromosome_ranks_list[[i]] <- Downs_MSNs_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==as.character(chromosome_labels_reduced[i,])),1]
}
names(chromosome_ranks_list) <- chromosome_labels_reduced$unique.chromosome_labels...c.27..26..2...1.
tmp <- melt(chromosome_ranks_list)
colnames(tmp) <- c("Rank","Chromosome")

tmp1 <- data.frame(chromosome_Downs_ranks[order(chromosome_Downs_ranks)])
tmp2 <- rownames(tmp1)
tmp2 <- tmp2[-which(rownames(tmp1)=="Un"|rownames(tmp1)=="X|Y"|rownames(tmp1)=="")]
tmp1 <- data.frame(tmp1[-which(rownames(tmp1)=="Un"|rownames(tmp1)=="X|Y"|rownames(tmp1)==""),])
tmp1[,2] <- as.factor(tmp2)
tmp1[,1] <- tmp1[,1]
colnames(tmp1) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","red","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))

ggplot(tmp1, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))+
  geom_joy(aes(tmp$Rank-7620.5,as.factor(tmp$Chromosome)),stat = "binline",bins=20,draw_baseline=F,scale=1) #geom_joy(stat = "binline", bins = 10, scale = 0.95, draw_baseline = FALSE)

# X chromomsome inactivation and escape (XCIE) enrichment for SCAx and Turner






# Annotate patient/control differences

tmp <- -Downs_MSNs_nodes_effects_lh[1:152,1]
names(tmp) <- c(1:152)
tmp_mean <- sapply(as.numeric(names(tmp[order(tmp,decreasing = F)])[1:50]), function(x) c(mean(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP=="DS")]),mean(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP!="DS")])))
tmp_se <- sapply(as.numeric(names(tmp[order(tmp,decreasing = F)])[1:50]), function(x) c(sd(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP=="DS")]),sd(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP!="DS")])))/sapply(as.numeric(names(tmp[order(tmp,decreasing = F)])[1:50]), function(x) c(sqrt(length(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP=="DS")])),sqrt(length(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP!="DS")]))))
tmp_mean <- melt(tmp_mean)
tmp_se <- melt(tmp_se)
png("/data/NIMH_CHP/PolySyn/Downs_nodalsimilarity_annotation_neg.png",width = 800,height = 1200)
ggplot(tmp_mean, aes(Var2,value,colour=as.factor(Var1)))+geom_hline(yintercept = 0,linetype="dashed")+geom_point()+
  geom_pointrange(aes(x=Var2, ymin=value-tmp_se$value, ymax=value+tmp_se$value),size=1)+#coord_flip()+
  theme_classic()+scale_color_manual(values = c("blue","black"))+theme(legend.position="none",text=element_text(size=48))+xlab("")+ylab("")
dev.off()

tmp <- -Downs_MSNs_nodes_effects_lh[1:152,1]
names(tmp) <- c(1:152)
tmp_mean <- sapply(as.numeric(names(tmp[order(tmp,decreasing = T)])[1:50]), function(x) c(mean(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP=="DS")]),mean(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP!="DS")])))
tmp_se <- sapply(as.numeric(names(tmp[order(tmp,decreasing = T)])[1:50]), function(x) c(sd(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP=="DS")]),sd(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP!="DS")])))/sapply(as.numeric(names(tmp[order(tmp,decreasing = T)])[1:50]), function(x) c(sqrt(length(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP=="DS")])),sqrt(length(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP!="DS")]))))
tmp_mean <- melt(tmp_mean)
tmp_se <- melt(tmp_se)
png("/data/NIMH_CHP/PolySyn/Downs_nodalsimilarity_annotation_pos.png",width = 800,height = 1200)
ggplot(tmp_mean, aes(Var2,value,colour=as.factor(Var1)))+geom_hline(yintercept = 0,linetype="dashed")+geom_point()+
  geom_pointrange(aes(x=Var2, ymin=value-tmp_se$value, ymax=value+tmp_se$value),size=1)+#coord_flip()+
  theme_classic()+scale_color_manual(values = c("red","black"))+theme(legend.position="none",text=element_text(size=48))+xlab("")+ylab("")
dev.off()

tmp <- Turner_MSNs_nodes_effects_lh[1:152,1]
names(tmp) <- c(1:152)
tmp_mean <- sapply(as.numeric(names(tmp[order(tmp,decreasing = F)])[1:50]), function(x) c(mean(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX=="TURNER")]),mean(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX!="TURNER")])))
tmp_se <- sapply(as.numeric(names(tmp[order(tmp,decreasing = F)])[1:50]), function(x) c(sd(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX=="TURNER")]),sd(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX!="TURNER")])))/sapply(as.numeric(names(tmp[order(tmp,decreasing = F)])[1:50]), function(x) c(sqrt(length(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX=="TURNER")])),sqrt(length(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX!="TURNER")]))))
tmp_mean <- melt(tmp_mean)
tmp_se <- melt(tmp_se)
ggplot(tmp_mean, aes(Var2,value,colour=as.factor(Var1)))+geom_hline(yintercept = 0,linetype="dashed")+geom_point()+
  geom_pointrange(aes(x=Var2, ymin=value-tmp_se$value, ymax=value+tmp_se$value),size=1)+coord_flip()+
  theme_classic()+scale_color_manual(values = c("blue","black"))+theme(legend.position="none",text=element_text(size=48))+xlab("")+ylab("")

tmp <- Turner_MSNs_nodes_effects_lh[1:152,1]
names(tmp) <- c(1:152)
tmp_mean <- sapply(as.numeric(names(tmp[order(tmp,decreasing = T)])[1:50]), function(x) c(mean(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX=="TURNER")]),mean(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX!="TURNER")])))
tmp_se <- sapply(as.numeric(names(tmp[order(tmp,decreasing = T)])[1:50]), function(x) c(sd(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX=="TURNER")]),sd(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX!="TURNER")])))/sapply(as.numeric(names(tmp[order(tmp,decreasing = T)])[1:50]), function(x) c(sqrt(length(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX=="TURNER")])),sqrt(length(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX!="TURNER")]))))
tmp_mean <- melt(tmp_mean)
tmp_se <- melt(tmp_se)
ggplot(tmp_mean, aes(Var2,value,colour=as.factor(Var1)))+geom_hline(yintercept = 0,linetype="dashed")+geom_point()+
  geom_pointrange(aes(x=Var2, ymin=value-tmp_se$value, ymax=value+tmp_se$value),size=1)+coord_flip()+
  theme_classic()+scale_color_manual(values = c("red","black"))+theme(legend.position="none",text=element_text(size=48))+xlab("")+ylab("")

All_MSNs_annotate <- matrix(0, ncol=4 ,nrow=6)
colnames(All_MSNs_annotate) <- c("MorePos","LessNeg","MoreNeg","LessPos")
rownames(All_MSNs_annotate) <- c("+X","+Y","+21","-X","-22q11.2","-11p13")
tmp <- SCA_MSNs_nodes_effects_an_lh[1:152,1]
tmp_mean <- t(sapply(1:152, function(x) c(mean(SCA_MSNs_nodes_lh[x,which(SCA_demo$SCan!=0)]),mean(SCA_MSNs_nodes_lh[x,which(SCA_demo$SCan==0)]))))
All_MSNs_annotate[1,1] <- length(which(tmp_mean[,2] > 0 & SCA_MSNs_nodes_effects_an_lh[1:152,1] > 0))/length(which(SCA_MSNs_nodes_effects_an_lh[1:152,1] > 0))
All_MSNs_annotate[1,2] <- 1-length(which(tmp_mean[,2] > 0 & SCA_MSNs_nodes_effects_an_lh[1:152,1] > 0))/length(which(SCA_MSNs_nodes_effects_an_lh[1:152,1] > 0))
All_MSNs_annotate[1,3] <- length(which(tmp_mean[,2] < 0 & SCA_MSNs_nodes_effects_an_lh[1:152,1] < 0))/length(which(SCA_MSNs_nodes_effects_an_lh[1:152,1] < 0))
All_MSNs_annotate[1,4] <- 1-length(which(tmp_mean[,2] < 0 & SCA_MSNs_nodes_effects_an_lh[1:152,1] < 0))/length(which(SCA_MSNs_nodes_effects_an_lh[1:152,1] < 0))
tmp <- SCA_MSNs_nodes_effects_an_lh[1:152,2]
tmp_mean <- t(sapply(1:152, function(x) c(mean(SCA_MSNs_nodes_lh[x,which(SCA_demo$SCan!=0)]),mean(SCA_MSNs_nodes_lh[x,which(SCA_demo$SCan==0)]))))
All_MSNs_annotate[2,1] <- length(which(tmp_mean[,2] > 0 & SCA_MSNs_nodes_effects_an_lh[1:152,2] > 0))/length(which(SCA_MSNs_nodes_effects_an_lh[1:152,2] > 0))
All_MSNs_annotate[2,2] <- 1-length(which(tmp_mean[,2] > 0 & SCA_MSNs_nodes_effects_an_lh[1:152,2] > 0))/length(which(SCA_MSNs_nodes_effects_an_lh[1:152,2] > 0))
All_MSNs_annotate[2,3] <- length(which(tmp_mean[,2] < 0 & SCA_MSNs_nodes_effects_an_lh[1:152,2] < 0))/length(which(SCA_MSNs_nodes_effects_an_lh[1:152,2] < 0))
All_MSNs_annotate[2,4] <- 1-length(which(tmp_mean[,2] < 0 & SCA_MSNs_nodes_effects_an_lh[1:152,2] < 0))/length(which(SCA_MSNs_nodes_effects_an_lh[1:152,2] < 0))
tmp <- -Downs_MSNs_nodes_effects_lh[1:152,1]
tmp_mean <- t(sapply(1:152, function(x) c(mean(Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP=="DS")]),mean(-Downs_MSNs_nodes_lh[x,which(Downs_demo$GROUP!="DS")]))))
All_MSNs_annotate[3,1] <- length(which(tmp_mean[,2] > 0 & -Downs_MSNs_nodes_effects_lh[1:152,1] > 0))/length(which(-Downs_MSNs_nodes_effects_lh[1:152,1] > 0))
All_MSNs_annotate[3,2] <- 1-length(which(tmp_mean[,2] > 0 & -Downs_MSNs_nodes_effects_lh[1:152,1] > 0))/length(which(-Downs_MSNs_nodes_effects_lh[1:152,1] > 0))
All_MSNs_annotate[3,3] <- length(which(tmp_mean[,2] < 0 & -Downs_MSNs_nodes_effects_lh[1:152,1] < 0))/length(which(-Downs_MSNs_nodes_effects_lh[1:152,1] < 0))
All_MSNs_annotate[3,4] <- 1-length(which(tmp_mean[,2] < 0 & -Downs_MSNs_nodes_effects_lh[1:152,1] < 0))/length(which(-Downs_MSNs_nodes_effects_lh[1:152,1] < 0))
tmp <- Turner_MSNs_nodes_effects_lh[1:152,1]
tmp_mean <- t(sapply(1:152, function(x) c(mean(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX=="TURNER")]),mean(Turner_MSNs_nodes_lh[x,which(Turner_demo$DX!="TURNER")]))))
All_MSNs_annotate[4,1] <- length(which(tmp_mean[,2] > 0 & Turner_MSNs_nodes_effects_lh[1:152,1] > 0))/length(which(Turner_MSNs_nodes_effects_lh[1:152,1] > 0))
All_MSNs_annotate[4,2] <- 1-length(which(tmp_mean[,2] > 0 & Turner_MSNs_nodes_effects_lh[1:152,1] > 0))/length(which(Turner_MSNs_nodes_effects_lh[1:152,1] > 0))
All_MSNs_annotate[4,3] <- length(which(tmp_mean[,2] < 0 & Turner_MSNs_nodes_effects_lh[1:152,1] < 0))/length(which(Turner_MSNs_nodes_effects_lh[1:152,1] < 0))
All_MSNs_annotate[4,4] <- 1-length(which(tmp_mean[,2] < 0 & Turner_MSNs_nodes_effects_lh[1:152,1] < 0))/length(which(Turner_MSNs_nodes_effects_lh[1:152,1] < 0))
tmp <- scale(VELOC_MSNs_nodes_effects_lh[1:152,1])
tmp_mean <- t(sapply(1:152, function(x) c(mean(VELOC_MSNs_nodes_lh[x,which(VELOC_demo$DX=="VELOC CHILD")]),mean(VELOC_MSNs_nodes_lh[x,which(VELOC_demo$DX!="VELOC CHILD")]))))
All_MSNs_annotate[5,1] <- length(which(tmp_mean[,2] > 0 & tmp[1:152,1] > 0))/length(which(tmp[1:152,1] > 0))
All_MSNs_annotate[5,2] <- 1-length(which(tmp_mean[,2] > 0 & tmp[1:152,1] > 0))/length(which(tmp[1:152,1] > 0))
All_MSNs_annotate[5,3] <- length(which(tmp_mean[,2] < 0 & tmp[1:152,1] < 0))/length(which(tmp[1:152,1] < 0))
All_MSNs_annotate[5,4] <- 1-length(which(tmp_mean[,2] < 0 & tmp[1:152,1] < 0))/length(which(tmp[1:152,1] < 0))
tmp <- WAGR_MSNs_nodes_effects_lh[1:152,1]
tmp_mean <- t(sapply(1:152, function(x) c(mean(WAGR_MSNs_nodes_lh_WAGR[x,which(WAGR_demo_WAGR$Group3_strict=="WAGR")]),mean(WAGR_MSNs_nodes_lh_WAGR[x,which(WAGR_demo_WAGR$Group3_strict!="WAGR")]))))
All_MSNs_annotate[6,1] <- length(which(tmp_mean[,2] > 0 & WAGR_MSNs_nodes_effects_lh[1:152,1] > 0))/length(which(WAGR_MSNs_nodes_effects_lh[1:152,1] > 0))
All_MSNs_annotate[6,2] <- 1-length(which(tmp_mean[,2] > 0 & WAGR_MSNs_nodes_effects_lh[1:152,1] > 0))/length(which(WAGR_MSNs_nodes_effects_lh[1:152,1] > 0))
All_MSNs_annotate[6,3] <- length(which(tmp_mean[,2] < 0 & WAGR_MSNs_nodes_effects_lh[1:152,1] < 0))/length(which(WAGR_MSNs_nodes_effects_lh[1:152,1] < 0))
All_MSNs_annotate[6,4] <- 1-length(which(tmp_mean[,2] < 0 & WAGR_MSNs_nodes_effects_lh[1:152,1] < 0))/length(which(WAGR_MSNs_nodes_effects_lh[1:152,1] < 0))

All_MSNs_annotate <- melt(All_MSNs_annotate)
colnames(All_MSNs_annotate) <- c("Group","Effect","Percent")
All_MSNs_annotate[,3] <- All_MSNs_annotate[,3]*100
png("/data/NIMH_CHP/PolySyn/All_MSNs_annotate.png",width = 1000, height = 600)
ggplot(All_MSNs_annotate, aes(Effect, Percent, fill=Effect, colour=Effect))+geom_bar(position=position_dodge(width=0.5),stat="identity")+
  theme_classic()+scale_fill_manual(values=c("red","red","blue","blue"))+scale_colour_manual(values=c("red","blue","blue","red"))+
  facet_wrap(~Group)+theme(text=element_text(size=24),axis.text.x=element_blank(),axis.ticks.x=element_blank())
dev.off()


# trying raincloud plots
source("/data/NIMH_CHP/PolySyn/RainCloudPlots/tutorial_R/R_rainclouds.R")
library(cowplot)

tmp <- data.frame(rank=rank(SCA_PLS_genesX_v2$loadings[,1]), chrom=chromosome_labels[,1])
for (i in levels(tmp$chrom)){
  tmp[which(tmp$chrom==i),3] <- median(tmp$rank[which(tmp$chrom==i)])-7620.5
  tmp[which(tmp$chrom==i),4] <- sd(tmp$rank[which(tmp$chrom==i)])/sqrt(sum(tmp$chrom==i))
}
colnames(tmp)[c(3,4)] <- c("median","se")
tmp <- tmp[-which(tmp$chrom=="" | tmp$chrom=="Un"),]
ggplot(
  tmp,aes(x=reorder(chrom,-median),y=rank-7620.5,colour=median,fill=median))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust=0.25, trim = FALSE, alpha=0.5)+
  geom_point(position = position_jitter(width = .5), size = .25, alpha=0.5)+coord_flip()+
  geom_point(data=tmp, aes(x=reorder(chrom,-median),y=median),size=3)+
  geom_errorbar(data=tmp, aes(x=reorder(chrom,-median),y=rank-7620.5,ymin=median-se,ymax=median+se),size=1)+
  ylab('')+xlab('')+theme_classic()+guides(fill = FALSE, colour = FALSE)+
  scale_colour_distiller(palette="RdBu")+scale_fill_distiller(palette="RdBu")+theme(text=element_text(size=36))+
  ylim(c(-7620.5,7620.5))
ggsave("/data/NIMH_CHP/PolySyn/SCAx_chromosome_ranks_PLS_raincloud.png",width = 10, height = 10)
dev.off()

tmp <- data.frame(rank=rank(SCA_PLS_genesY$loadings[,1]), chrom=chromosome_labels[,1])
for (i in levels(tmp$chrom)){
  tmp[which(tmp$chrom==i),3] <- median(tmp$rank[which(tmp$chrom==i)])-7620.5
  tmp[which(tmp$chrom==i),4] <- sd(tmp$rank[which(tmp$chrom==i)])/sqrt(sum(tmp$chrom==i))
}
colnames(tmp)[c(3,4)] <- c("median","se")
tmp <- tmp[-which(tmp$chrom=="" | tmp$chrom=="Un"),]
ggplot(
  tmp,aes(x=reorder(chrom,-median),y=rank-7620.5,colour=median,fill=median))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust=0.25, trim = FALSE, alpha=0.5)+
  geom_point(position = position_jitter(width = .5), size = .25, alpha=0.5)+coord_flip()+
  geom_point(data=tmp, aes(x=reorder(chrom,-median),y=median),size=3)+
  geom_errorbar(data=tmp, aes(x=reorder(chrom,-median),y=rank-7620.5,ymin=median-se,ymax=median+se),size=1)+
  ylab('')+xlab('')+theme_classic()+guides(fill = FALSE, colour = FALSE)+
  scale_colour_distiller(palette="RdBu")+scale_fill_distiller(palette="RdBu")+theme(text=element_text(size=36))+
  ylim(c(-7620.5,7620.5))
ggsave("/data/NIMH_CHP/PolySyn/SCAy_chromosome_ranks_PLS_raincloud.png",width = 10, height = 10)
dev.off()

tmp <- data.frame(rank=rank(Downs_PLS_genes2$loadings[,1]), chrom=chromosome_labels[,1])
pal <- list()
for (i in levels(tmp$chrom)){
  tmp[which(tmp$chrom==i),3] <- median(tmp$rank[which(tmp$chrom==i)])-7620.5
  tmp[which(tmp$chrom==i),4] <- sd(tmp$rank[which(tmp$chrom==i)])/sqrt(sum(tmp$chrom==i))
  pal[which(tmp$chrom==i)] <- "gray"
}
colnames(tmp)[c(3,4)] <- c("median","se")
pal <- pal[-which(tmp$chrom=="" | tmp$chrom=="Un")]
tmp <- tmp[-which(tmp$chrom=="" | tmp$chrom=="Un"),]
pal[which(tmp$chrom=="21")] <- "red"
ggplot(
  tmp,aes(x=reorder(chrom,-median),y=rank-7620.5,colour=median,fill=median))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust=0.25, trim = FALSE, alpha=0.5)+
  geom_point(position = position_jitter(width = .5), size = .25, alpha=0.5)+coord_flip()+
  geom_point(data=tmp, aes(x=reorder(chrom,-median),y=median),size=3)+
  geom_errorbar(data=tmp, aes(x=reorder(chrom,-median),y=rank-7620.5,ymin=median-se,ymax=median+se),size=1)+
  ylab('')+xlab('')+theme_classic()+guides(fill = FALSE, colour = FALSE)+
  scale_colour_distiller(palette="RdBu")+scale_fill_distiller(palette="RdBu")+theme(text=element_text(size=36))+
  ylim(c(-7620.5,7620.5))
ggsave("/data/NIMH_CHP/PolySyn/Downs_chromosome_ranks_PLS_raincloud.png",width = 10, height = 10)
dev.off()

tmp <- data.frame(rank=rank(Turner_PLS_genes$loadings[,1]), chrom=chromosome_labels[,1])
for (i in levels(tmp$chrom)){
  tmp[which(tmp$chrom==i),3] <- median(tmp$rank[which(tmp$chrom==i)])-7620.5
  tmp[which(tmp$chrom==i),4] <- sd(tmp$rank[which(tmp$chrom==i)])/sqrt(sum(tmp$chrom==i))
}
colnames(tmp)[c(3,4)] <- c("median","se")
tmp <- tmp[-which(tmp$chrom=="" | tmp$chrom=="Un"),]
ggplot(
  tmp,aes(x=reorder(chrom,-median),y=rank-7620.5,colour=median,fill=median))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust=0.25, trim = FALSE, alpha=0.5)+
  geom_point(position = position_jitter(width = .5), size = .25, alpha=0.5)+coord_flip()+
  geom_point(data=tmp, aes(x=reorder(chrom,-median),y=median),size=3)+
  geom_errorbar(data=tmp, aes(x=reorder(chrom,-median),y=rank-7620.5,ymin=median-se,ymax=median+se),size=1)+
  ylab('')+xlab('')+theme_classic()+guides(fill = FALSE, colour = FALSE)+
  scale_colour_distiller(palette="RdBu")+scale_fill_distiller(palette="RdBu")+theme(text=element_text(size=36),axis.ticks = )+
  ylim(c(-7620.5,7620.5))
ggsave("/data/NIMH_CHP/PolySyn/Turner_chromosome_ranks_PLS_raincloud.png",width = 10, height = 10)
dev.off()



# Circos plots
library(RCircos)

# par(mai=c(0.25, 0.25, 0.25, 0.25));
# plot.new();
# plot.window(c(-2.5,2.5), c(-2.5, 2.5));
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Ch

# Chromosome rank tables

tmp <- data.frame(chromosome_an_ranks[order(chromosome_an_ranks[,1]),1],chromosome_an_ranks_se[order(chromosome_an_ranks[,1]),1])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)==""),],table(chromosome_labels)[-c(1,24)][tmp1])
tmp[,3] <- as.factor(tmp1)
tmp[,1] <- (-tmp[,1]+15043)-7521.5
colnames(tmp) <- c("Rank","SE","Chromosome","N")
write.csv(tmp, "SCAx_chromosome_ranks_PLS.csv",row.names=F)

tmp <- data.frame(chromosome_an_ranks[order(chromosome_an_ranks[,2]),2],chromosome_an_ranks_se[order(chromosome_an_ranks[,2]),2])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)==""),],table(chromosome_labels)[-c(1,24)][tmp1])
tmp[,3] <- as.factor(tmp1)
tmp[,1] <- (-tmp[,1]+15043)-7521.5
colnames(tmp) <- c("Rank","SE","Chromosome","N")
write.csv(tmp, "SCAy_chromosome_ranks_PLS.csv",row.names=F)

tmp <- data.frame(chromosome_Turner_ranks[order(chromosome_Turner_ranks[,1]),])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)==""),],table(chromosome_labels)[-c(1,24)][tmp1])
tmp[,3] <- as.factor(tmp1)
tmp[,1] <- (-tmp[,1]+15043)-7521.5
colnames(tmp) <- c("Rank","SE","Chromosome","N")
write.csv(tmp, "Turner_chromosome_ranks_PLS.csv",row.names=F)

tmp <- data.frame(chromosome_Downs_ranks[order(chromosome_Downs_ranks[,1]),])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)==""),],table(chromosome_labels)[-c(1,24)][tmp1])
tmp[,3] <- as.factor(tmp1)
tmp[,1] <- (-tmp[,1]+15043)-7521.5
colnames(tmp) <- c("Rank","SE","Chromosome","N")
write.csv(tmp, "Downs_chromosome_ranks_PLS.csv",row.names=F)






















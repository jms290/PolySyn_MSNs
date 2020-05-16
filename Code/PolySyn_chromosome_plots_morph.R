# chromosome rank plots for Polysyndrome project - cortical thickness

# SCAx

SCA_PLS_genesX_CT <- plsr(SCA_CT_nodes_effects_an_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
SCA_PLS_genesY_CT <- plsr(SCA_CT_nodes_effects_an_lh[1:152,2]~as.matrix(gene_expression),ncomp=2)
SCA_CT_nodes_effects_an_AIBS_PLS_ranks <- data.frame(rank(SCA_PLS_genesX_CT$loadings[,1]),rank(SCA_PLS_genesY_CT$loadings[,1]))
chromosome_an_ranks_CT <- matrix(0, ncol=dim(SCA_MSNs_edges_effects_an_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_an_ranks_CT) <- t(unique(chromosome_labels))
chromosome_an_ranks_CT_se <- matrix(0, ncol=dim(SCA_MSNs_edges_effects_an_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_an_ranks_CT_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_an_ranks_CT[i,1] <- median(SCA_CT_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_an_ranks_CT[i,2] <- median(SCA_CT_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2])
  chromosome_an_ranks_CT_se[i,1] <- sd(SCA_CT_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(SCA_CT_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
  chromosome_an_ranks_CT_se[i,2] <- sd(SCA_CT_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2])/sqrt(length(SCA_CT_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2]))
}

SCA_PLS_genesX_GM <- plsr(SCA_GM_nodes_effects_an_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
SCA_PLS_genesY_GM <- plsr(SCA_GM_nodes_effects_an_lh[1:152,2]~as.matrix(gene_expression),ncomp=2)
SCA_GM_nodes_effects_an_AIBS_PLS_ranks <- data.frame(rank(SCA_PLS_genesX_GM$loadings[,1]),rank(SCA_PLS_genesY_GM$loadings[,1]))
chromosome_an_ranks_GM <- matrix(0, ncol=dim(SCA_MSNs_edges_effects_an_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_an_ranks_GM) <- t(unique(chromosome_labels))
chromosome_an_ranks_GM_se <- matrix(0, ncol=dim(SCA_MSNs_edges_effects_an_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_an_ranks_GM_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_an_ranks_GM[i,1] <- median(SCA_GM_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_an_ranks_GM[i,2] <- median(SCA_GM_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2])
  chromosome_an_ranks_GM_se[i,1] <- sd(SCA_GM_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(SCA_GM_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
  chromosome_an_ranks_GM_se[i,2] <- sd(SCA_GM_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2])/sqrt(length(SCA_GM_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2]))
}

SCA_PLS_genesX_SA <- plsr(SCA_SA_nodes_effects_an_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
SCA_PLS_genesY_SA <- plsr(SCA_SA_nodes_effects_an_lh[1:152,2]~as.matrix(gene_expression),ncomp=2)
SCA_SA_nodes_effects_an_AIBS_PLS_ranks <- data.frame(rank(SCA_PLS_genesX_SA$loadings[,1]),rank(SCA_PLS_genesY_SA$loadings[,1]))
chromosome_an_ranks_SA <- matrix(0, ncol=dim(SCA_MSNs_edges_effects_an_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_an_ranks_SA) <- t(unique(chromosome_labels))
chromosome_an_ranks_SA_se <- matrix(0, ncol=dim(SCA_MSNs_edges_effects_an_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_an_ranks_SA_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_an_ranks_SA[i,1] <- median(SCA_SA_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_an_ranks_SA[i,2] <- median(SCA_SA_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2])
  chromosome_an_ranks_SA_se[i,1] <- sd(SCA_SA_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(SCA_SA_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
  chromosome_an_ranks_SA_se[i,2] <- sd(SCA_SA_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2])/sqrt(length(SCA_SA_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2]))
}

SCA_PLS_genesX_MC <- plsr(SCA_MC_nodes_effects_an_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
SCA_PLS_genesY_MC <- plsr(SCA_MC_nodes_effects_an_lh[1:152,2]~as.matrix(gene_expression),ncomp=2)
SCA_MC_nodes_effects_an_AIBS_PLS_ranks <- data.frame(rank(SCA_PLS_genesX_MC$loadings[,1]),rank(SCA_PLS_genesY_MC$loadings[,1]))
chromosome_an_ranks_MC <- matrix(0, ncol=dim(SCA_MSNs_edges_effects_an_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_an_ranks_MC) <- t(unique(chromosome_labels))
chromosome_an_ranks_MC_se <- matrix(0, ncol=dim(SCA_MSNs_edges_effects_an_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_an_ranks_MC_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_an_ranks_MC[i,1] <- median(SCA_MC_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_an_ranks_MC[i,2] <- median(SCA_MC_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2])
  chromosome_an_ranks_MC_se[i,1] <- sd(SCA_MC_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(SCA_MC_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
  chromosome_an_ranks_MC_se[i,2] <- sd(SCA_MC_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2])/sqrt(length(SCA_MC_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2]))
}

SCA_PLS_genesX_IC <- plsr(SCA_IC_nodes_effects_an_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
SCA_PLS_genesY_IC <- plsr(SCA_IC_nodes_effects_an_lh[1:152,2]~as.matrix(gene_expression),ncomp=2)
SCA_IC_nodes_effects_an_AIBS_PLS_ranks <- data.frame(rank(SCA_PLS_genesX_IC$loadings[,1]),rank(SCA_PLS_genesY_IC$loadings[,1]))
chromosome_an_ranks_IC <- matrix(0, ncol=dim(SCA_MSNs_edges_effects_an_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_an_ranks_IC) <- t(unique(chromosome_labels))
chromosome_an_ranks_IC_se <- matrix(0, ncol=dim(SCA_MSNs_edges_effects_an_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_an_ranks_IC_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_an_ranks_IC[i,1] <- median(SCA_IC_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_an_ranks_IC[i,2] <- median(SCA_IC_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2])
  chromosome_an_ranks_IC_se[i,1] <- sd(SCA_IC_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(SCA_IC_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
  chromosome_an_ranks_IC_se[i,2] <- sd(SCA_IC_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2])/sqrt(length(SCA_IC_nodes_effects_an_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),2]))
}

tmp <- data.frame(chromosome_an_ranks_CT[order(chromosome_an_ranks_CT[,1]),1])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
pal[which(tmp[,2]=="X")-1] <- "red"
png("/data/NIMH_CHP/PolySyn/SCAx_chromosome_ranks_PLS_CT.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_an_ranks_GM[order(chromosome_an_ranks_GM[,1]),1])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "red","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAx_chromosome_ranks_PLS_GM.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_an_ranks_SA[order(chromosome_an_ranks_SA[,1]),1])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","red",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAx_chromosome_ranks_PLS_SA.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_an_ranks_MC[order(chromosome_an_ranks_MC[,1]),1])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","red","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAx_chromosome_ranks_PLS_MC.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_an_ranks_IC[order(chromosome_an_ranks_IC[,1]),1])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","red","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAx_chromosome_ranks_PLS_IC.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

# SCAy

tmp <- data.frame(chromosome_an_ranks_CT[order(chromosome_an_ranks_CT[,2]),2])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","gray","gray","red","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAy_chromosome_ranks_PLS_CT.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_an_ranks_GM[order(chromosome_an_ranks_GM[,2]),2])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","red","gray","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAy_chromosome_ranks_PLS_GM.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_an_ranks_SA[order(chromosome_an_ranks_SA[,2]),2])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","red","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAy_chromosome_ranks_PLS_SA.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_an_ranks_MC[order(chromosome_an_ranks_MC[,2]),2])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","red","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAy_chromosome_ranks_PLS_MC.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_an_ranks_IC[order(chromosome_an_ranks_IC[,2]),2])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","red","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/SCAy_chromosome_ranks_PLS_IC.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

# Turner

Turner_PLS_genes_CT <- plsr(Turner_CT_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
Turner_CT_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(Turner_PLS_genes_CT$loadings[,1]))
#Turner_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(Turner_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(Turner_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
chromosome_ranks_se <- matrix(0, ncol=dim(Turner_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(Turner_CT_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_ranks_se[i,1] <- sd(Turner_CT_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(Turner_CT_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
}
chromosome_Turner_ranks_CT <- data.frame(chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),],chromosome_ranks_se[order(chromosome_ranks[,1], decreasing = F),])

Turner_PLS_genes_GM <- plsr(Turner_GM_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
Turner_GM_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(Turner_PLS_genes_GM$loadings[,1]))
#Turner_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(Turner_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(Turner_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
chromosome_ranks_se <- matrix(0, ncol=dim(Turner_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(Turner_GM_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_ranks_se[i,1] <- sd(Turner_GM_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(Turner_GM_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
}
chromosome_Turner_ranks_GM <- data.frame(chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),],chromosome_ranks_se[order(chromosome_ranks[,1], decreasing = F),])

Turner_PLS_genes_SA <- plsr(Turner_SA_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
Turner_SA_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(Turner_PLS_genes_SA$loadings[,1]))
#Turner_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(Turner_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(Turner_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
chromosome_ranks_se <- matrix(0, ncol=dim(Turner_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(Turner_SA_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_ranks_se[i,1] <- sd(Turner_SA_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(Turner_SA_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
}
chromosome_Turner_ranks_SA <- data.frame(chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),],chromosome_ranks_se[order(chromosome_ranks[,1], decreasing = F),])

Turner_PLS_genes_MC <- plsr(Turner_MC_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
Turner_MC_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(Turner_PLS_genes_MC$loadings[,1]))
#Turner_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(Turner_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(Turner_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
chromosome_ranks_se <- matrix(0, ncol=dim(Turner_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(Turner_MC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_ranks_se[i,1] <- sd(Turner_MC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(Turner_MC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
}
chromosome_Turner_ranks_MC <- data.frame(chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),],chromosome_ranks_se[order(chromosome_ranks[,1], decreasing = F),])

Turner_PLS_genes_IC <- plsr(Turner_IC_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
Turner_IC_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(Turner_PLS_genes_IC$loadings[,1]))
#Turner_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(Turner_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(Turner_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
chromosome_ranks_se <- matrix(0, ncol=dim(Turner_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(Turner_IC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_ranks_se[i,1] <- sd(Turner_IC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(Turner_IC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
}
chromosome_Turner_ranks_IC <- data.frame(chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),],chromosome_ranks_se[order(chromosome_ranks[,1], decreasing = F),])


tmp <- data.frame(chromosome_Turner_ranks_CT[order(chromosome_Turner_ranks_CT)])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","blue","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Turner_chromosome_ranks_PLS_CT.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_Turner_ranks_GM[order(chromosome_Turner_ranks_GM)])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","blue","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Turner_chromosome_ranks_PLS_GM.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_Turner_ranks_SA[order(chromosome_Turner_ranks_SA)])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","blue","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Turner_chromosome_ranks_PLS_SA.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_Turner_ranks_MC[order(chromosome_Turner_ranks_MC)])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","blue","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Turner_chromosome_ranks_PLS_MC.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_Turner_ranks_IC[order(chromosome_Turner_ranks_IC)])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","blue","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Turner_chromosome_ranks_PLS_IC.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

# Downs

Downs_PLS_genes_CT <- plsr(-Downs_CT_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
Downs_CT_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(Downs_PLS_genes_CT$loadings[,1]))
#Downs_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(Downs_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(Downs_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
chromosome_ranks_se <- matrix(0, ncol=dim(Downs_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(Downs_CT_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_ranks_se[i,1] <- sd(Downs_CT_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(Downs_CT_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
}
chromosome_Downs_ranks_CT <- data.frame(chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),],chromosome_ranks_se[order(chromosome_ranks[,1], decreasing = F),])

Downs_PLS_genes_GM <- plsr(-Downs_GM_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
Downs_GM_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(Downs_PLS_genes_GM$loadings[,1]))
#Downs_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(Downs_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(Downs_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
chromosome_ranks_se <- matrix(0, ncol=dim(Downs_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(Downs_GM_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_ranks_se[i,1] <- sd(Downs_GM_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(Downs_GM_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
}
chromosome_Downs_ranks_GM <- data.frame(chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),],chromosome_ranks_se[order(chromosome_ranks[,1], decreasing = F),])

Downs_PLS_genes_SA <- plsr(-Downs_SA_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
Downs_SA_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(Downs_PLS_genes_SA$loadings[,1]))
#Downs_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(Downs_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(Downs_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
chromosome_ranks_se <- matrix(0, ncol=dim(Downs_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(Downs_SA_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_ranks_se[i,1] <- sd(Downs_SA_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(Downs_SA_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
}
chromosome_Downs_ranks_SA <- data.frame(chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),],chromosome_ranks_se[order(chromosome_ranks[,1], decreasing = F),])

Downs_PLS_genes_MC <- plsr(-Downs_MC_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
Downs_MC_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(Downs_PLS_genes_MC$loadings[,1]))
#Downs_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(Downs_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(Downs_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
chromosome_ranks_se <- matrix(0, ncol=dim(Downs_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(Downs_MC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_ranks_se[i,1] <- sd(Downs_MC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(Downs_MC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
}
chromosome_Downs_ranks_MC <- data.frame(chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),],chromosome_ranks_se[order(chromosome_ranks[,1], decreasing = F),])

Downs_PLS_genes_IC <- plsr(-Downs_IC_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
Downs_IC_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(Downs_PLS_genes_IC$loadings[,1]))
#Downs_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(Downs_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(Downs_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
chromosome_ranks_se <- matrix(0, ncol=dim(Downs_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks_se) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(Downs_IC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
  chromosome_ranks_se[i,1] <- sd(Downs_IC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])/sqrt(length(Downs_IC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1]))
}
chromosome_Downs_ranks_IC <- data.frame(chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),],chromosome_ranks_se[order(chromosome_ranks[,1], decreasing = F),])

tmp <- data.frame(chromosome_Downs_ranks_CT[order(chromosome_Downs_ranks_CT)])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","red","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Downs_chromosome_ranks_PLS_CT.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_Downs_ranks_GM[order(chromosome_Downs_ranks_GM)])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","gray","red","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Downs_chromosome_ranks_PLS_GM.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_Downs_ranks_SA[order(chromosome_Downs_ranks_SA)])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray",
         "gray","gray","red","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Downs_chromosome_ranks_PLS_SA.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_Downs_ranks_MC[order(chromosome_Downs_ranks_MC)])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","red","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Downs_chromosome_ranks_PLS_MC.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

tmp <- data.frame(chromosome_Downs_ranks_IC[order(chromosome_Downs_ranks_IC)])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="X|Y"|rownames(tmp)==""),])
tmp[,2] <- as.factor(tmp1)
tmp[,1] <- tmp[,1]
colnames(tmp) <- c("Rank","Chromosome")
pal <- c("gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","red","gray","gray",
         "gray","gray","gray","gray","gray","gray","gray","gray","gray","gray","gray")
png("/data/NIMH_CHP/PolySyn/Downs_chromosome_ranks_PLS_IC.png",width = 800,height = 1200)
ggplot(tmp, aes(reorder(Chromosome,-Rank),Rank-7620.5))+geom_bar(stat="identity",fill=pal)+coord_flip()+
  theme_classic()+ylab("") + xlab("")+theme(text=element_text(size=48))
dev.off()

# VCFS 

VELOC_PLS_genes_CT <- plsr(VELOC_CT_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
VELOC_CT_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(VELOC_PLS_genes_CT$loadings[,1]))
#VELOC_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(VELOC_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(VELOC_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(VELOC_CT_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
}
chromosome_VELOC_ranks_CT <- chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),]

VELOC_PLS_genes_GM <- plsr(VELOC_GM_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
VELOC_GM_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(VELOC_PLS_genes_GM$loadings[,1]))
#VELOC_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(VELOC_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(VELOC_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(VELOC_GM_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
}
chromosome_VELOC_ranks_GM <- chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),]

VELOC_PLS_genes_SA <- plsr(VELOC_SA_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
VELOC_SA_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(VELOC_PLS_genes_SA$loadings[,1]))
#VELOC_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(VELOC_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(VELOC_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(VELOC_SA_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
}
chromosome_VELOC_ranks_SA <- chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),]

VELOC_PLS_genes_MC <- plsr(VELOC_MC_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
VELOC_MC_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(VELOC_PLS_genes_MC$loadings[,1]))
#VELOC_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(VELOC_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(VELOC_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(VELOC_MC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
}
chromosome_VELOC_ranks_MC <- chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),]

VELOC_PLS_genes_IC <- plsr(VELOC_IC_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
VELOC_IC_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(VELOC_PLS_genes_IC$loadings[,1]))
#VELOC_MSNs_nodes_effects_AIBS_ranks <- data.frame(rank(abs(VELOC_MSNs_nodes_effects_AIBS[,1])))
chromosome_ranks <- matrix(0, ncol=dim(VELOC_MSNs_nodes_effects_AIBS)[2], nrow=length(t(unique(chromosome_labels))))
rownames(chromosome_ranks) <- t(unique(chromosome_labels))
for (i in 1:length(t(unique(chromosome_labels)))){
  chromosome_ranks[i,1] <- median(VELOC_IC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels==t(unique(chromosome_labels))[i]),1])
}
chromosome_VELOC_ranks_IC <- chromosome_ranks[order(chromosome_ranks[,1], decreasing = F),]

tmp <- data.frame(median(VELOC_CT_nodes_effects_AIBS_PLS_ranks[match(VELOC_genes_critical_region$V1,colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(VELOC_CT_nodes_effects_AIBS_PLS_ranks[,1])[1:66])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(VELOC_CT_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="22"),1])[1:66])
}
png("/data/NIMH_CHP/PolySyn/VELOC_chromosome_ranks_PLS_CT.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),size=2)+geom_density(aes(rand_ranks2-7620.5),colour="gray",size=2)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+ylab("")+
  xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp <- data.frame(median(VELOC_GM_nodes_effects_AIBS_PLS_ranks[match(VELOC_genes_critical_region$V1,colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(VELOC_GM_nodes_effects_AIBS_PLS_ranks[,1])[1:66])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(VELOC_GM_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="22"),1])[1:66])
}
png("/data/NIMH_CHP/PolySyn/VELOC_chromosome_ranks_PLS_GM.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),size=2)+geom_density(aes(rand_ranks2-7620.5),colour="gray",size=2)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+ylab("")+
  xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp <- data.frame(median(VELOC_SA_nodes_effects_AIBS_PLS_ranks[match(VELOC_genes_critical_region$V1,colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(VELOC_SA_nodes_effects_AIBS_PLS_ranks[,1])[1:66])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(VELOC_SA_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="22"),1])[1:66])
}
png("/data/NIMH_CHP/PolySyn/VELOC_chromosome_ranks_PLS_SA.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),size=2)+geom_density(aes(rand_ranks2-7620.5),colour="gray",size=2)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+ylab("")+
  xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp <- data.frame(median(VELOC_MC_nodes_effects_AIBS_PLS_ranks[match(VELOC_genes_critical_region$V1,colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(VELOC_MC_nodes_effects_AIBS_PLS_ranks[,1])[1:66])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(VELOC_MC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="22"),1])[1:66])
}
png("/data/NIMH_CHP/PolySyn/VELOC_chromosome_ranks_PLS_MC.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),size=2)+geom_density(aes(rand_ranks2-7620.5),colour="gray",size=2)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+ylab("")+
  xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp <- data.frame(median(VELOC_IC_nodes_effects_AIBS_PLS_ranks[match(VELOC_genes_critical_region$V1,colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(VELOC_IC_nodes_effects_AIBS_PLS_ranks[,1])[1:66])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(VELOC_IC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="22"),1])[1:66])
}
png("/data/NIMH_CHP/PolySyn/VELOC_chromosome_ranks_PLS_IC.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),size=2)+geom_density(aes(rand_ranks2-7620.5),colour="gray",size=2)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+ylab("")+
  xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

# WAGR

WAGR_PLS_genes_CT <- plsr(WAGR_CT_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
WAGR_CT_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(WAGR_PLS_genes_CT$loadings[,1]))
WAGR_PLS_genes_GM <- plsr(WAGR_GM_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
WAGR_GM_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(WAGR_PLS_genes_GM$loadings[,1]))
WAGR_PLS_genes_SA <- plsr(WAGR_SA_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
WAGR_SA_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(WAGR_PLS_genes_SA$loadings[,1]))
WAGR_PLS_genes_MC <- plsr(WAGR_MC_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
WAGR_MC_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(WAGR_PLS_genes_MC$loadings[,1]))
WAGR_PLS_genes_IC <- plsr(WAGR_IC_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
WAGR_IC_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(WAGR_PLS_genes_IC$loadings[,1]))

tmp <- data.frame(median(WAGR_CT_nodes_effects_AIBS_PLS_ranks[match(WAGR_genes$V1,colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(WAGR_CT_nodes_effects_AIBS_PLS_ranks[,1])[1:59])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(WAGR_CT_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="11"),1])[1:59])
}
png("/data/NIMH_CHP/PolySyn/WAGR_chromosome_ranks_PLS_CT.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),size=2)+geom_density(aes(rand_ranks2-7620.5),colour="gray",size=2)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+ylab("")+
  xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp <- data.frame(median(WAGR_GM_nodes_effects_AIBS_PLS_ranks[match(WAGR_genes$V1,colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(WAGR_GM_nodes_effects_AIBS_PLS_ranks[,1])[1:59])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(WAGR_GM_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="11"),1])[1:59])
}
png("/data/NIMH_CHP/PolySyn/WAGR_chromosome_ranks_PLS_GM.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),size=2)+geom_density(aes(rand_ranks2-7620.5),colour="gray",size=2)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+ylab("")+
  xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp <- data.frame(median(WAGR_SA_nodes_effects_AIBS_PLS_ranks[match(WAGR_genes$V1,colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(WAGR_SA_nodes_effects_AIBS_PLS_ranks[,1])[1:59])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(WAGR_SA_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="11"),1])[1:59])
}
png("/data/NIMH_CHP/PolySyn/WAGR_chromosome_ranks_PLS_SA.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),size=2)+geom_density(aes(rand_ranks2-7620.5),colour="gray",size=2)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+ylab("")+
  xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp <- data.frame(median(WAGR_MC_nodes_effects_AIBS_PLS_ranks[match(WAGR_genes$V1,colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(WAGR_MC_nodes_effects_AIBS_PLS_ranks[,1])[1:59])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(WAGR_MC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="11"),1])[1:59])
}
png("/data/NIMH_CHP/PolySyn/WAGR_chromosome_ranks_PLS_MC.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),size=2)+geom_density(aes(rand_ranks2-7620.5),colour="gray",size=2)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+ylab("")+
  xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

tmp <- data.frame(median(WAGR_IC_nodes_effects_AIBS_PLS_ranks[match(WAGR_genes$V1,colnames(gene_expression)),],na.rm=T))
rand_ranks <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks[i,] <- median(sample(WAGR_IC_nodes_effects_AIBS_PLS_ranks[,1])[1:59])
}
rand_ranks2 <- matrix(0, ncol=1, nrow=10000)
for (i in 1:10000){
  rand_ranks2[i,] <- median(sample(WAGR_IC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="11"),1])[1:59])
}
png("/data/NIMH_CHP/PolySyn/WAGR_chromosome_ranks_PLS_IC.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(rand_ranks-7620.5),size=2)+geom_density(aes(rand_ranks2-7620.5),colour="gray",size=2)+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+ylab("")+
  xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))
dev.off()

# PAX6

PAX6_PLS_genes_CT <- plsr(PAX6_CT_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
PAX6_CT_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(PAX6_PLS_genes_CT$loadings[,1]))
PAX6_PLS_genes_GM <- plsr(PAX6_GM_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
PAX6_GM_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(PAX6_PLS_genes_GM$loadings[,1]))
PAX6_PLS_genes_SA <- plsr(PAX6_SA_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
PAX6_SA_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(PAX6_PLS_genes_SA$loadings[,1]))
PAX6_PLS_genes_MC <- plsr(PAX6_MC_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
PAX6_MC_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(PAX6_PLS_genes_MC$loadings[,1]))
PAX6_PLS_genes_IC <- plsr(PAX6_IC_nodes_effects_lh[1:152,1]~as.matrix(gene_expression),ncomp=2)
PAX6_IC_nodes_effects_AIBS_PLS_ranks <- data.frame(rank(PAX6_PLS_genes_IC$loadings[,1]))

tmp <- data.frame(PAX6_CT_nodes_effects_AIBS_PLS_ranks[which(colnames(gene_expression)=="PAX6"),])
# tmp1 <- c(SCA_MSNs_nodes_effects_an_AIBS_PLS_ranks[which(colnames(gene_expression)=="PAX6"),1:2],
#                    Turner_MSNs_nodes_effects_AIBS_PLS_ranks[which(colnames(gene_expression)=="PAX6"),],
#                    Downs_MSNs_nodes_effects_AIBS_PLS_ranks[which(colnames(gene_expression)=="PAX6"),],
#                    VELOC_MSNs_nodes_effects_AIBS_ranks[which(colnames(gene_expression)=="PAX6"),],
#                    WAGR_MSNs_nodes_effects_AIBS_ranks[which(colnames(gene_expression)=="PAX6"),])
png("/data/NIMH_CHP/PolySyn/PAX6_chromosome_ranks_PLS_CT.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(PAX6_CT_nodes_effects_AIBS_PLS_ranks-7620.5),size=2)+
  geom_density(aes(PAX6_CT_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="11"),]-7620.5),size=2,colour="gray")+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))+scale_x_continuous(breaks=c(-5000,0,5000))
dev.off()

tmp <- data.frame(PAX6_GM_nodes_effects_AIBS_PLS_ranks[which(colnames(gene_expression)=="PAX6"),])
png("/data/NIMH_CHP/PolySyn/PAX6_chromosome_ranks_PLS_GM.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(PAX6_GM_nodes_effects_AIBS_PLS_ranks-7620.5),size=2)+
  geom_density(aes(PAX6_GM_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="11"),]-7620.5),size=2,colour="gray")+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))+scale_x_continuous(breaks=c(-5000,0,5000))
dev.off()

tmp <- data.frame(PAX6_SA_nodes_effects_AIBS_PLS_ranks[which(colnames(gene_expression)=="PAX6"),])
png("/data/NIMH_CHP/PolySyn/PAX6_chromosome_ranks_PLS_SA.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(PAX6_SA_nodes_effects_AIBS_PLS_ranks-7620.5),size=2)+
  geom_density(aes(PAX6_SA_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="11"),]-7620.5),size=2,colour="gray")+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))+scale_x_continuous(breaks=c(-5000,0,5000))
dev.off()

tmp <- data.frame(PAX6_MC_nodes_effects_AIBS_PLS_ranks[which(colnames(gene_expression)=="PAX6"),])
png("/data/NIMH_CHP/PolySyn/PAX6_chromosome_ranks_PLS_MC.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(PAX6_MC_nodes_effects_AIBS_PLS_ranks-7620.5),size=2)+
  geom_density(aes(PAX6_MC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="11"),]-7620.5),size=2,colour="gray")+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))+scale_x_continuous(breaks=c(-5000,0,5000))
dev.off()

tmp <- data.frame(PAX6_IC_nodes_effects_AIBS_PLS_ranks[which(colnames(gene_expression)=="PAX6"),])
png("/data/NIMH_CHP/PolySyn/PAX6_chromosome_ranks_PLS_IC.png",width = 800,height = 400)
ggplot()+theme_classic()+geom_density(aes(PAX6_MC_nodes_effects_AIBS_PLS_ranks-7620.5),size=2)+
  geom_density(aes(PAX6_MC_nodes_effects_AIBS_PLS_ranks[which(chromosome_labels=="11"),]-7620.5),size=2,colour="gray")+
  geom_vline(aes(xintercept=tmp-7620.5), colour="blue",size=3)+
  ylab("")+xlab("")+theme(text=element_text(size=36))+
  scale_y_continuous(expand = c(0,0))+scale_x_continuous(breaks=c(-5000,0,5000))
dev.off()





# Chromosome rank tables

tmp <- data.frame(chromosome_an_ranks_CT[order(chromosome_an_ranks[,1]),1],chromosome_an_ranks_CT_se[order(chromosome_an_ranks[,1]),1],
                  chromosome_an_ranks_GM[order(chromosome_an_ranks[,1]),1],chromosome_an_ranks_GM_se[order(chromosome_an_ranks[,1]),1],
                  chromosome_an_ranks_SA[order(chromosome_an_ranks[,1]),1],chromosome_an_ranks_SA_se[order(chromosome_an_ranks[,1]),1],
                  chromosome_an_ranks_MC[order(chromosome_an_ranks[,1]),1],chromosome_an_ranks_MC_se[order(chromosome_an_ranks[,1]),1],
                  chromosome_an_ranks_IC[order(chromosome_an_ranks[,1]),1],chromosome_an_ranks_IC_se[order(chromosome_an_ranks[,1]),1])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)==""),],table(chromosome_labels)[-c(1,24)][tmp1])
tmp[,11] <- as.factor(tmp1)
tmp[,seq(1,10,2)] <- (-tmp[,seq(1,10,2)]+15043)-7521.5
colnames(tmp) <- c("Rank_CT","SE_CT","Rank_GM","SE_GM","Rank_SA","SE_SA","Rank_MC","SE_MC","Rank_IC","SE_IC","Chromosome","N")
write.csv(tmp, "SCAx_chromosome_ranks_PLS_features.csv",row.names=F)

tmp <- data.frame(chromosome_an_ranks_CT[order(chromosome_an_ranks[,2]),2],chromosome_an_ranks_CT_se[order(chromosome_an_ranks[,2]),2],
                  chromosome_an_ranks_GM[order(chromosome_an_ranks[,2]),2],chromosome_an_ranks_GM_se[order(chromosome_an_ranks[,2]),2],
                  chromosome_an_ranks_SA[order(chromosome_an_ranks[,2]),2],chromosome_an_ranks_SA_se[order(chromosome_an_ranks[,2]),2],
                  chromosome_an_ranks_MC[order(chromosome_an_ranks[,2]),2],chromosome_an_ranks_MC_se[order(chromosome_an_ranks[,2]),2],
                  chromosome_an_ranks_IC[order(chromosome_an_ranks[,2]),2],chromosome_an_ranks_IC_se[order(chromosome_an_ranks[,2]),2])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)==""),],table(chromosome_labels)[-c(1,24)][tmp1])
tmp[,11] <- as.factor(tmp1)
tmp[,seq(1,10,2)] <- (-tmp[,seq(1,10,2)]+15043)-7521.5
colnames(tmp) <- c("Rank_CT","SE_CT","Rank_GM","SE_GM","Rank_SA","SE_SA","Rank_MC","SE_MC","Rank_IC","SE_IC","Chromosome","N")
write.csv(tmp, "SCAy_chromosome_ranks_PLS_features.csv",row.names=F)


tmp <- data.frame(chromosome_Turner_ranks_CT[rownames(chromosome_Turner_ranks)[order(chromosome_Turner_ranks[,1])],],
                  chromosome_Turner_ranks_GM[rownames(chromosome_Turner_ranks)[order(chromosome_Turner_ranks[,1])],],
                  chromosome_Turner_ranks_SA[rownames(chromosome_Turner_ranks)[order(chromosome_Turner_ranks[,1])],],
                  chromosome_Turner_ranks_MC[rownames(chromosome_Turner_ranks)[order(chromosome_Turner_ranks[,1])],],
                  chromosome_Turner_ranks_IC[rownames(chromosome_Turner_ranks)[order(chromosome_Turner_ranks[,1])],])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="NA")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="NA"),],table(chromosome_labels)[-c(1,24)][tmp1])
tmp[,11] <- as.factor(tmp1)
tmp[,seq(1,10,2)] <- (-tmp[,seq(1,10,2)]+15043)-7521.5
colnames(tmp) <- c("Rank_CT","SE_CT","Rank_GM","SE_GM","Rank_SA","SE_SA","Rank_MC","SE_MC","Rank_IC","SE_IC","Chromosome","N")
write.csv(tmp, "Turner_chromosome_ranks_PLS_features.csv",row.names=F)

tmp <- data.frame(chromosome_Downs_ranks_CT[rownames(chromosome_Downs_ranks)[order(chromosome_Downs_ranks[,1])],],
                  chromosome_Downs_ranks_GM[rownames(chromosome_Downs_ranks)[order(chromosome_Downs_ranks[,1])],],
                  chromosome_Downs_ranks_SA[rownames(chromosome_Downs_ranks)[order(chromosome_Downs_ranks[,1])],],
                  chromosome_Downs_ranks_MC[rownames(chromosome_Downs_ranks)[order(chromosome_Downs_ranks[,1])],],
                  chromosome_Downs_ranks_IC[rownames(chromosome_Downs_ranks)[order(chromosome_Downs_ranks[,1])],])
tmp1 <- rownames(tmp)
tmp1 <- tmp1[-which(rownames(tmp)=="Un"|rownames(tmp)=="NA")]
tmp <- data.frame(tmp[-which(rownames(tmp)=="Un"|rownames(tmp)=="NA"),],table(chromosome_labels)[-c(1,24)][tmp1])
tmp[,11] <- as.factor(tmp1)
tmp[,seq(1,10,2)] <- (-tmp[,seq(1,10,2)]+15043)-7521.5
colnames(tmp) <- c("Rank_CT","SE_CT","Rank_GM","SE_GM","Rank_SA","SE_SA","Rank_MC","SE_MC","Rank_IC","SE_IC","Chromosome","N")
write.csv(tmp, "Downs_chromosome_ranks_PLS_features.csv",row.names=F)





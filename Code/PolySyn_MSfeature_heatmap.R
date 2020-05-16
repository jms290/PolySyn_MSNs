# For NatComms rev1

# Housekeeping

setwd("~/Desktop/NatComms_revisions/")

library(ggplot2)
library(reshape2)

# Load in single feature contrast maps (not necessary but good to store with these other data)

SCAx_morph_nodes_effects_lh <- data.frame(CT=read.table("SCAx_CT_DXvNV.txt")[1:152,],
                                           GM=read.table("SCAx_GM_DXvNV.txt")[1:152,],
                                           SA=read.table("SCAx_SA_DXvNV.txt")[1:152,],
                                           MC=read.table("SCAx_MC_DXvNV.txt")[1:152,],
                                           IC=read.table("SCAx_IC_DXvNV.txt")[1:152,])
SCAy_morph_nodes_effects_lh <- data.frame(CT=read.table("SCAy_CT_DXvNV.txt")[1:152,],
                                           GM=read.table("SCAy_GM_DXvNV.txt")[1:152,],
                                           SA=read.table("SCAy_SA_DXvNV.txt")[1:152,],
                                           MC=read.table("SCAy_MC_DXvNV.txt")[1:152,],
                                           IC=read.table("SCAy_IC_DXvNV.txt")[1:152,])
Downs_morph_nodes_effects_lh <- data.frame(CT=read.table("Downs_CT_DXvNV.txt")[1:152,],
                                           GM=read.table("Downs_GM_DXvNV.txt")[1:152,],
                                           SA=read.table("Downs_SA_DXvNV.txt")[1:152,],
                                           MC=read.table("Downs_MC_DXvNV.txt")[1:152,],
                                           IC=read.table("Downs_IC_DXvNV.txt")[1:152,])
Turner_morph_nodes_effects_lh <- data.frame(CT=read.table("Turner_CT_DXvNV.txt")[1:152,],
                                           GM=read.table("Turner_GM_DXvNV.txt")[1:152,],
                                           SA=read.table("Turner_SA_DXvNV.txt")[1:152,],
                                           MC=read.table("Turner_MC_DXvNV.txt")[1:152,],
                                           IC=read.table("Turner_IC_DXvNV.txt")[1:152,])
VELOC_morph_nodes_effects_lh <- data.frame(CT=read.table("VELOC_CT_DXvNV.txt")[1:152,],
                                           GM=read.table("VELOC_GM_DXvNV.txt")[1:152,],
                                           SA=read.table("VELOC_SA_DXvNV.txt")[1:152,],
                                           MC=read.table("VELOC_MC_DXvNV.txt")[1:152,],
                                           IC=read.table("VELOC_IC_DXvNV.txt")[1:152,])
WAGR_morph_nodes_effects_lh <- data.frame(CT=read.table("WAGR_CT_DXvNV.txt")[1:152,],
                                           GM=read.table("WAGR_GM_DXvNV.txt")[1:152,],
                                           SA=read.table("WAGR_SA_DXvNV.txt")[1:152,],
                                           MC=read.table("WAGR_MC_DXvNV.txt")[1:152,],
                                           IC=read.table("WAGR_IC_DXvNV.txt")[1:152,])

# read in tables if leave-one-feature-out (LOFO) MS change maps

SCAx_LOFO <- read.table("SCAx_MSNs_LOFO.txt")
SCAy_LOFO <- read.table("SCAy_MSNs_LOFO.txt")
Downs_LOFO <- read.table("Downs_MSNs_LOFO.txt")
Turner_LOFO <- read.table("Turner_MSNs_LOFO.txt")
VELOC_LOFO <- read.table("VELOC_MSNs_LOFO.txt")
WAGR_LOFO <- read.table("WAGR_MSNs_LOFO.txt")

# make heatmap of correlations between LOFO maps and empirical MS change map

hm <- data.frame("+X"=cor(SCAx_LOFO)[1,],"-X"=cor(Turner_LOFO)[1,],
                 "+Y"=cor(SCAy_LOFO)[1,],"+21"=cor(Downs_LOFO)[1,],
           "-22q11"=cor(VELOC_LOFO)[1,],"-11p13"=cor(WAGR_LOFO)[1,])
colnames(hm) <- c("+X","-X","+Y","+21","-22q11","-11p13")
hm <- melt(t(hm))
colnames(hm) <- c("Group","Feature","r")
hm_mins <- rep("",36)
hm_mins[hm$r %in% sapply(levels(hm$Group),function(x) min(hm$r[which(hm$Group==x)]))] <- "#"
ggplot(hm,aes(Group,Feature,fill=r))+geom_tile()+geom_text(label=hm_mins,size=7)+theme_minimal()+
  theme(text=element_text(size=24))+scale_fill_gradient2(low="blue",mid="white",high="red")+
  ylab("Feature (left out)")


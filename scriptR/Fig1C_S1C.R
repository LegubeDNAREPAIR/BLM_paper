source("src/Functions.R")
library("BSgenome.Hsapiens.UCSC.hg19") ;
library("tidyverse")
seqlens <- seqlengths( BSgenome.Hsapiens.UCSC.hg19 );
library( "rtracklayer" ) ;
library( "plyranges" ) ;
#Figure gène les plus proches des top 20
load("../data/asi.closed.gene.rda")
load("../data/one.closed.genes.rda")
bless80 <- read_bed("../data/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
res.bless <- "../data/BLESS80_ordering_by_ratioBLM2kb_BLESS500bp.csv" %>% read_csv()
ens.to.seqname <- names(seqlens)[1:24]
newseqlens <- seqlens[names(seqlens) %in% ens.to.seqname]
names(newseqlens) <- names(ens.to.seqname)

res.plot <- res.bless %>% select(Top_ratio,ratio.BLM.BLESS,ratio.RAD51.XRCC4) %>% filter(Top_ratio %in% c("BLM-high","BLM-low"))%>%  reshape2::melt()
##FIG S1C
p <- ggplot(res.plot,aes(x=Top_ratio,y=value)) + stat_boxplot(geom ='errorbar') + geom_boxplot() + facet_wrap(~variable,scales = "free_y")+ theme_classic() + labs(y='BLM/BLESS ratio',x='Categories',title = "")
pdf("Ratio_BLMvsBLESS.pdf",height=12,width=20)
print(p)
dev.off()


blm.high.BLMBLESS <- res.plot %>% filter(Top_ratio == "BLM-high" & variable == "ratio.BLM.BLESS") %>% select(value)
blm.low.BLMBLESS <- res.plot %>% filter(Top_ratio == "BLM-low" & variable == "ratio.BLM.BLESS") %>% select(value)

blm.high.RAD51 <- res.plot %>% filter(Top_ratio == "BLM-high" & variable == "ratio.RAD51.XRCC4") %>% select(value)
blm.low.RAD51 <- res.plot %>% filter(Top_ratio == "BLM-low" & variable == "ratio.RAD51.XRCC4") %>% select(value)

wilcox.test(blm.high.BLMBLESS$value,blm.low.BLMBLESS$value)$p.value
wilcox.test(blm.high.RAD51$value,blm.low.RAD51$value)$p.value

#Profile TOP High/LOW BLM for BLM
BLM.high.sites <- res.bless[which(res.bless$Top_ratio == "BLM-high"),]$SITES
BLM.high.profile <- computeProfile(bed = bless80[which(bless80$name %in% BLM.high.sites)],wig = myWigBLM$pOHT,w = 10000,span=200,seqlens = seqlens)
BLM.low.sites <- res.bless[which(res.bless$Top_ratio == "BLM-low"),]$SITES
BLM.low.profile <- computeProfile(bed = bless80[which(bless80$name %in% BLM.low.sites)],wig = myWigBLM$pOHT,w = 10000,span=200,seqlens = seqlens)


#Cat POLII

bless80.gene.old <- bless80[which(bless80$name %in% names(asi.closed.genes))]
bless80.genes <- one.closed.genes[which(names(one.closed.genes) %in% bless80.gene.old$name)]
bless80.n <- unlist(bless80.genes) %>% names() %>% lapply(.,function(x){unlist(strsplit(x,"\\."))[1]}) %>% unlist()
bless80.genes <- unlist(bless80.genes)
bless80.genes$SITES <- bless80.n

BLM.high.genes <- bless80.genes[which(bless80.genes$SITES %in% BLM.high.sites)]
BLM.low.genes <- bless80.genes[which(bless80.genes$SITES %in% BLM.low.sites)]

files <- c("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Ser2P_Pol2/PROCESSED/BY_US/WIGGLE/H3LNWBGX3_Ser2P_Pol2_DIVA_17s001630-1-1_Clouaire_lane117s001630_sequence_normalized.bw","/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Tot-Pol2/PROCESSED/BY_US/WIGGLE/H3LNWBGX3_Tot_Pol2_DIVA_17s001632-1-1_Clouaire_lane117s001632_sequence_normalized.bw")
names(files) <- c("Ser2P_Pol2_DIVA","Tot_Pol2_DIVA")

aparcourir <- c(files,"H3K36me3"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K36me3/BGI_RUN1_201203/PROCESSED/BIGWIG/h3k36me3_mOHT_normalized_hg19.bw")
aparcourir <- aparcourir %>% map(import.bw,as="RleList")
monjoliprof <- list()
for(f in names(aparcourir)){
  wig <- aparcourir[[f]]
  bed <- BLM.high.genes
  
  monjoliprof[[f]][["high"]][["genebody"]] <- computeProfileGenesENS(bed,wig,seqlens=seqlens,fun="mean",get="all")
  monjoliprof[[f]][["high"]][["TSS"]] <- computeProfileTSSENS(bed,wig,seqlens=seqlens,w1=5000,w2=10000,span=200,fun="mean",get="all")
  bed <- BLM.low.genes
  
  monjoliprof[[f]][["low"]][["genebody"]] <- computeProfileGenesENS(bed,wig,seqlens=seqlens,fun="mean",get="all")
  monjoliprof[[f]][["low"]][["TSS"]] <- computeProfileTSSENS(bed,wig,seqlens=seqlens,w1=5000,w2=10000,span=200,fun="mean",get="all")
  
}
##GENEBODY
data = rbind(data.frame("Windows"=1:130,"value"=colMeans(monjoliprof[["Ser2P_Pol2_DIVA"]][["high"]][["genebody"]]),"Type"="BLM-High",Histone = "Ser2P_Pol2_DIVA"),
             data.frame("Windows"=1:130,"value"=colMeans(monjoliprof[["Ser2P_Pol2_DIVA"]][["low"]][["genebody"]]),"Type"="BLM-Low",Histone = "Ser2P_Pol2_DIVA"),
             data.frame("Windows"=1:130,"value"=colMeans(monjoliprof[["Tot_Pol2_DIVA"]][["high"]][["genebody"]]),"Type"="BLM-High",Histone = "Tot_Pol2_DIVA"),
             data.frame("Windows"=1:130,"value"=colMeans(monjoliprof[["Tot_Pol2_DIVA"]][["low"]][["genebody"]]),"Type"="BLM-Low",Histone = "Tot_Pol2_DIVA"),
             data.frame("Windows"=1:130,"value"=colMeans(monjoliprof[["H3K36me3"]][["high"]][["genebody"]]),"Type"="BLM-High",Histone = "H3K36me3"),
             data.frame("Windows"=1:130,"value"=colMeans(monjoliprof[["H3K36me3"]][["low"]][["genebody"]]),"Type"="BLM-Low",Histone = "H3K36me3"))


for(f in unique(data$Histone)){
  dat.plot <- data %>% filter(Histone == f)
  p<- ggplot(dat.plot,aes(Windows,value,colour= Type)) +
    geom_line() +
    labs(list(title = "", x = "Position", y = "Average CHIP-seq count")) +
    #geom_ribbon(aes(ymin=0, ymax=value,fill="#e74c3c"),alpha= 0.5) +
    geom_vline(xintercept = 15,linetype="longdash") +
    geom_vline(xintercept = 115,linetype="longdash") +
    scale_x_continuous(name = 'Position',
                       breaks = c(1,15,115,130),
                       labels = c("TSS-3000bp", 'TSS', 'TES', 'TES+3000bp')
    ) +
    #facet_wrap(~Type,ncol=2) +
    theme_classic()
  pdf(paste("Profile_on_Gene",f,"BLM-High_BLM_Low.pdf",sep="_"),height=12,width=18)
  print(p)
  dev.off()
}
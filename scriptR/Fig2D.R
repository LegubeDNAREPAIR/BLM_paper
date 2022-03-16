#Scatterplot
source("src/Functions.R")
require(tidyverse)
require(plyranges)
require(rtracklayer)
seqlens <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 %>% seqlengths()
bless80 <- read_bed("../data/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
bw.files <- c("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/XRCC4_RAD51/EMBL2_64HK1AAXX_RAD51_XRCC4_4h/PROCESSED/mapping/BIGWIG/64HK1AAXX_RAD51_12s004588-1-1_Legube_lane812s004588_sequence_normalized.bw",
              "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/XRCC4_RAD51/EMBL3_HJKVVBGX5_RAD51_XRCC4_1h_4h_24h/PROCESSED/mapping/BIGWIG/HJKVVBGX5_RAD51_24H_17s006149-1-1_Clouaire_lane117s006149_sequence_normalized.bw",
              "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/BLM/BLM_4H_BGI_RUN1_201203/4H/PROCESSED/ALIGNED/WIGGLE/BLM_normalized_01_02_2018.bw",
              "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/BLM/24H/HKK2NBGX7/PROCESSED/BIGWIG/HKK2NBGX7_BLM_24H_18s002882-1-1_Clouaire_lane118s002882_sequence_normalized.bw"
              
)
Cond <- c( paste("RAD51",c("4H","24H"),sep="_"), paste("BLM",c("4H","24H"),sep="_"))


wigs <- lapply(bw.files,import.bw,as = "RleList")
names(wigs) <- Cond

#V4
df.scat <- do.call("cbind",lapply(names(wigs),function(z){
  wig <- wigs[[z]]
  do.call("cbind",lapply(c(5000,20000),function(w){
    compute1ValPerSite(bless80,wig = wig,w = w,seqlens = seqlens)
  }))
}))
colnames(df.scat) <- paste(rep(names(wigs),each=2),rep(c(5000,20000),length(names(wigs))),sep="_")



pdf("Correlation_BLM_High_Low_RAD51_BLM_4H_24H.pdf",height=12,width=12)
df.scat %>% as.data.frame %>% ggplot(aes(x=RAD51_4H_5000,y=BLM_4H_5000)) + geom_point() + theme_classic() + geom_smooth(method="lm")
df.scat %>% as.data.frame %>% ggplot(aes(x=RAD51_24H_20000,y=BLM_24H_20000)) + geom_point() + theme_classic() + geom_smooth(method="lm")


matrix.count.cors <- cor(df.scat,method="spearman")
corrplot(  matrix.count.cors, type = "upper", method = "circle", order = "original", addCoef.col = "black", diag = FALSE, 
           tl.col = "black", tl.cex = 0.8, tl.srt = 20, title = "",number.cex = 0.8)

dev.off()
library(Rsamtools) ; # ScanBam
library(GenomicAlignments) ; # GAlignement
library("BSgenome.Hsapiens.UCSC.hg19") ;
library( "rtracklayer" ) ; #Export bw
source("src/Functions.R")
seqlens <- seqlengths( Hsapiens );

asi = import.bed("../data/ASIsites_hg19.bed")
bless80 <- import.bed("../data/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
#S1A
files <- "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/BLM/BLM_4H_BGI_RUN1_201203/4H/PROCESSED/ALIGNED/WIGGLE/BLM_normalized_01_02_2018.bw"

myWigBLM <- readData( dir = "", fnpOHT = files, fnmOHT = NULL );

files <- "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/GAMMA/U2OS/BGI_RUN1_201203/PROCESSED/ALIGNED/WIGGLE/hg19/GAM.clean_normalized_01022018.bw"

myWigGAMMA <- readData( dir = "", fnpOHT = files, fnmOHT = NULL );

w <- 50000
span <- 1000
BLM.prof <- computeProfile(bless80,wig = myWigBLM$pOHT,w = 50000,span=1000,seqlens = seqlens)
GAM.prof <- computeProfile(bless80,wig = myWigGAMMA$pOHT,w = 50000,span=1000,seqlens = seqlens)
data = rbind(data.frame("Windows"=seq(-w,w-span+1,span),"value"=colMeans(BLM.prof),"Type"="bless",Histone = "BLM"),
             data.frame("Windows"=seq(-w,w-span+1,span),"value"=colMeans(GAM.prof),"Type"="bless",Histone = "GAMMA"))
p <- ggplot(na.omit(data),aes(Windows,value,colour = Histone)) +
  labs(list(title = "", x = "", y = "")) +
  geom_line()+
  #facet_wrap(~ Histone,ncol = 3,scales = "free_x") +
  #scale_colour_manual(values=whichcolor) +
  theme_classic()

pdf("FIG_prof_BLM_GAMMA_50kb_span1000_pOHT.pdf",height=8,width=14)
print(p)
dev.off()

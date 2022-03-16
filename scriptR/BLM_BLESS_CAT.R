source("src/Functions.R")
library("BSgenome.Hsapiens.UCSC.hg19") ;
library("tidyverse")
seqlens <- seqlengths( BSgenome.Hsapiens.UCSC.hg19 );
library( "rtracklayer" ) ;
library( "plyranges" ) ;
#Figure gène les plus proches des top 20

BLESS.dir <- c("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/BLESS/CHIP-SEQ/Jun2015/PROCESSED/ALIGNED_PAIRED/WIGGLE/trimmed_BLESS_U2OSpTam.sorted_fragments_800bp_rmdups_normalized.bw",
               "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/BLESS/CHIP-SEQ/Jun2015/PROCESSED/ALIGNED_PAIRED/WIGGLE/trimmed_BLESS_U2OS-Tam.sorted_fragments_800bp_rmdups_normalized.bw")
myWigBLESS = readData( dir = "", fnpOHT = BLESS.dir[1], fnmOHT = BLESS.dir[2] );

w = 500
span = 1
profs = list()
profs[["BLESS_80"]][["pOHT"]] <- compute1ValPerSite(bed = bless80, wig = myWigBLESS[["pOHT"]], w=w, seqlens = seqlens )
#profs[["BLESS_80"]][["mOHT"]] <- computeProfile(bed = bless80, wig = myWigBLESS[["mOHT"]], w=w, span=span, seqlens = oldseqlens , method="mean")


myWigBLM <- readData( dir = "", fnpOHT = "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/BLM/BLM_4H_BGI_RUN1_201203/4H/PROCESSED/ALIGNED/WIGGLE/BLM_normalized_01_02_2018.bw", fnmOHT = NULL );

w = 2000
profs[["BLM"]][["pOHT"]] <-  compute1ValPerSite(bed = bless80, wig = myWigBLM[["pOHT"]], w=w, seqlens = seqlens )

#ADD RAD51 and XRCC4
w = 5000

wXRCC4 <- 1000
fnmp = list()
# fnmp[["RAD51"]]=c("/mnt/volumes/NAS1/DATA/ChIP-Seq/RAD51/U2OS/EMBL2/PROCESSED/WIGGLE/hg19/RAD51_hg19_normalized.bw")
fnmp[["RAD51"]]=c("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/XRCC4_RAD51/EMBL2_64HK1AAXX_RAD51_XRCC4_4h/PROCESSED/OLD PROCESS/WIGGLE/hg19/RAD51_hg19_normalized.bw")
#ADD XRCC4 location
fnmp[["XRCC4"]]=c("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/XRCC4_RAD51/EMBL2_64HK1AAXX_RAD51_XRCC4_4h/PROCESSED/OLD PROCESS/WIGGLE/hg19/XRCC4_hg19_normalized.bw")


wigRAD51 = readData( dir = "", fnpOHT = fnmp[["RAD51"]][1], fnmOHT = NULL );
profRAD51 = compute1ValPerSite(bed = bless80, wig = wigRAD51[[1]], w = w, seqlens = seqlens )


wigXRCC4 = readData( dir = "", fnpOHT = fnmp[["XRCC4"]][1], fnmOHT = NULL );
profXRCC4 = compute1ValPerSite(bed = bless80, wig = wigXRCC4[[1]], w = wXRCC4,seqlens = seqlens)


res.bless <- data.frame("SITES"=bless80$name,"BLM"=profs[["BLM"]][["pOHT"]],"BLESS"=profs[["BLESS_80"]][["pOHT"]],"RAD51"=profRAD51,"XRCC4"=profXRCC4)
res.bless <- res.bless %>% mutate(ratio.BLM.BLESS = BLM/BLESS) %>% mutate(ratio.RAD51.XRCC4 = RAD51/XRCC4) %>% arrange(desc(ratio.BLM.BLESS))

res.bless$Top_ratio= "NA"
res.bless[1:20,]$Top_ratio = "BLM-high"
res.bless[61:80,]$Top_ratio = "BLM-low"

write_csv(res.bless,"../data/BLESS80_ordering_by_ratioBLM2kb_BLESS500bp.csv")

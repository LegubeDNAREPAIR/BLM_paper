require(tidyverse)
require(cowplot)
require(plyranges)
require(rtracklayer)
require(parallel)
require(gridExtra)
require(grid)
require(BSgenome.Hsapiens.UCSC.hg19.masked)

bin.var <- function (x, bins = 4, method = c("intervals", "proportions", "natural"), labels = FALSE) {
  method <- match.arg(method)
  if (length(x) < bins) {
    stop()
  }
  x <- if (method == "intervals")
    cut(x, bins, labels = labels)
  else if (method == "proportions")
    cut(x, quantile(x, probs = seq(0, 1, 1/bins), na.rm = TRUE), include.lowest = TRUE, labels = labels)
  else {
    xx <- na.omit(x)
    breaks <- c(min(xx), tapply(xx, KMeans(xx, bins)$cluster,
                                max))
    cut(x, breaks, include.lowest = TRUE, labels = labels)
  }
  as.factor(x)
}

seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
bless80 = read_bed("../data/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")

files <- c(
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Tot-Pol2/PROCESSED/BY_US/WIGGLE/H3LNWBGX3_Tot_Pol2_DIVA_17s001632-1-1_Clouaire_lane117s001632_sequence_normalized.bw",
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Ser2P_Pol2/PROCESSED/BY_US/WIGGLE/H3LNWBGX3_Ser2P_Pol2_DIVA_17s001630-1-1_Clouaire_lane117s001630_sequence_normalized.bw"
)


mes_wig <- c("Tot_Pol2_mOHT","Ser2P_Pol2_mOHT")
names(files) <- mes_wig
files <- files %>% purrr::map(import.bw,as="RleList")  %>% setNames(mes_wig)

mes_windows <- c(5000)
boxplots.dat <- lapply(mes_windows,function(one_win){
  lapply(names(files),function(one.w){
    message(one.w)
    # PhDfunc::Get1val(Name=one.w,G4.files[[one.w]],bless80 %>%
    #                    anchor_center() %>%
    #                    mutate(width=one_win))
    tibble(rowname = bless80$name,
           value = PhDfunc::compute1ValperPeak(bless80,files[[one.w]],w=one_win,seqlens = seqlens),
           wig = one.w)
  }) %>% bind_rows()
}) %>% setNames(mes_windows)%>% bind_rows(.id = "Type")


myWig <- import.bw("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/BLM/BLM_4H_BGI_RUN1_201203/4H/PROCESSED/ALIGNED/WIGGLE/BLM_normalized_01_02_2018.bw",as="RleList")
BLM.prof.2 <- PhDfunc::compute1ValperPeak(bless80,myWig,w=2000,seqlens = seqlens)

final.dat.80 <- tibble(rowname = bless80$name,
                       value = BLM.prof.2)


GetCatPOLII <- function(mydat){
  cat.name <- c("Low","Mid1","Mid2","High")
  mydat  %>% mutate(Cat = bin.var(x = mydat$value,bins = 4,method = "proportions",labels = cat.name))
}

boxplots.dat <- boxplots.dat %>% group_by(wig) %>% nest() %>% mutate(data = purrr::map(data,GetCatPOLII)) %>% unnest(data)

boxplots.dat <- boxplots.dat %>% left_join(final.dat.80,by="rowname",suffix = c("_PolII", "_BLM"))


p <- boxplots.dat %>% mutate(Cat = forcats::fct_relevel(Cat,c("High","Mid2","Mid1","Low"))) %>% ggplot(aes(x=Cat,y=value_BLM)) + stat_boxplot(geom ='errorbar') + geom_boxplot() + theme_classic() +facet_wrap(~wig,ncol=1)


sub_dat <- boxplots.dat %>% spread(key = Cat,value = value_BLM) %>% group_by(wig) %>% nest()

res <- lapply(sub_dat$data,function(my_dat){
  my_cat <- list(
    Low = na.omit(my_dat$Low) %>% as.vector(),
    Mid2 = na.omit(my_dat$Mid2) %>% as.vector(),
    Mid1 = na.omit(my_dat$Mid1) %>% as.vector(),
    High = na.omit(my_dat$High) %>% as.vector()
  )
  
  mes_test <- combn(names(my_cat),2)
  lapply(1:ncol(mes_test),function(i){
    Cat1 <- mes_test[1,i]
    Cat2 <- mes_test[2,i]
    
    tibble(Cat1=Cat1,Cat2=Cat2,p.value = wilcox.test(my_cat[[Cat1]],my_cat[[Cat2]])$p.value)
    
    
  }) %>% bind_rows()
  
}) %>% setNames(sub_dat$wig) %>% bind_rows(.id="wig")

pdf("PLOTS/Boxplot_2kb_polCat_BLM.pdf",height=12,width=12)
print(p)
grid.newpage()
grid.table(res %>% arrange(p.value))
dev.off()


require(tidyverse)
require(cowplot)
require(plyranges)
require(rtracklayer)
require(parallel)
require(gridExtra)
require(BSgenome.Hsapiens.UCSC.hg19.masked)
seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
bless80 = read_bed("../data/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
res.bless <- "../data/BLESS80_ordering_by_ratioBLM2kb_BLESS500bp.csv" %>% read_csv()
res.bless$Top_ratio= "NA"
res.bless[1:20,]$Top_ratio = "BLM-high"
res.bless[61:80,]$Top_ratio = "BLM-low"
#PROFIL G4 LOOP
G4.files <- list.files(path = "/mnt/NAS/DATA/DATA_FROM_OTHER_LABS_and_DATABASES/G4/G4_GSE99205/PROCESSED/mapping/BIGWIG",pattern = "_normalized.bw",full.names = T) %>% 
  setNames(str_remove(basename(.),"_normalized.bw")) %>% purrr::map(import.bw,as="RleList")
BLM.high.sites <- res.bless[which(res.bless$Top_ratio == "BLM-high"),]$SITES
BLM.low.sites <- res.bless[which(res.bless$Top_ratio == "BLM-low"),]$SITES
mes_windows <- c(500,1000,2000,5000,10000,20000,50000)
data.G4 <- mclapply(mes_windows,function(one_win){
  lapply(names(G4.files),function(one.w){
    message(one.w)
    # PhDfunc::Get1val(Name=one.w,G4.files[[one.w]],bless80 %>%
    #                    anchor_center() %>%
    #                    mutate(width=one_win))
    tibble(rowname = bless80$name,
           value = PhDfunc::compute1ValperPeak(bless80,G4.files[[one.w]],w=one_win,seqlens = seqlens),
           wig = one.w)
  }) %>% bind_rows()
},mc.cores=length(mes_windows)) %>% setNames(mes_windows)%>% bind_rows(.id = "Type")




data.G4 <- data.G4 %>% mutate(Group = case_when(
  rowname %in% BLM.high.sites ~ "BLM-high",
  rowname %in% BLM.low.sites ~ "BLM-low",
  TRUE ~ "None"
)) %>% filter(wig != "SRR5586990")
p <- data.G4 %>% filter(Group != "None") %>% ggplot(aes(x=wig,y=value,fill=Group)) + stat_boxplot(geom ='errorbar') + geom_boxplot() + facet_wrap(~Type,ncol=1,scales="free_y") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
my.test <- . %>%
  wilcox.test(value~Group,data=.) %>%  .$p.value
table_moi_ca <- data.G4 %>%
  filter(Group != "None") %>%
  group_by(Type,wig) %>% nest() %>% mutate(wilcox.test.df = purrr::map_dbl(data,my.test)) %>% 
  dplyr::select(Type,wig,wilcox.test.df)


p <- data.G4 %>% filter(Group != "None",Type =="1000") %>% ggplot(aes(x=wig,y=value,fill=Group)) + stat_boxplot(geom ='errorbar') + geom_boxplot() + facet_wrap(~Type,ncol=1,scales="free_y") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("PLOTS/FIG_boxplot_superpose_G4_BLM-High_BLM-Low_pOHT_1kb_GSE99205.pdf",height=6,width=12)
print(p)

dev.off()






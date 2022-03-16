##FOR FIG 3D 3E

require(dplyr)
require(rtracklayer)
require(ggplot2)
require(plyranges)
library("BSgenome.Hsapiens.UCSC.hg19")
seqlens = seqlengths( Hsapiens );
#Load AsiSI locations
asi = import.bed("../data/ASIsites_hg19.bed")
bless80 = import.bed("../data/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
HR = import.bed("../data/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed")
NHEJ = import.bed("../data/BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed")
Random80 = import.bed("../data/80random.bed")
Random30 = import.bed("../data/30random.bed")
uncut <- asi[!asi$name %in% bless80$name]

list.sites <- list("cut"=bless80,"uncut"=uncut,"HR"=HR,"NHEJ"=NHEJ,"Random80"=Random80,"Random30"=Random30)

windows = list("500kb"=500000,"20kb"=20000,"5kb"=5000,"500bp"=500)
spans = list("500kb"=10000,"20kb"=200,"5kb"=50,"500bp"=1)


wig.List <- list(
  "HKKWHBGX7_DRIP_C1"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HKKWHBGX7/PROCESSED/mapping/BIGWIG/HKKWHBGX7_DRIP_C1_18s002883-1-1_Clouaire_lane118s002883_sequence_normalized.bw",
  "HKKWHBGX7_DRIP_C1_24H"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HKKWHBGX7/PROCESSED/mapping/BIGWIG/HKKWHBGX7_DRIP_C1_24H_18s002884-1-1_Clouaire_lane118s002884_sequence_normalized.bw",
  "HLFCGBGX7_DRIP_BLM"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HLFCGBGX7/PROCESSED/mapping/BIGWIG/HLFCGBGX7_DRIP_BLM_18s002887-1-1_Clouaire_lane118s002887_sequence_normalized.bw",
  "HLFCGBGX7_DRIP_BLM_24H"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HLFCGBGX7/PROCESSED/mapping/BIGWIG/HLFCGBGX7_DRIP_BLM_24H_18s002888-1-1_Clouaire_lane118s002888_sequence_normalized.bw",
  "HKKWHBGX7_DRIP_STX"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HKKWHBGX7/PROCESSED/mapping/BIGWIG/HKKWHBGX7_DRIP_STX_18s002885-1-1_Clouaire_lane118s002885_sequence_normalized.bw",
  "HKKWHBGX7_DRIP_STX_24H"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HKKWHBGX7/PROCESSED/mapping/BIGWIG/HKKWHBGX7_DRIP_STX_24H_18s002886-1-1_Clouaire_lane118s002886_sequence_normalized.bw",
  "HLFCGBGX7_DRIP_BLMSTX"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HLFCGBGX7/PROCESSED/mapping/BIGWIG/HLFCGBGX7_DRIP_BLMSTX_18s002889-1-1_Clouaire_lane118s002889_sequence_normalized.bw",
  "HLFCGBGX7_DRIP_BLMSTX_24H"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HLFCGBGX7/PROCESSED/mapping/BIGWIG/HLFCGBGX7_DRIP_BLMSTX_24H_18s002890-1-1_Clouaire_lane118s002890_sequence_normalized.bw",
  
  
  "DRIP_BLM_vs_CTRL_mOHT"="/media/HDD_ROCHER/PROJET_INGE/DRIP-SEQ_BAMCOMPARE/DRIP_BLM_vs_CTRL_mOHT.bw",
  "DRIP_BLM_vs_CTRL_pOHT24H"="/media/HDD_ROCHER/PROJET_INGE/DRIP-SEQ_BAMCOMPARE/DRIP_BLM_vs_CTRL_pOHT24H.bw",
  "DRIP_STX_vs_CTRL_mOHT_HKKWHBGX7"="/media/HDD_ROCHER/PROJET_INGE/DRIP-SEQ_BAMCOMPARE/DRIP_STX_vs_CTRL_mOHT_HKKWHBGX7.bw",
  "DRIP_STX_vs_CTRL_HKKWHBGX7"="/media/HDD_ROCHER/PROJET_INGE/DRIP-SEQ_BAMCOMPARE/DRIP_STX_vs_CTRL_HKKWHBGX7.bw",
  "DRIP_BLMSTX_vs_CTRL_mOHT"="/media/HDD_ROCHER/PROJET_INGE/DRIP-SEQ_BAMCOMPARE/DRIP_BLMSTX_vs_CTRL_mOHT.bw",
  "DRIP_BLMSTX_vs_CTRL_pOHT24H"="/media/HDD_ROCHER/PROJET_INGE/DRIP-SEQ_BAMCOMPARE/DRIP_BLMSTX_vs_CTRL_pOHT24H.bw"
)

wigs <- lapply(wig.List,import.bw,as="RleList") %>% setNames(names(wig.List))
#process BW
computeProfile = function( bed, wig, w = 20000, span = 200, seqlens ,method="mean"){
  if( class( wig ) != "SimpleRleList" ){
    stop( "ERROR : unknown class of wig, please provide a SimpleRleList" );
  }
  mat = NULL;
  for( i in 1:length( bed ) ){
    message( i, "/", length( bed ) );
    bedi = bed[i, ];
    chr = as.character( seqnames( bedi ) );
    cov = wig[[chr]];
    
    
    center = start( bedi ) + 4;
    stW = center - w;
    edW = center + w;
    
    if( span == 1 ){
      vm = as.numeric( Views( cov, start = stW, end = edW )[[1]] )
    }else{
      sts = seq( stW, edW - span + 1, span );
      eds = seq( stW + span - 1, edW, span );
      v = Views( cov, start = sts, end = eds );
      if(method =="sum"){
        vm =  sum( v );    
      }else {
        vm =  mean( v );
      }
      
      vm[sts >= seqlens[chr] | sts > length( cov )] = 0;
      vm[sts < 1] = 0;
    }
    mat = rbind( mat, vm );
  }
  #rv = colMeans( mat );
  return( mat );
}



dat.boxplot.all <- lapply(names(windows),function(one_w){
  lapply(wigs,function(wig){
    lapply(names(list.sites),function(one_n){
      x <- list.sites[[one_n]] %>% anchor_center() %>% mutate(width = 2 * windows[[one_w]])
      PhDfunc::Get1val(one_n,one.w = wig,x = x)
    }) %>% setNames(names(list.sites))
  }) %>% bind_rows(.id = "Condition") %>% mutate(win = one_w)
}) %>% bind_rows()

prof.dat <- mclapply(names(windows),function(one_w){
  lapply(wigs,function(wig){
    lapply(names(list.sites), function(x){
      subdom <- list.sites[[x]]
      d1 <- computeProfile(bed=subdom,wig = wig,seqlens = seqlens,w = windows[[one_w]],span=spans[[one_w]]) %>% colMeans()
      data.frame(Value=d1,Windows=seq(-windows[[one_w]],windows[[one_w]]-spans[[one_w]]+1,spans[[one_w]]),Type=x)
    }) %>%bind_rows()
    
  }) %>%bind_rows(.id = "wig") %>% mutate(win = one_w)
},mc.cores=8) %>% bind_rows()

cutcolors <- c("#FDBECD","black","#BEBEBE")
HRNHEJ = c("#F5AB35","#049372","#BEBEBE")

#cutvsuncutvsrandom

for(mywin in names(windows)){
  pdf(str_c(mywin,"_profiles_drip-seq.pdf"),width=16,height=10)
  cc <- prof.dat %>% 
    filter(win == mywin) %>% 
    filter(!str_detect(wig,"vs")) %>% 
    mutate(Condition = ifelse(str_detect(wig,"24H"),"pOHT_24h","mOHT")) %>% 
    mutate(wig = str_remove(wig,"_24H"))
  ccboxplot <- dat.boxplot.all %>% 
    filter(win == mywin) %>% 
    filter(!str_detect(Condition,"vs")) %>% 
    dplyr::rename(wig = "Condition") %>% 
    mutate(Condition = ifelse(str_detect(wig,"24H"),"pOHT_24h","mOHT")) %>% 
    mutate(wig = str_remove(wig,"_24H"))
  p1 <- cc %>% filter(Type %in% c("cut","uncut","Random80"))%>%ggplot(aes(Windows,Value,colour = Condition)) +
    labs(list(title = "", x = "", y = "")) +
    geom_line()+
    facet_grid(wig~ Type,scales = "free_x") +
    theme_classic() +ggtitle(mywin)
  print(p1)
  
  p2 <- cc %>% filter(Type %in% c("HR","NHEJ","Random30"))%>%ggplot(aes(Windows,Value,colour = Condition)) +
    labs(list(title = "", x = "", y = "")) +
    geom_line()+
    facet_grid(wig~ Type,scales = "free_x") +
    theme_classic() +ggtitle(mywin)
  print(p2)
  
  
  b1 <- ccboxplot %>% filter(Type %in% c("cut","uncut","Random80"))%>%ggplot(aes(Condition,value,fill = Type)) +
    labs(list(title = "", x = "", y = "")) +
    geom_boxplot()+
    facet_grid(wig~ Type) +
    theme_classic() +
    ggtitle(mywin)
  print(b1)
  
  b2 <- ccboxplot %>%filter(Type %in% c("HR","NHEJ","Random30"))%>%ggplot(aes(Condition,value,fill = Type)) +
    labs(list(title = "", x = "", y = "")) +
    geom_boxplot()+
    facet_grid(wig~ Type) +
    theme_classic() +
    ggtitle(mywin)
  print(b2)
  
  dev.off()
}


for(mywin in names(windows)){
  pdf(str_c(mywin,"_bamcompare_profiles_drip-seq.pdf"),width=16,height=10)
  cc <- prof.dat %>% 
    filter(win == mywin) %>% 
    filter(str_detect(wig,"vs")) 
  ccboxplot <- dat.boxplot.all %>% 
    filter(win == mywin) %>% 
    filter(str_detect(Condition,"vs")) 
  p1 <- cc %>% filter(Type %in% c("cut","uncut","Random80"))%>%ggplot(aes(Windows,Value,colour = Type)) +
    labs(list(title = "", x = "", y = "")) +
    geom_line()+
    facet_wrap(~wig,ncol=1,scales = "free_x") +
    scale_color_manual(values=cutcolors) +
    theme_classic() +ggtitle(mywin)
  print(p1)
  
  p2 <- cc %>% filter(Type %in% c("HR","NHEJ","Random30"))%>%ggplot(aes(Windows,Value,colour = Type)) +
    labs(list(title = "", x = "", y = "")) +
    geom_line()+
    facet_wrap(~wig,ncol=1,scales = "free_x") +
    scale_color_manual(values=HRNHEJ) +
    theme_classic() +ggtitle(mywin)
  print(p2)
  
  
  b1 <- ccboxplot %>% filter(Type %in% c("cut","uncut","Random80"))%>%ggplot(aes(Type,value,fill = Type)) +
    labs(list(title = "", x = "", y = "")) +
    geom_boxplot()+
    facet_grid(~ Condition) +
    theme_classic() +
    scale_fill_manual(values=cutcolors) +
    ggtitle(mywin)
  print(b1)
  
  b2 <-ccboxplot %>% filter(Type %in% c("HR","NHEJ","Random30"))%>%ggplot(aes(Type,value,fill = Type)) +
    labs(list(title = "", x = "", y = "")) +
    geom_boxplot()+
    facet_grid(~ Condition) +
    theme_classic() +
    scale_fill_manual(values=HRNHEJ) +
    ggtitle(mywin)
  print(b2)
  
  dev.off()
}

#DIfferent kind of plots bigwig on same tracks


#cutvsuncutvsrandom

for(mywin in names(windows)){
  pdf(str_c(mywin,"_profiles_drip-seq_08_03_2021.pdf"),width=16,height=10)
  cc <- prof.dat %>% 
    filter(win == mywin) %>% 
    filter(!str_detect(wig,"vs")) %>% 
    mutate(Condition = ifelse(str_detect(wig,"24H"),"pOHT_24h","mOHT")) %>% 
    mutate(wig = str_remove(wig,"_24H"))
  ccboxplot <- dat.boxplot.all %>% 
    filter(win == mywin) %>% 
    filter(!str_detect(Condition,"vs")) %>% 
    dplyr::rename(wig = "Condition") %>% 
    mutate(Condition = ifelse(str_detect(wig,"24H"),"pOHT_24h","mOHT")) %>% 
    mutate(wig = str_remove(wig,"_24H"))
  p1 <- cc %>% filter(Type %in% c("cut","uncut","Random80"))%>%ggplot(aes(Windows,Value,colour = wig)) +
    labs(list(title = "", x = "", y = "")) +
    geom_line()+
    facet_grid(Condition~ Type,scales = "free_x") +
    theme_classic() +ggtitle(mywin) + theme(legend.position = "bottom")
  print(p1)
  
  p2 <- cc %>% filter(Type %in% c("HR","NHEJ","Random30"))%>%ggplot(aes(Windows,Value,colour = wig)) +
    labs(list(title = "", x = "", y = "")) +
    geom_line()+
    facet_grid(Condition~ Type,scales = "free_x") +
    theme_classic() +ggtitle(mywin) + theme(legend.position = "bottom")
  print(p2)
  
  
  b1 <- ccboxplot %>% filter(Type %in% c("cut","uncut","Random80"))%>%ggplot(aes(wig,value,fill = wig)) +
    labs(list(title = "", x = "", y = "")) +
    geom_boxplot()+
    facet_grid(Condition~ Type) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(mywin)
  print(b1)
  
  b2 <- ccboxplot %>%filter(Type %in% c("HR","NHEJ","Random30"))%>%ggplot(aes(wig,value,fill = wig)) +
    labs(list(title = "", x = "", y = "")) +
    geom_boxplot()+
    facet_grid(Condition~ Type) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(mywin)
  print(b2)
  
  dev.off()
}


for(mywin in names(windows)){
  pdf(str_c(mywin,"_bamcompare_profiles_drip-seq_08_03_2021.pdf"),width=16,height=10)
  cc <- prof.dat %>% 
    filter(win == mywin) %>% 
    filter(str_detect(wig,"vs")) 
  ccboxplot <- dat.boxplot.all %>% 
    filter(win == mywin) %>% 
    filter(str_detect(Condition,"vs")) 
  p1 <- cc %>% filter(Type %in% c("cut","uncut","Random80"))%>%ggplot(aes(Windows,Value,colour = wig)) +
    labs(list(title = "", x = "", y = "")) +
    geom_line()+
    facet_wrap(~Type,nrow=1,scales = "free_x") +
    theme_classic() +ggtitle(mywin) + theme(legend.position = "bottom")
  print(p1)
  
  p2 <- cc %>% filter(Type %in% c("HR","NHEJ","Random30"))%>%ggplot(aes(Windows,Value,colour = wig)) +
    labs(list(title = "", x = "", y = "")) +
    geom_line()+
    facet_wrap(~Type,nrow=1,scales = "free_x") +
    theme_classic() +ggtitle(mywin) + theme(legend.position = "bottom")
  print(p2)
  
  
  b1 <- ccboxplot %>% filter(Type %in% c("cut","uncut","Random80"))%>%ggplot(aes(Condition,value,fill = Condition)) +
    labs(list(title = "", x = "", y = "")) +
    geom_boxplot()+
    facet_grid(~ Type) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(mywin)
  print(b1)
  
  b2 <-ccboxplot %>% filter(Type %in% c("HR","NHEJ","Random30"))%>%ggplot(aes(Condition,value,fill = Condition)) +
    labs(list(title = "", x = "", y = "")) +
    geom_boxplot()+
    facet_grid(~ Type) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(mywin)
  print(b2)
  
  dev.off()
}







#PROFILS GENES

require(tidyverse)
require(plyranges)
require(rtracklayer)
mydir <- "MATRIX/"
files <- list.files(mydir)
dataplot <- lapply(files,function(onefile){
  read_tsv(str_c(mydir,onefile),skip = 1,col_names = F) %>% dplyr::select(-X1:-X6) %>% colMeans()  %>%
    enframe() %>% mutate(seq = 1:320)%>%
    mutate(wig = str_remove(onefile,".tab.gz")) %>% separate(wig,into = c("Bed","File"),sep="\\.")
}) %>% bind_rows()


dataplot.pm <- dataplot %>% filter(!str_detect(Bed,"vs")) %>% 
  mutate(Condition = ifelse(str_detect(Bed,"24H"),"pOHT_24h","mOHT")) %>% 
  mutate(Bed = str_remove(Bed,"_24H"))

p1 <- dataplot.pm %>% ggplot(aes(seq,value,col=Condition,fill=Condition)) +
  geom_line() +
  geom_ribbon(aes(ymin=0,ymax=value,x=seq),alpha=0.4) + 
  geom_vline(xintercept = 60,linetype="longdash") +
  geom_vline(xintercept = 260,linetype="longdash") +
  scale_x_continuous(name = 'Position',
                     breaks = c(1,60,260,320),
                     labels = c("TSS-3kb", 'TSS', 'TES', 'TES+3kb')
  ) +
  facet_grid(Bed~File) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom") + ggtitle("MOHT vs pOHT for Genes Cat")

dataplot.pm <- dataplot %>% filter(str_detect(Bed,"vs"))

p2 <- dataplot.pm %>% ggplot(aes(seq,value,col=File,fill=File)) +
  geom_line() +
  geom_ribbon(aes(ymin=0,ymax=value,x=seq),alpha=0.4) + 
  geom_vline(xintercept = 60,linetype="longdash") +
  geom_vline(xintercept = 260,linetype="longdash") +
  scale_x_continuous(name = 'Position',
                     breaks = c(1,60,260,320),
                     labels = c("TSS-3kb", 'TSS', 'TES', 'TES+3kb')
  ) +
  facet_grid(Bed~File) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom") + ggtitle("Bamcompare for Genes Cat")



pdf("Profile_drip-seq_Genes.pdf",width=12,height=18)
print(p1)
print(p2)
dev.off()

#UPDATE AU 09/03/21

require(dplyr)
require(rtracklayer)
require(ggplot2)
require(plyranges)
library("BSgenome.Hsapiens.UCSC.hg19")
seqlens = seqlengths( Hsapiens );
#Load AsiSI locations
asi = import.bed("../data/ASIsites_hg19.bed")
bless80 = import.bed("../data/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
HR = import.bed("../data/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed")


list.sites <- list("cut"=bless80,"HR"=HR)

windows = list("50kb"=50000)
spans = list("50kb"=1000)


wig.List <- list(
  "HKKWHBGX7_DRIP_DIvA"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HKKWHBGX7/PROCESSED/mapping/BIGWIG/HKKWHBGX7_DRIP_C1_18s002883-1-1_Clouaire_lane118s002883_sequence_normalized.bw",
  "HKKWHBGX7_DRIP_pOHT24H"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HKKWHBGX7/PROCESSED/mapping/BIGWIG/HKKWHBGX7_DRIP_C1_24H_18s002884-1-1_Clouaire_lane118s002884_sequence_normalized.bw",
  "HYG77BGX2_DRIPseq_DIvA"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/HYG77BGX2_DRIPseq_SETX/PROCESSED/WIGGLE/hg19/HYG77BGX2_DRIP2_C1_DIVA_17s002460-1-1_Clouaire_lane117s002460_normalized_hg19_nodups.bw",
  "HYG77BGX2_DRIPseq_pOHT4H"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/HYG77BGX2_DRIPseq_SETX/PROCESSED/WIGGLE/hg19/HYG77BGX2_DRIP2_C1_OHT_17s002461-1-1_Clouaire_lane117s002461_normalized_hg19_nodups.bw",
  "RAD51_pOHT4H"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/XRCC4_RAD51/EMBL3_HJKVVBGX5_RAD51_XRCC4_1h_4h_24h/PROCESSED/mapping/BIGWIG/HJKVVBGX5_RAD51_4H_17s006148-1-1_Clouaire_lane117s006148_sequence_normalized.bw",
  "RAD51_pOHT24H"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/XRCC4_RAD51/EMBL3_HJKVVBGX5_RAD51_XRCC4_1h_4h_24h/PROCESSED/mapping/BIGWIG/HJKVVBGX5_RAD51_24H_17s006149-1-1_Clouaire_lane117s006149_sequence_normalized.bw"
)

wigs <- lapply(wig.List,import.bw,as="RleList") %>% setNames(names(wig.List))



prof.dat <- mclapply(names(windows),function(one_w){
  lapply(wigs,function(wig){
    lapply(names(list.sites), function(x){
      subdom <- list.sites[[x]]
      d1 <- computeProfile(bed=subdom,wig = wig,seqlens = seqlens,w = windows[[one_w]],span=spans[[one_w]]) %>% colMeans()
      data.frame(Value=d1,Windows=seq(-windows[[one_w]],windows[[one_w]]-spans[[one_w]]+1,spans[[one_w]]),Type=x)
    }) %>%bind_rows()
    
  }) %>%bind_rows(.id = "wig") %>% mutate(win = one_w)
},mc.cores=8) %>% bind_rows()

cc <- prof.dat %>% 
  mutate(Condition = str_extract(wig,"DIvA|pOHT24H|pOHT4H")) %>% 
  mutate(wig = str_remove(wig,"_DIvA|_pOHT24H|_pOHT4H"))

ccbless <- cc %>% filter(Type == "cut",wig != "RAD51") %>% ggplot(aes(x=Windows,y=Value,col=Condition)) + geom_line() + facet_wrap(~wig,scales="free")
RAD51.bless <- cc %>% filter(Type == "cut",wig == "RAD51") %>% ggplot(aes(x=Windows,y=Value,col=Condition)) + geom_line() + facet_wrap(~wig,scales="free")

pbless <- cowplot::plot_grid(RAD51.bless,ccbless,nrow=2)

cchr <- cc %>% filter(Type == "HR",wig != "RAD51") %>% ggplot(aes(x=Windows,y=Value,col=Condition)) + geom_line() + facet_wrap(~wig,scales="free")
RAD51.hr <- cc %>% filter(Type == "HR",wig == "RAD51") %>% ggplot(aes(x=Windows,y=Value,col=Condition)) + geom_line() + facet_wrap(~wig,scales="free")

phr <- cowplot::plot_grid(RAD51.hr,cchr,nrow=2)

pdf("Profile_drip-seq_rad51_4h24_bless.pdf",width=12,height=6)
print(pbless)
dev.off()
pdf("Profile_drip-seq_rad51_4h24_hr.pdf",width=12,height=6)
print(phr)
dev.off()

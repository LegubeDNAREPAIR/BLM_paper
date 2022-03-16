computeProfileAugmented <- function (bed, wig, w = 20000, span = 200, seqlens, row.method = rowMeans, 
          col.method = colMeans, heatmap = F) 
{
  if (class(wig) != "SimpleRleList") {
    stop("ERROR : unknown class of wig, please provide a SimpleRleList")
  }
  s.asi <- sortSeqlevels(bed)
  names(s.asi) <- bed$name
  res <- do.call("rbind", lapply(split(s.asi, seqnames(s.asi)), 
                                 function(x) {
                                   cov <- wig[[unique(as.character(seqnames(x)))]]
                                   center <- x %>% as_tibble() %>% mutate(center = ifelse(width == 
                                                                                            1, start, start(x) + (width(x)/2))) %>% pull(center)
                                   Views(cov, start = center - w, end = center + w) %>% 
                                     as.matrix()
                                 }))
  if (span != 1) {
    sts = seq(1, ncol(res) - span + 1, span)
    eds = seq(span, ncol(res), span)
    sub.res <- do.call("cbind", lapply(1:length(sts), function(j) {
      row.method(res[, sts[j]:eds[j]])
    }))
  }
  else {
    sub.res <- res
  }
  if (heatmap == T) {
    colnames(sub.res) <- seq(1, ncol(res) - span + 1, span)
    rownames(sub.res) <- bed$name
    sub.res <- sub.res %>% reshape2::melt() %>% setNames(c("Name", 
                                                           "Windows", "value"))
    return(sub.res)
  }
  return(tibble(Window = seq(-w, w - span + 1, span), Value = col.method(sub.res)))
}
require(tidyverse)
require(cowplot)
require(plyranges)
require(rtracklayer)
require(BSgenome.Hsapiens.UCSC.hg19.masked)
seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
mes_bw <- c("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/BLM/BLM_4H_BGI_RUN1_201203/4H/PROCESSED/ALIGNED/WIGGLE/blm_normalized_hg19.bw",
            "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/GAMMA/U2OS/BGI_RUN1_201203/PROCESSED/ALIGNED/WIGGLE/hg19/GAM.clean_normalized_01022018.bw"
)%>% setNames(c("BLM","GAMMA")) %>% map(import.bw,as="RleList")

mes_beds <- list(
  bless80 = read_bed("../dataBLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
)

res.prof <- lapply(mes_beds,function(one_bed){
  lapply(names(mes_bw),function(one.w){
    message(one.w)
    computeProfileAugmented(wig=mes_bw[[one.w]],bed=one_bed,w=50000,span=200,seqlens = seqlens,heatmap = T) %>% mutate(wig = one.w) 
  }) %>% bind_rows()
}) %>% bind_rows(.id = "Type")
#Sort by GAMMA
gamma.values <- res.prof %>%
  filter(wig =="GAMMA") %>%
  group_by(Name) %>% summarise(sumVal = sum(value)) %>% arrange((sumVal))
background = "black"
highvalues = "#c23616"

plot_func <- function(my_data,label){
  min.treshold <- 0
  max.treshold <- quantile(my_data$value,0.99)
  message(max.treshold)
  my_data %>%
    mutate(value = ifelse(value > max.treshold,max.treshold,value)) %>% 
    ggplot(aes(x=Windows,y=Name,fill=value)) + geom_tile() +theme_classic(base_size=18) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),legend.position="bottom") +ggtitle(label) +
    scale_fill_gradient(low = background,high = highvalues,na.value=background,limits=c(min.treshold, max.treshold))
}


f.res.prof <- res.prof %>% mutate(Name = factor(Name,levels=gamma.values$Name)) %>% group_by(wig) %>% nest() %>% 
  mutate(plots = map2(data,wig, plot_func))

p<- cowplot::plot_grid(plotlist = f.res.prof$plots,ncol=2)
pdf("PLOTS/FIG1D_sorted_by_gamma_bless80.pdf",width=12,height = 18)
print(p)
dev.off()
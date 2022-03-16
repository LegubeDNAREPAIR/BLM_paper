##Fig 4D 4E S1G S6B

ccprof <- function(bed, wig, w = 20000, span = 200, seqlens, method = "mean") 
{
  if (class(wig) != "SimpleRleList") {
    stop("ERROR : unknown class of wig, please provide a SimpleRleList")
  }
  mat = NULL
  for (i in 1:length(bed)) {
    message(i, "/", length(bed))
    bedi = bed[i, ]
    chr = as.character(seqnames(bedi))
    cov = wig[[chr]]
    center = start(bedi) + 4
    stW = center - w
    edW = center + w
    if (span == 1) {
      vm = as.numeric(Views(cov, start = stW, end = edW)[[1]])
    }
    else {
      sts = seq(stW, edW - span + 1, span)
      eds = seq(stW + span - 1, edW, span)
      v = Views(cov, start = sts, end = eds)
      if (method == "sum") {
        vm = sum(v)
      }
      else {
        vm = mean(v)
      }
      vm[sts >= seqlens[chr] | sts > length(cov)] = 0
      vm[sts < 1] = 0
    }
    mat = rbind(mat, vm)
  }
  return(mat)
}

bless80 <- read_bed("../data/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
HR <-  read_bed("../data/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed")
NHEJ <-  read_bed("../data/BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed")

FuncManage2 <- function(x,my.w = 2000,my.span = 50){
  dir_x <- dirname(x)
  my.bw <- glue::glue("{x}") %>% import.bw(as="RleList")
  x <- basename(x)
  res <- ccprof(bed=bless80,wig = my.bw,w=my.w,span=my.span,seqlens=seqlens)
  colnames(res) <- seq(-my.w,my.w-my.span+1,my.span)
  rownames(res) <- glue::glue("{bless80$name}")
  
  res <- res %>% reshape2::melt() %>% dplyr::rename(Windows = "Var2",Name = "Var1") %>% mutate(name = x)
  return(res)
}

bw.files <- "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/EdU-Seq/Clouaire_HVHVTBGXG" %>% 
  list.files(recursive = T,pattern=".bw",full.names = T) 

my.w <- 5000
my.span <- 50
res_f <- mclapply(bw.files,FuncManage2,my.w=my.w,my.span=my.span,mc.cores=8) %>% bind_rows()

res_f <- res_f %>% mutate(Cat = case_when(
  Name %in% HR$name ~ "HR",
  Name %in% NHEJ$name ~ "NHEJ",
  TRUE ~ "Other"
))



res2 <- res_f %>% filter(!str_detect(name,"bamCompare"))
res2 <- res2 %>% group_by(name,Windows) %>% summarise(value = mean(value)) %>% 
  mutate(name = str_remove(name,".+Aude"))
p2 <- res2 %>% ggplot(aes(x=Windows,y=value,col=name)) + geom_line() + theme(legend.position="bottom") + scale_color_brewer(palette="Set2") +
  geom_vline(xintercept = 0,linetype="dashed")

pdf(glue::glue("Edu_for_gaelle_10_11_2021_{my.w}_{my.span}_bless80.pdf"),width = 12,height=8)
print(p2)
dev.off()


bw.files.1 <- "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/EdU-Seq/Clouaire_HVHVTBGXG" %>% 
  list.files(recursive = T,pattern=".bw",full.names = T) 
bw.files.2 <- c(
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/BLM/BLM_4H_BGI_RUN1_201203/4H/PROCESSED/ALIGNED/WIGGLE/BLM_normalized_01_02_2018.bw",
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Tot-Pol2/PROCESSED/BY_US/WIGGLE/H3LNWBGX3_Tot_Pol2_DIVA_17s001632-1-1_Clouaire_lane117s001632_sequence_normalized.bw",
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Ser2P_Pol2/PROCESSED/BY_US/WIGGLE/H3LNWBGX3_Ser2P_Pol2_DIVA_17s001630-1-1_Clouaire_lane117s001630_sequence_normalized.bw",
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/HYG77BGX2_DRIPseq_SETX/PROCESSED/WIGGLE/hg19/HYG77BGX2_DRIP2_C1_OHT_17s002461-1-1_Clouaire_lane117s002461_hg19_nodups.bw"
)

Get1val <- function (Name, one.w, x) 
{
  require(magrittr)
  lapply(split(x, droplevels(seqnames(x))), function(zz) {
    message(unique(as.character(seqnames(zz))))
    cov <- one.w[[unique(as.character(seqnames(zz)))]]
    score <- IRanges::Views(cov, start = start(zz), end = end(zz)) %>% 
      mean()
    tibble::tibble(wig = Name, value = score, rowname = zz$name)
  }) %>% bind_rows()
}
FuncManage3 <- function(x,my.w = 2000){
  dir_x <- dirname(x)
  my.bw <- glue::glue("{x}") %>% import.bw(as="RleList")
  mon_bed <- bless80 %>% anchor_center() %>% mutate(width=my.w*2) %>% as_granges()
  res <- Get1val(Name=basename(x),one.w = my.bw,x=mon_bed)
  return(res)
}

my.w <- 5000
res_f_1 <- mclapply(bw.files.1,FuncManage3,my.w=my.w,mc.cores=8) %>% bind_rows()
my.w <- 1000
res_f_2 <- mclapply(bw.files.2,FuncManage3,my.w=my.w,mc.cores=8) %>% bind_rows()
res_f <- rbind(res_f_1,res_f_2) %>%
  filter(!str_detect(wig,"bamCompare")) %>% 
  mutate(wig = str_remove(wig,"HVHVTBGXG_EdU-IP_47_20s003953-1-1_Clouaire_lane1")) %>% 
  mutate(wig = str_remove(wig,"_hg19_nodups.bw|_sequence_normalized.bw|.bw")) %>% 
  spread(key = wig,value=value)

p1_scat <- res_f %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_matrix(cols=vars(`BLM_normalized_01_02_2018`),rows=vars(`H3LNWBGX3_Tot_Pol2_DIVA_17s001632-1-1_Clouaire_lane117s001632`))


p3_scat <- res_f %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_matrix(cols=vars(`BLM_normalized_01_02_2018`),rows=vars(`H3LNWBGX3_Ser2P_Pol2_DIVA_17s001630-1-1_Clouaire_lane117s001630`))



pdf(glue::glue("BLM_scatter_for_gaelle_10_11_2021_bless80.pdf"),width = 6,height=6)
print(p1_scat)
print(p3_scat)
dev.off()






bw.files.3 <- c("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/BGI_POLII-S2P/PROCESSED/ALIGNED/WIGGLE/polII_normalized_hg19.bw",
                "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Tot-Pol2/PROCESSED/BY_US/WIGGLE/H3LNWBGX3_Tot_Pol2_DIVA_17s001632-1-1_Clouaire_lane117s001632_sequence_normalized.bw",
                "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Ser2P_Pol2/PROCESSED/BY_US/WIGGLE/H3LNWBGX3_Ser2P_Pol2_DIVA_17s001630-1-1_Clouaire_lane117s001630_sequence_normalized.bw",
                "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Ser2P_Pol2/PROCESSED/BY_US/WIGGLE/H3LNWBGX3_Ser2P_Pol2_OHT_17s001631-1-1_Clouaire_lane117s001631_sequence_normalized.bw",
                "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/HYG77BGX2_DRIPseq_SETX/PROCESSED/WIGGLE/hg19/HYG77BGX2_DRIP2_C1_OHT_17s002461-1-1_Clouaire_lane117s002461_hg19_nodups.bw",
                "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/bigwigCompare/HYG77BGX2_DRIPseq_SETX/HYG77BGX2_CTRL_pOHT_diff_CTRL_mOHT_bigwigCompareDiff.bw"
)

my.w <- 2000
res_f_3 <- mclapply(bw.files.3,FuncManage3,my.w=my.w,mc.cores=8) %>% bind_rows()

#Edu
#pm 1KB
pEdu1_scat <- res_f %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_matrix(cols=vars(`EdU47AudemimsiC20hOHT`),rows=vars(`H3LNWBGX3_Tot_Pol2_DIVA_17s001632-1-1_Clouaire_lane117s001632`))


pEdu2_scat <- res_f %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_matrix(cols=vars(`EdU47AudemimsiC20hOHT`),rows=vars(`H3LNWBGX3_Ser2P_Pol2_DIVA_17s001630-1-1_Clouaire_lane117s001630`))


pdf(glue::glue("EdU_pm1kb_scatter_for_gaelle_10_11_2021_bless80.pdf"),width = 6,height=6)

print(pEdu1_scat)
print(pEdu2_scat)
dev.off()

#pm2KB

res_f2kb <- rbind(res_f_1,res_f_3) %>%
  filter(!str_detect(wig,"bamCompare")) %>% 
  mutate(wig = str_remove(wig,"HVHVTBGXG_EdU-IP_47_20s003953-1-1_Clouaire_lane1")) %>% 
  mutate(wig = str_remove(wig,"_hg19_nodups.bw|_sequence_normalized.bw|.bw")) %>% 
  spread(key = wig,value=value)

pEdu3_scat <- res_f2kb %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_matrix(cols=vars(`EdU47AudemimsiC20hOHT`),rows=vars(`H3LNWBGX3_Tot_Pol2_DIVA_17s001632-1-1_Clouaire_lane117s001632`))


pEdu4_scat <- res_f2kb %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_matrix(cols=vars(`EdU47AudemimsiC20hOHT`),rows=vars(`H3LNWBGX3_Ser2P_Pol2_DIVA_17s001630-1-1_Clouaire_lane117s001630`))

pdf(glue::glue("EdU_pm2kb_scatter_for_gaelle_10_11_2021_bless80.pdf"),width = 6,height=6)

print(pEdu3_scat)
print(pEdu4_scat)
dev.off()

#DRIP
#pm 1kb
pDRIP_scat <- res_f %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_matrix(cols = vars(`HYG77BGX2_DRIP2_C1_OHT_17s002461-1-1_Clouaire_lane117s002461`),rows=vars(`EdU47AudemimsiC20hOHT`))

pDRIP2_scat <- res_f %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_matrix(cols = vars(`HYG77BGX2_CTRL_pOHT_diff_CTRL_mOHT_bigwigCompareDiff`),rows=vars(`EdU47AudemimsiC20hOHT`))



pdf(glue::glue("DRIP_pm1kb_scatter_for_gaelle_10_11_2021_bless80.pdf"),width = 6,height=6)
print(pDRIP_scat)
print(pDRIP2_scat)
dev.off()

#pm 2kb
pDRIP3_scat <- res_f2kb %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_matrix(cols = vars(`HYG77BGX2_DRIP2_C1_OHT_17s002461-1-1_Clouaire_lane117s002461`),rows=vars(`EdU47AudemimsiC20hOHT`))

pDRIP4_scat <- res_f2kb %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_matrix(cols = vars(`HYG77BGX2_CTRL_pOHT_diff_CTRL_mOHT_bigwigCompareDiff`),rows=vars(`EdU47AudemimsiC20hOHT`))



pdf(glue::glue("DRIP_pm2kb_scatter_for_gaelle_10_11_2021_bless80.pdf"),width = 6,height=6)
print(pDRIP3_scat)
print(pDRIP4_scat)
dev.off()

cc.corr <- res_f %>% dplyr::select(-rowname) %>% correlate()

cc.corr.blm <- cc.corr %>% focus(`BLM_normalized_01_02_2018`) %>% filter(str_detect(rowname,"H3LNWBGX3|polII"))

cc.corr.DRIP <- cc.corr %>% focus(`EdU47AudemimsiC20hOHT`) %>% filter(!str_detect(rowname,"EdU47|BLM|polII"))


cc.corr.2kb <- res_f2kb %>% dplyr::select(-rowname) %>% correlate()
cc.corr.DRIP2kb <- cc.corr.2kb %>% focus(`EdU47AudemimsiC20hOHT`) %>% filter(!str_detect(rowname,"EdU47|BLM|polII"))



# scatterplot Edu (-5/+5) vs DRIP OHT -5/+5
# scatterplot Edu SETXmim20h(-5/+5) vs DRIPOHT SETX -5/+5
bw.files <- c(
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/EdU-Seq/Clouaire_HVHVTBGXG/PROCESSED/mapping/BIGWIG/HVHVTBGXG_EdU-IP_47_20s003953-1-1_Clouaire_lane1EdU47AudemimsiSETX20hOHT_sequence_normalized.bw",
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/EdU-Seq/Clouaire_HVHVTBGXG/PROCESSED/mapping/BIGWIG/HVHVTBGXG_EdU-IP_47_20s003953-1-1_Clouaire_lane1EdU47AudemimsiC20hOHT_sequence_normalized.bw",
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/HYG77BGX2_DRIPseq_SETX/PROCESSED/WIGGLE/hg19/HYG77BGX2_DRIP2_C1_OHT_17s002461-1-1_Clouaire_lane117s002461_normalized_hg19_nodups.bw",
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/HYG77BGX2_DRIPseq_SETX/PROCESSED/WIGGLE/hg19/HYG77BGX2_DRIP2_STX_OHT_17s002463-1-1_Clouaire_lane117s002463_normalized_hg19_nodups.bw"
)

Get1val <- function (Name, one.w, x) 
{
  require(magrittr)
  lapply(split(x, droplevels(seqnames(x))), function(zz) {
    message(unique(as.character(seqnames(zz))))
    cov <- one.w[[unique(as.character(seqnames(zz)))]]
    score <- IRanges::Views(cov, start = start(zz), end = end(zz)) %>% 
      mean()
    tibble::tibble(wig = Name, value = score, rowname = zz$name)
  }) %>% bind_rows()
}
FuncManage3 <- function(x,my.w = 2000){
  dir_x <- dirname(x)
  my.bw <- glue::glue("{x}") %>% import.bw(as="RleList")
  mon_bed <- bless80 %>% anchor_center() %>% mutate(width=my.w*2) %>% as_granges()
  res <- Get1val(Name=basename(x),one.w = my.bw,x=mon_bed)
  return(res)
}

my.w <- 5000
res_f <- mclapply(bw.files,FuncManage3,my.w=my.w,mc.cores=8) %>% bind_rows()
res_f <- res_f %>%
  mutate(wig = str_remove(wig,"HVHVTBGXG_EdU-IP_47_20s003953-1-1_Clouaire_lane1")) %>% 
  mutate(wig = str_remove(wig,"_hg19_nodups.bw|_17s[0-9]+.+_normalized_hg19_nodups.bw|.bw")) %>% 
  spread(key = wig,value=value)

require(ggforce)
pDRIP_scat <- res_f %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_matrix(cols = vars(`HYG77BGX2_DRIP2_C1_OHT`),rows=vars(`EdU47AudemimsiC20hOHT_sequence_normalized`))

pDRIP2_scat <- res_f %>% 
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_matrix(cols = vars(`HYG77BGX2_DRIP2_STX_OHT`),rows=vars(`EdU47AudemimsiSETX20hOHT_sequence_normalized`))



pdf(glue::glue("DRIP_pm1kb_scatter_for_gaelle_01_12_2021_bless80.pdf"),width = 6,height=6)
print(pDRIP_scat)
print(pDRIP2_scat)
dev.off()

require(corrr)
cc.corr <- res_f %>% dplyr::select(-rowname) %>% correlate()

cc.corr %>% focus(`HYG77BGX2_DRIP2_STX_OHT`) %>% filter(rowname == "EdU47AudemimsiSETX20hOHT_sequence_normalized")
cc.corr %>% focus(`HYG77BGX2_DRIP2_C1_OHT`) %>% filter(rowname == "EdU47AudemimsiC20hOHT_sequence_normalized")

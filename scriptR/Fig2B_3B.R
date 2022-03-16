require(reshape2)
require(tidyverse)
require(rtracklayer)
matrices.5k <- list.files("/mnt/NAS2/DOSSIERS_PERSO/VINCENT/PROJET_RAPH_13_04_18/MATRIX",pattern="_5000.gz",full.names = T)
w <- 5000
bin <- 10
col <- seq(-w,w-bin+1,bin)

#Get data
HR.sites <- import.bed("/mnt/NAS/DATA/AsiSI/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed") 
NHEJ.sites <- import.bed("/mnt/NAS/DATA/AsiSI/BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed") 
res.bless <- read.csv("../data/BLESS80_ordering_by_ratioBLM2kb_BLESS500bp.csv",h=T) %>% arrange(desc(ratio.BLM.BLESS))
res.bless$Top_ratio= "NA"
res.bless[1:20,]$Top_ratio = "BLM-high"
res.bless[61:80,]$Top_ratio = "BLM-low"

#Compute prof

for(vw in c(500,2000,5000)){
  interval.w <- seq(-vw,vw-bin,bin)
  col.of.interest <- col %in% interval.w
  tmp.SETX = read.table(gzfile(matrices.5k[grep("SETX",matrices.5k)][2]),skip=1,h=F) %>% column_to_rownames(var="V4") %>% select(-V1,-V2,-V3,-V5,-V6)  %>% .[,col.of.interest] 
  tmp.BLM = read.table(gzfile(matrices.5k[grep("BLM",matrices.5k)]),skip=1,h=F) %>% column_to_rownames(var="V4") %>% select(-V1,-V2,-V3,-V5,-V6)  %>% .[,col.of.interest] 
  tmp.RAD51 = read.table(gzfile(matrices.5k[grep("RAD51",matrices.5k)]),skip=1,h=F) %>% column_to_rownames(var="V4") %>% select(-V1,-V2,-V3,-V5,-V6)  %>% .[,col.of.interest] 
  
  colnames(tmp.SETX) = colnames(tmp.BLM) = colnames(tmp.RAD51) = as.character(interval.w)
  
  #Correlation
  dat.plot <- cbind(SETX=rowSums(tmp.SETX),BLM=rowSums(tmp.BLM),RAD51=rowSums(tmp.RAD51)) %>% 
    as.data.frame() %>% tibble::rownames_to_column(var="SITES") %>% 
    filter(SITES %in% res.bless$SITES) %>%
    mutate(Type = ifelse(SITES %in% HR.sites$name,"HR",
                         ifelse(SITES %in% NHEJ.sites$name,"NHEJ","NA")))
  cor.res.p <- paste("Pearson correlation",cor(dat.plot$SETX,dat.plot$BLM),sep=": ")
  cor.res.s <- paste("Spearman correlation",cor(dat.plot$SETX,dat.plot$BLM,method = "spearman"),sep=": ")
  p <- ggplot(dat.plot,aes(x=SETX,y=BLM,col=Type)) + geom_point() + theme_classic() + scale_color_manual(values = c("red","black","blue")) + ggtitle(cor.res.p)
  pdf(paste("ScatterplotSETXoverBLM",vw,"80bless.pdf"),height=12,width=12)
  print(p)
  dev.off()
  cor.res.p <- paste("Pearson correlation",cor(dat.plot$SETX,dat.plot$RAD51),sep=": ")
  cor.res.s <- paste("Spearman correlation",cor(dat.plot$SETX,dat.plot$RAD51,method = "spearman"),sep=": ")
  
  #Profile
  ##BLM/HIGH/LOW
  BLM.high.BLM <- tmp.BLM %>% .[rownames(tmp.BLM) %in% res.bless[which(res.bless$Top_ratio == "BLM-high"),]$SITES,] %>% colSums()
  BLM.low.BLM <- tmp.BLM %>% .[rownames(tmp.BLM) %in% res.bless[which(res.bless$Top_ratio == "BLM-low"),]$SITES,] %>% colSums()
  BLM.high.RAD51 <- tmp.RAD51 %>% .[rownames(tmp.RAD51) %in% res.bless[which(res.bless$Top_ratio == "BLM-high"),]$SITES,] %>% colSums()
  BLM.low.RAD51 <- tmp.RAD51 %>% .[rownames(tmp.RAD51) %in% res.bless[which(res.bless$Top_ratio == "BLM-low"),]$SITES,] %>% colSums()
  
  
  
  dat.plot <- rbind( data.frame("value"=BLM.high.BLM,"Windows"=interval.w,"Type"="BLM-high","File"="BLM"),
                     data.frame("value"=BLM.low.BLM,"Windows"=interval.w,"Type"="BLM-low","File"="BLM"),
                     data.frame("value"=BLM.high.RAD51,"Windows"=interval.w,"Type"="BLM-high","File"="RAD51"),
                     data.frame("value"=BLM.low.RAD51,"Windows"=interval.w,"Type"="BLM-low","File"="RAD51")
  )
  
  pdf(paste("Prof_RAD51_BLM",vw,"80bless_BLM_HIGHLOW.pdf",sep="_"),height=12,width=12)
  p <-ggplot(dat.plot,aes(Windows,value,colour=File)) + geom_line() + theme_classic()+ facet_wrap(~Type,nrow=2)
  print(p)
  dev.off()
  
  
  
}



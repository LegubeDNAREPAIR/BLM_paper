#Profile TOP High/LOW BLM for BLM
BLM.high.sites <- res.bless[which(res.bless$Top_ratio == "BLM-high"),]$SITES
BLM.high.profile <- computeProfile(bed = bless80[which(bless80$name %in% BLM.high.sites)],wig = myWigBLM$pOHT,w = 10000,span=200,seqlens = seqlens)
BLM.low.sites <- res.bless[which(res.bless$Top_ratio == "BLM-low"),]$SITES
BLM.low.profile <- computeProfile(bed = bless80[which(bless80$name %in% BLM.low.sites)],wig = myWigBLM$pOHT,w = 10000,span=200,seqlens = seqlens)

w = 10000
span = 200
data = rbind(data.frame("Windows"=seq(-w,w-span+1,span),"value"=colMeans(BLM.high.profile),"Type"="BLM-High",Histone = "BLM"),
             data.frame("Windows"=seq(-w,w-span+1,span),"value"=colMeans(BLM.low.profile),"Type"="BLM-Low",Histone = "BLM"))


p <- ggplot(na.omit(data),aes(Windows,value,col=Type)) +
  labs(list(title = "", x = "", y = "")) +
  geom_line()+
  #facet_wrap(~ Type,nrow = 2,scales = "free_x") +
  theme_classic()
pdf("FIG_prof_BLM_BLM-High_BLM-Low_10kb_pOHT.pdf",height=12,width=16)
print(p + facet_wrap(~ Type,nrow = 2,scales = "free_x"))
dev.off()
pdf("FIG_prof_superpose_BLM_BLM-High_BLM-Low_10kb_pOHT.pdf",height=12,width=16)
print(p)
dev.off()

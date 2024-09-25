
##fig3a
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
gexpr <- subset(gexpr,lineage=="mgeLin")
gexpr$group <- gexpr$lobe
gexpr$group[which(gexpr$region == "MGE")] <- "MGE"
gexpr$group[which(gexpr$region == "CGE")] <- "CGE"
gexpr$group[which(gexpr$region == "LGE")] <- "LGE"
gexpr <- subset(gexpr,group %in% c("GE","NCX"),invert=T)
dim(gexpr)
unique(gexpr$group)

mydata <- FetchData(gexpr,vars = c("group","species"))
mydata <- mydata %>% group_by(species,group) %>% summarise(counts=n())
mydata <- mydata %>% group_by(species) %>% mutate(sums=sum(counts))
mydata$Freq <- round(mydata$counts/mydata$sums,4)

groupRef <- c("MGE","CGE","LGE","FC","MSC","TC","OC","Insula","HIP","AMY","STR")
mydata$group <- factor(mydata$group,levels = groupRef)
mydata <- mydata[order(mydata$group),]


cols=c("#f768a1","#c51b8a","#7a0177","#fdd49e","#fdbb84","#fc8d59","#e34a33","#b30000","#f6e8c3","#d8b365","#8c510a")
names(cols)=groupRef

mydata1=mydata[mydata$species=="Human",]
mydata2=mydata[mydata$species=="Monkey",]


p1=ggplot(mydata1,aes(x="",y=Freq,fill=group))+geom_bar(stat = "identity")+labs(title = "Human",x = '', y = '')+
  scale_fill_manual(values =cols[mydata1$group])+coord_polar(theta = "y")+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"),panel.background = element_blank(),legend.position = "none"
  )+
  geom_text(aes(y=cumsum(rev(Freq))-rev(Freq)/2,x=1.5,label=rev(percent(mydata1$Freq))),size=2)

p2=ggplot(mydata2,aes(x="",y=Freq,fill=group))+geom_bar(stat = "identity")+labs(title = "Monkey",x = '', y = '')+
  scale_fill_manual(values =cols[mydata2$group])+coord_polar(theta = "y")+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,face = "bold"),panel.background = element_blank())+
  geom_text(aes(y=cumsum(rev(Freq))-rev(Freq)/2,x=1.5,label=rev(percent(mydata2$Freq))),size=2)



pdf("fig3A_Lobe_pie.pdf",6,6)
par(mar=c(0.1,0.1,0.1,0.1))
p1+p2
dev.off()

####################################################################
##fig3b
library(Seurat)
library(dplyr)
library(ggplot2)
library(corrplot)
library(RColorBrewer)
library(openxlsx)
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
gexpr <- subset(gexpr,lineage=="mgeLin")
gexpr$group <- gexpr$lobe
gexpr$group[which(gexpr$region == "MGE")] <- "MGE"
gexpr$group[which(gexpr$region == "CGE")] <- "CGE"
gexpr$group[which(gexpr$region == "LGE")] <- "LGE"
gexpr <- subset(gexpr,group %in% c("GE","NCX"),invert=T)
gexpr <- subset(gexpr,labsite != "Kriegstein2022")
groupRef <- c("MGE","CGE","LGE","FC","MSC","TC","OC","Insula","HIP","AMY","STR")
unique(gexpr$group)
gexpr$group <- factor(gexpr$group,levels = groupRef)
dim(gexpr)


pdf("fig3B_cor.pdf",10,5)
par(omi=c(0.1,0.1,0.1,0.1))
par(mfrow=c(1,2))

set.seed(1234)
gexpr5 = subset (gexpr, species=="Human")
Idents(gexpr5) <- "group"
recMarkers <- FindAllMarkers(gexpr5,only.pos= T, min.pct = 0.1, logfc.threshold = log(1.25),verbose = F,max.cells.per.ident = 2000)
write.xlsx(recMarkers,file = "2.fig3B_Human_cor.xlsx")
recMarkers <- subset(recMarkers[grep("^RP[L|S]",recMarkers$gene, ignore.case = FALSE,invert=TRUE),],subset=p_val_adj < 0.05)
region_av <- AverageExpression(gexpr5,group.by = "group",features = unique(recMarkers$gene),assays = "RNA") 
region_av <- region_av$RNA
gene_cell_exp <- t(scale(t(region_av),scale = T,center = T))

dat1_cor <- cor(as.matrix(gene_cell_exp))
p1=corrplot(dat1_cor,order = "original", type = "full", addrect = 4,tl.pos = "lt", tl.cex=.8,col = rev(brewer.pal(n=8, name="RdYlBu")),
            tl.col = "black",cl.cex = 0.7,cl.pos = "r",cl.length = 5,title="Human",
            cl.ratio=0.1,mar=c(3, 0, 2, 0)) 


gexpr6 = subset (gexpr, species=="Monkey")
Idents(gexpr6) <- "group"
recMarkers <- FindAllMarkers(gexpr6,only.pos= T, min.pct = 0.1, logfc.threshold = log(1.25),verbose = F,max.cells.per.ident = 2000)
write.xlsx(recMarkers,file = "2.fig3B_Monkey_cor.xlsx")
recMarkers <- subset(recMarkers[grep("^RP[L|S]",recMarkers$gene, ignore.case = FALSE,invert=TRUE),],subset=p_val_adj < 0.05)
region_av <- AverageExpression(gexpr6,group.by = "group",features = unique(recMarkers$gene),assays = "RNA") 
region_av <- region_av$RNA
gene_cell_exp <- t(scale(t(region_av),scale = T,center = T))

dat2_cor <- cor(as.matrix(gene_cell_exp))
p2=corrplot(dat2_cor,order = "original", type = "full", addrect = 4,tl.pos = "lt", tl.cex=.8,col = rev(brewer.pal(n=8, name="RdYlBu")),
            tl.col = "black",cl.cex = 0.7,cl.pos = "r",cl.length = 5,title="Monkey",
            cl.ratio=0.1,mar=c(3, 0, 2, 0)) 


dev.off()

####################################################################
##fig3c and d
library(clusterProfiler)
library(stringr)
library(org.Hs.eg.db)
# library(enrichplot)
library(GseaVis)
mydata1 <- read.xlsx("3.fig3CD_venn_organized.xlsx")
mydata1_H <- mydata1[mydata1$species=="Human",]
mydata1_M <- mydata1[mydata1$species=="Monkey",]

CTX_H=unique(mydata1_H$gene[which(mydata1_H$ident2=="CTX" & mydata1_H$cluster=="CTX")])
CGE_H=unique(mydata1_H$gene[which(mydata1_H$ident2=="CGE"& mydata1_H$cluster=="CGE")])
LGE_H=unique(mydata1_H$gene[which(mydata1_H$ident2=="LGE"& mydata1_H$cluster=="LGE")])

CTX_M=unique(mydata1_M$gene[which(mydata1_M$ident2=="CTX" & mydata1_M$cluster=="CTX")])
CGE_M=unique(mydata1_M$gene[which(mydata1_M$ident2=="CGE"& mydata1_M$cluster=="CGE")])
LGE_M=unique(mydata1_M$gene[which(mydata1_M$ident2=="LGE"& mydata1_M$cluster=="LGE")])

library(VennDiagram)
p=venn.diagram(
  x = list(CTX_H, CGE_H, LGE_H),
  category.names = c("CTX enriched v.s. MGE" , "CGE enriched v.s. MGE" , "LGE enriched v.s. MGE"),
  output=F,
  filename = NULL,
  
  lwd = 2, 
  lty = "blank",  
  fill = myCol,  
  
  cex = .5,  
  fontface = "bold", 
  fontfamily = "sans", 
  
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",  
  cat.pos = c(-27, 27, 180),  
  cat.dist = c(0.055, 0.055, 0.055),  
  cat.fontfamily = "sans",
  rotation = 1,  
  main = "Human",
  main.just = c(0.5, 1),
  main.fontface = "bold"
)

p1=venn.diagram(
  x = list(CTX_M, CGE_M, LGE_M),
  category.names = c("CTX enriched v.s. MGE" , "CGE enriched v.s. MGE" , "LGE enriched v.s. MGE"),
  output=F,
  filename = NULL,
  
  lwd = 2,
  lty = "blank",  
  fill = myCol,  
  
  cex = .5,  
  fontface = "bold",  
  fontfamily = "sans",  
  
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",  
  cat.pos = c(-27, 27, 180),  
  cat.dist = c(0.055, 0.055, 0.055),  
  cat.fontfamily = "sans",
  rotation = 1, 
  main = "Monkey",
  main.just = c(0.5, 1),
  main.fontface = "bold"
)
pdf("3.fig3CD_venn.pdf",6,3)
par(omi=c(0.1,0.1,0.1,0.1))
par(mfrow=c(1,2))
# grid.draw(p)
require(gridExtra)
grid.arrange(p, p1,ncol = 2, nrow = 1)
dev.off()









#####-------------------
library(clusterProfiler)
library(stringr)
library(org.Hs.eg.db)
# library(enrichplot)
library(GseaVis)

geneset <- read.gmt("c5.all.v2023.1.Hs.entrez.gmt")


##for human 
CTX_H=unique(mydata1_H$gene[which(mydata1_H$ident2=="CTX" & mydata1_H$cluster=="CTX")])
MGE_H=unique(mydata1_H$gene[which(mydata1_H$ident2=="CTX" & mydata1_H$cluster=="MGE")])

intersect(CTX_H,MGE_H)

marker1 <- mydata1_H[mydata1_H$ident2 %in% c("CTX") & mydata1_H$cluster=="MGE" & mydata1_H$gene %in% c(MGE_H),] 
marker2 <- mydata1_H[mydata1_H$ident2=="CTX" & mydata1_H$cluster=="CTX" & mydata1_H$gene %in% c(CTX_H),]
marker <- rbind(marker1,marker2)
marker$avg_log2FC[which(marker$cluster %in% c("MGE"))] <- -marker$avg_log2FC[which(marker$cluster %in% c("MGE"))]
table(marker$avg_log2FC>0)

gs <-bitr(marker$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
setdiff(marker$gene,gs$SYMBOL)
marker<-cbind(marker[match(gs[,1],marker$gene),],gs)
geneList = marker$avg_log2FC
names(geneList) = marker$ENTREZID
geneList = sort(geneList,decreasing = T)  


set.seed(1)

egmt_h <- GSEA(geneList, TERM2GENE=geneset,verbose=F,eps = 0,seed = 0)  
egmt_h<- setReadable(egmt_h,OrgDb=org.Hs.eg.db, keyType = "ENTREZID")  
y=data.frame(egmt_h)
y$Description = unlist(lapply(str_split(y$Description,pattern = "_",n=2,simplify = F), function(x) tolower(x[2])))
y$Description <- gsub("_"," ",x = y$Description)
egmt_h@result <- y
write.xlsx(y,file = "3.fig3CD_human_GSEA.xlsx")
paths <- c("GOCC_SOMATODENDRITIC_COMPARTMENT","GOBP_SYNAPTIC_SIGNALING","GOCC_SYNAPSE","GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING","GOCC_NEURON_PROJECTION",
           "GOBP_CELL_DIVISION","GOBP_FOREBRAIN_DEVELOPMENT","GOBP_CELL_CYCLE")
# p2 <- gseaplot2(egmt_h, geneSetID = paths, subplots=c(1,2), pvalue_table = F,rel_heights=c(1, .2),base_size = 11) 
p1=gseaNb(object = egmt_h,geneSetID = paths,curveCol = c("#8C5929","#E8B789","#E39042","#634E3B","#B37134",
                                                         # "#E08F80","#E67457",
                                                         "#CF5978","#B756C4","#8353DB"),
          # "#849C49","#516912","#CFC476","#B1B5D8","#49679C"
          subPlot = 1,
          termWidth = 35,lineSize = 0.8,base_size = 8
          # legend.position = c(0.8,0.8)
)+theme(aspect.ratio = 1/1)
# +labs(title = "Human")+theme(plot.title = element_text(hjust = 0.5,vjust = 0.5))





##for monkey 
CTX_M=unique(mydata1_M$gene[which(mydata1_M$ident2=="CTX" & mydata1_M$cluster=="CTX")])
MGE_M=unique(mydata1_M$gene[which(mydata1_M$ident2=="CTX" & mydata1_M$cluster=="MGE")])

intersect(CTX_M,MGE_M)

marker1 <- mydata1_M[mydata1_M$ident2 %in% c("CTX") & mydata1_M$cluster %in% c("MGE") & mydata1_M$gene %in% c(MGE_M),] 
marker2 <- mydata1_M[mydata1_M$ident2=="CTX" & mydata1_M$cluster=="CTX" & mydata1_M$gene %in% c(CTX_M),]
marker <- rbind(marker1,marker2)
marker$avg_log2FC[which(marker$cluster %in% c("MGE"))] <- -marker$avg_log2FC[which(marker$cluster %in% c("MGE"))]

gs <-bitr(marker$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
setdiff(marker$gene,gs$SYMBOL)
marker<-cbind(marker[match(gs[,1],marker$gene),],gs)
geneList = marker$avg_log2FC
names(geneList) = marker$ENTREZID
geneList = sort(geneList,decreasing = T)  

set.seed(1)
egmt1 <- GSEA(geneList, TERM2GENE=geneset,verbose=F,eps = 0,seed = 0)  
egmt1<- setReadable(egmt1,OrgDb=org.Hs.eg.db, keyType = "ENTREZID") 
y=data.frame(egmt1)
y$Description = unlist(lapply(str_split(y$Description,pattern = "_",n=2,simplify = F), function(x) tolower(x[2])))
y$Description <- gsub("_"," ",x = y$Description)
egmt1@result <- y
write.xlsx(y,file = "3.fig3CD_monkey_GSEA.xlsx")
paths <- c("GOBP_REGULATION_OF_MEMBRANE_POTENTIAL","GOCC_SYNAPSE","GOCC_SOMATODENDRITIC_COMPARTMENT","GOCC_POSTSYNAPSE","GOCC_DENDRITIC_TREE",
           "GOMF_RNA_POLYMERASE_II_TRANSCRIPTION_REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING")
# p2 <- gseaplot2(egmt1, geneSetID = paths, subplots=c(1,2), pvalue_table = F,rel_heights=c(1, .2),base_size = 11) # rel_heights 副图的相对高度
p2=gseaNb(object = egmt1,geneSetID = paths,curveCol = c("#8C5929","#E8B789","#E39042","#634E3B","#B37134",
                                                        "#E67457"),
          # "#849C49","#516912","#CFC476","#B1B5D8","#49679C"
          subPlot = 1,
          termWidth = 35,lineSize = 0.8,base_size = 8
          # legend.position = c(0.8,0.8)
)+theme(aspect.ratio = 1/1)
# +labs(title = "Human")+theme(plot.title = element_text(hjust = 0.5,vjust = 0.5))


save(egmt1,egmt_h,file = "3.fig3CD_GSEA.RData")

pdf('3.fig3CD_GSEA.pdf',width=5,height=4)
par(omi=c(0.1,0.1,0.1,0.1))
p1
p2
dev.off()



####################################################################
##fig3e
set.seed(1234)
sample_recMarkers <- NULL
sample_recMarkers <- read.xlsx("6.fig3G_volcano.xlsx")
sample_recMarkers$avg_log2FC[which(sample_recMarkers$cluster=="Monkey")]= -(sample_recMarkers$avg_log2FC[which(sample_recMarkers$cluster=="Monkey")])
sample_recMarkers$p_val_adj[which(sample_recMarkers$p_val_adj==0)] <- 1e-323  
mydata1 = sample_recMarkers
mydata1$log10_p_val_adj=log10(mydata1$p_val_adj)
## set label and up down ns 
cut_off_pvalue = 0.01
cut_off_logFC = 1
mydata1$threshold="NS";
mydata1[which(mydata1$avg_log2FC  > cut_off_logFC & -mydata1$log10_p_val_adj > (-log10(cut_off_pvalue))),]$threshold="Human up";
mydata1[which(mydata1$avg_log2FC  < (-cut_off_logFC) & -mydata1$log10_p_val_adj > (-log10(cut_off_pvalue))),]$threshold="Macaque up";
mydata1$threshold=factor(mydata1$threshold, levels=c('Human up','NS','Macaque up'))


library(ggrepel)
p_list=NULL
n=0
for (i in 1:length(c("MAF","NKX2_1"))) {
  group.use=c("MAF","NKX2_1")[i]
  cat(group.use,"\n")
  mydata2=mydata1[mydata1$group==group.use,]
  lobeRef = unique(mydata2$lobe)
  for (j in 1:length(lobeRef)) {
    lobe.use = lobeRef[j]
    mydata3=mydata2[mydata2$lobe==lobe.use,]
    mydata_label=mydata3[mydata3$threshold %in% c("Human up","Macaque up")&abs(mydata3$avg_log2FC)>2,]
    mydata_label=mydata_label %>% arrange(desc(avg_log2FC),desc(log10_p_val_adj)) %>% group_by(threshold) %>% top_n(n = 8,wt = abs(avg_log2FC))
    mydata_label=rbind(mydata3[mydata3$gene %in% c("FOS","COX7C"),],mydata_label)
    p <- ggplot(mydata3,aes(avg_log2FC, -log10_p_val_adj))+
      geom_hline(yintercept = 2, linetype = "dashed", color = "#999999")+
      geom_vline(xintercept = c(-cut_off_logFC,cut_off_logFC), linetype = "dashed", color = "#999999")+
      geom_point(aes(size=-log10_p_val_adj,color=threshold))+
      scale_color_manual(name="",labels=c("Human up",'NS',"Macaque up"),
                         values=c("#fc4e07", "grey","#00afbb" ) )+
      scale_size_continuous(range = c(0.1,3))+
      theme_bw()+
      theme(panel.grid = element_blank())+
      theme(aspect.ratio = 1/1)+
      scale_x_continuous(limits = c(-ceiling(max(abs(mydata3$avg_log2FC))),ceiling(max(abs(mydata3$avg_log2FC)))))+
      ylim(0,ceiling(max(-mydata3$log10_p_val_adj*1.2)))+
      guides (color = guide_legend (override.aes = list(shape = 15)))+
      geom_text_repel(data = subset(mydata_label,avg_log2FC>0),mapping = aes(label = gene), size = 3,fontface="bold",color="#fc4e07",
                      max.overlaps = 20,force = 1.5,seed = 1234)+
      geom_text_repel(data = subset(mydata_label,avg_log2FC<0),mapping = aes(label = gene), size = 3,fontface="bold",color="#00afbb",
                      max.overlaps = 20,force = 1.5,seed = 1234)+
      xlab("Log2FC")+
      ylab("-Log10 P value adujusted)")+labs(title = paste(group.use,lobe.use,sep = "_"))+theme(plot.title = element_text(hjust = 0.5))+
      annotate('text',x=3,y=ceiling(max(abs(mydata3$log10_p_val_adj))*1.17),label=paste0("Human ",table(mydata3$threshold)[[1]]),size=3,color='black')+
      annotate('text',x=-3,y=ceiling(max(abs(mydata3$log10_p_val_adj))*1.17),label=paste0("Macaque ",table(mydata3$threshold)[[3]]),size=3,color='black')+
      annotate("segment", x=-1,y=ceiling(max(abs(mydata3$log10_p_val_adj))*1.13), xend = -5, yend = ceiling(max(abs(mydata3$log10_p_val_adj))*1.13),
               colour = "#00afbb", linewidth = 1.5, arrow = arrow(length = unit(.2,"cm"),angle = 30))+
      annotate("segment", x=1,y=ceiling(max(abs(mydata3$log10_p_val_adj))*1.13), xend = 5, yend = ceiling(max(abs(mydata3$log10_p_val_adj))*1.13),
               colour = "#fc4e07", linewidth = 1.5, arrow = arrow(length = unit(.2,"cm"),angle = 30))
    
    n=n+1
    p_list[[n]]=p
    
  }
}


pdf("6.fig3E_volcano.pdf",20,10)
par(mar=c(0.1,0.1,0.1,0.1))
patchwork::wrap_plots(p_list, nrow = 2, ncol = 4)
dev.off()


####################################################################
##fig 3g
cellex <- read.csv("MAF_cellex_data.csv",row.names = 1)
species_transAgeRef=c("Human_GW18","Human_GW20","Monkey_GW16","Monkey_GW20","Monkey_GW21",
                      "Monkey_GW26","Monkey_GW33")

result_final <- NULL
for (i in 1:length(species_transAgeRef)) {
  species_transAge.use=species_transAgeRef[i]
  mydata1=cellex[cellex$species_transAge==species_transAge.use,]
  for (j in 1:100) {
    mydata2=mydata1[mydata1$order==j,]
    result <- as.data.frame(t(apply(mydata2[,1:4], 2, function(x) length(which(x>0.2)))))
    result$order=j
    result$species_transAge=species_transAge.use
    result_final <- rbind(result_final,result)
  }
  
}

result_3=NULL
for (i in 1:length(species_transAgeRef)) {
  species_transAge.use=species_transAgeRef[i]
  result_1=result_final[result_final$species_transAge==species_transAge.use,]
  result_2=as.data.frame(sapply(result_1[,1:4], function(x) c("Stand dev" = sd(x),
                                                              "Mean"= mean(x,na.rm=TRUE))))
  result_2=rownames_to_column(result_2,var = "values")
  result_2$species_transAge=species_transAge.use
  result_3=rbind(result_3,result_2)
}
result_3[,2:5] <- round(result_3[,2:5],1)

dat_MAF_1 <- result_3[result_3$species_transAge %in% c("Monkey_GW16","Monkey_GW20","Monkey_GW21",
                                                       "Monkey_GW26","Monkey_GW33") & result_3$values=="Mean",]
dat_MAF_1 <- dat_MAF_1[,-1]
dat_MAF_1 <- dat_MAF_1[,c(5,1:4)]
dat_MAF_1$species_transAge[dat_MAF_1$species_transAge %in% c("Monkey_GW20","Monkey_GW21")] <- "Monkey_GW20&GW21"
dat_MAF_1 = dat_MAF_1 %>% group_by(species_transAge) %>% mutate(FC=mean(FC),MSC=mean(MSC),TC=mean(TC),OC=mean(OC)) %>% distinct(species_transAge,.keep_all = T)

dat_MAF_1$species_transAge <- factor(dat_MAF_1$species_transAge,levels = c("Monkey_GW16","Monkey_GW20&GW21",
                                                                           "Monkey_GW26","Monkey_GW33"))
dat_MAF_1 <- dat_MAF_1[order(dat_MAF_1$species_transAge),]

library(ggradar)
p2 <- ggradar(dat_MAF_1,
              centre.y = 0,
              grid.min = 0,
              grid.mid = 100,
              grid.max = 500,
              values.radar = c("0", "100", "500"),
              background.circle.colour = "white",
              axis.label.size=4,
              legend.position = "bottom",
              group.line.width = 0.8,
              group.point.size = 0,
              grid.label.size = 3,
              gridline.min.colour = 'grey80',
              gridline.mid.colour = 'grey80',
              gridline.max.colour = 'blue',
              # gridline.min.linetype = 'solid',
              legend.text.size=8,
              plot.title = "Monkey"
)+theme(plot.title = element_text(hjust = 0.5,size = 12))



dat_MAF_2 <- result_3[result_3$species_transAge %in% c("Human_GW18","Human_GW20") & result_3$values=="Mean",]
dat_MAF_2 <- dat_MAF_2[,-1]
dat_MAF_2 <- dat_MAF_2[,c(5,1:4)]
dat_MAF_2$species_transAge <- factor(dat_MAF_2$species_transAge,levels = c("Human_GW18","Human_GW20"))
dat_MAF_2 <- dat_MAF_2[order(dat_MAF_2$species_transAge),]

p1 <- ggradar(dat_MAF_2,
              centre.y = 0,
              grid.min = 50,
              grid.mid = 100,
              grid.max = 400,
              values.radar = c("50", "100", "400"),
              background.circle.colour = "white",
              axis.label.size=4,
              legend.position = "bottom",
              group.line.width = 0.8,
              group.point.size = 0,
              grid.label.size = 3,
              gridline.min.colour = 'grey80',
              gridline.mid.colour = 'grey80',
              gridline.max.colour = 'blue',
              gridline.min.linetype = 'solid',
              legend.text.size=8,
              plot.title = "Human",
              group.colours=c("#f8aea7","#3b9ab0")
)+theme(plot.title = element_text(hjust = 0.5,size = 12))

pdf("5.fig3F_MAF_radar.pdf",8,4)
par(omi=c(0.1,0.1,0.1,0.1))
p1+p2+plot_layout(ncol=2,guides = "collect")&theme(legend.position = "bottom")
dev.off()




####################################################################
##fig 3f
mat=read.xlsx("4.fig3E_Shnnon_Entropy.CTX_mgeLin.xlsx")
mydata <- melt(mat,id=c("transAge","species_region"))
colnames(mydata) <- c("transAge","species_region","region","entropy")
mydata <- mydata[!is.na(mydata$entropy),]
library(tidyr)
library(tidyverse)
mydata <- separate(mydata,col = species_region,into=c("species","Cortex"),remove = F)

mydata$transAge <- factor(mydata$transAge,levels = ageRef)
mydata$species <- factor(mydata$species,levels = c("Human","Monkey"))
mydata <- mydata[order(mydata$species,mydata$transAge),]


p_m=ggplot(data=mydata, aes(transAge, entropy,fill=species,color="color"),show.legend = F) +
  labs(x="Age", y="Entropy") +
  geom_point(shape=21) +
  # ylim(0,3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title=element_text(hjust = 0.5))+
  geom_smooth(level=0.3,aes(group=species,color=species,fill=species,alpha=0.5),linewidth=1.5,span=1)+
  scale_color_manual(name = "Species",values=c("#ea218e","#00aeee"),limits=c("Human","Monkey"))+
  scale_fill_manual(name = "Species",values=c("#ea218e","#00aeee"),limits=c("Human","Monkey"))+
  guides(alpha='none')

pdf('4.fig3E_shannon_update.pdf',6,3)
par(omi=c(1,0.1,0.1,0.1))
p_m
dev.off()


meta.use=FetchData(gexpr.m_h,vars = c("transAge","subtype"))
df_1=meta.use %>% group_by(transAge,subtype) %>% mutate(count1=n()) %>% unique()
p=ggplot(df_1, aes(x=transAge, y=count1,fill=subtype)) + 
  geom_bar(stat="identity", position = "fill",width = 0.9)+
  scale_fill_manual(values = cols[unique(df_1$subtype)])+
  xlab("")+
  ylab("")+
  # ggtitle(species.use)+
  theme(panel.background=element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
  )+
  theme(legend.position = "right",legend.title=element_text(size=12),legend.text=element_text(size=8))
pdf('4.fig3E_stacked_legend.pdf',12,10)
par(omi=c(1,0.1,0.1,0.1))
p
dev.off()





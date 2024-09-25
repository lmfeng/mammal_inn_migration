library(circlize)
library(dplyr)
library(tibble)
library(ggsci)


#####################################-----------------------
####
gexpr <- readRDS("final_data/mammal_InN_migration.final.p8_01.rds")
Marker=c("VIM","MKI67","LHX6","NR2F2","PROX1","MEIS2")
targetGeneExpr.df=data.frame(t(as.matrix(gexpr@assays$RNA@data[Marker,])))

all(rownames(targetGeneExpr.df)==rownames(gexpr@meta.data))

maxGene <- function(x){
  ifelse(max(x)==0,"Unsign",names(which.max(x)))   
}

MarkerGroup=apply(targetGeneExpr.df, 1, maxGene)
gexpr$MarkerGroup <- MarkerGroup
gexpr1 <- subset(gexpr,MarkerGroup!="Unsign")
gexpr1$MarkerGroup <- factor(gexpr1$MarkerGroup,levels = Marker)
myCol2=c("#845EC2","#D65DB1","#FF6F91","#FF9671","#FFC75F","#F9F871")
names(myCol2) <- Marker
# Idents(gexpr.m) <- gexpr.m$lineage
p2 <- DimPlot(gexpr1,reduction ="scANVI_umap",pt.size = 0.01,label = F,shuffle = T,group.by = "MarkerGroup",cols = myCol2,raster = T)+
  theme(legend.position = "right",aspect.ratio = 1/1,axis.line = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),legend.key.width=unit(3,'mm'),plot.title = element_text(size = 8))

pdf("featureplot_five_genes.pdf",6,4)
p2
dev.off()


gexpr <- readRDS("final_data/mammal_InN_migration.final.p8_01.rds")

gexpr$stage="11"
gexpr$stage[which(gexpr$transAge %in% c("GW07","GW08","GW09","GW10","GW11","GW12","GW13"))] <- "1st"
gexpr$stage[which(gexpr$transAge %in% c("GW14","GW15","GW16","GW18","GW19","GW20","GW21","GW22","GW25","GW26","GW27","GW2trimester"))] <- "2nd"
gexpr$stage[which(gexpr$transAge %in% c("GW29","GW30","GW33","GW3trimester"))] <- "3rd"

p2=DimPlot(gexpr,group.by = "lobe",reduction ="scANVI_umap",split.by = "stage",raster = T,cols = c("#f8766d","#db8e00","#aea200","#64b200","#00bd5c","#00c1a7","#00bade","#b385ff","#ef67eb","#ff63b6"),order = rev(c("GE","FC","MSC","TC","OC","Insula","NCX","HIP","AMY","STR")),
           shuffle = T,pt.size = 0.01) &
  theme(legend.position = "right",panel.border = element_blank(),aspect.ratio = 1/1,
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),axis.line = element_blank(),plot.title = element_blank())


pdf("fig1_splitbyStage_umap.pdf",8,6)
par(omi=c(0.1,0.1,0.1,0.1))
p2
dev.off()






#################################################################################################################
gexpr <- readRDS("final_data/mammal_InN_migration.final.p8_01.rds")

set.seed(1234)

p1=DimPlot(gexpr,group.by = "lineage",reduction ="scVI",raster = T,cols = c(alpha(c("#984ea3","#377eb8","#ff7f00","#e41a1c","#00A087FF"),0.8)),order = rev(c("RGC","inIPC","mgeLin","cgeLin","lgeLin")),
           shuffle = T,pt.size = 0.01) &
  theme(legend.position = "right",panel.border = element_blank(),aspect.ratio = 1/1,
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),axis.line = element_blank(),plot.title = element_blank())


pdf("umap_lineage.pdf",6,6)
par(omi=c(0.1,0.1,0.1,0.1))
p1
dev.off()

#################################################################################################################
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
dim(gexpr)

library(gplots)




mydata=FetchData(gexpr,vars = c("labsite","species","transAge"))
mydata=mydata %>% group_by(labsite,species,transAge) %>% mutate(counts=n()) %>% unique()
mydata=reshape2::dcast(mydata,labsite+species~transAge,value.var = "counts")


rownames(mydata) = as.character(mydata[,1])
mydata = mydata[, c(-1, -2)]
mydata=mydata[c("KriegsteinLab", "Kriegstein2022","WangLab","Linnarsson2022","RakicLab","PollenLab"),]
mydata=mydata[,c("GW07","GW08","GW09","GW10","GW11","GW12","GW13","GW14","GW15" , "GW16" , "GW18" ,"GW19" ,"GW20","GW21","GW22","GW25" , "GW26" ,"GW27",
                 "GW2trimester","GW29","GW30" ,"GW33","GW3trimester")]
mydata = as.matrix(log10(mydata))


###
res = mydata
##remove mouse data
res = res[1:6, ]

###----plot
pdf('mammal.age.match.pdf',7,4)
par(omi = c(0.1, 0.1, 0.1, 0.1))
pairs.breaks <- seq(min(res,na.rm=T), max(res,na.rm=T), length.out=101)
mycol <- colorpanel(n=100,low="white",high="brown");
heatmap.2(res, breaks=pairs.breaks,col=mycol,main="",Rowv=F,Colv=F,dendrogram="none"
          ,na.rm=TRUE,trace="none",na.color="white",density.info="none"
          ,key.xlab = NA
)
dev.off()






##fig4a
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
gexpr <- subset(gexpr,lineage=="cgeLin")
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


pdf("fig4A_Lobe_pie.pdf",6,6)
par(mar=c(0.1,0.1,0.1,0.1))
p1+p2
dev.off()

############################################################################################
##fig4b
library(Seurat)
library(dplyr)
library(ggplot2)
library(corrplot)
library(RColorBrewer)
library(openxlsx)
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
gexpr <- subset(gexpr,lineage=="cgeLin")
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

pdf("2.fig4B_cor.pdf",10,5)
par(omi=c(0.1,0.1,0.1,0.1))
par(mfrow=c(1,2))

set.seed(1234)
gexpr5 = subset (gexpr, species=="Human")
Idents(gexpr5) <- "group"
recMarkers <- FindAllMarkers(gexpr5,only.pos= T, min.pct = 0.1, logfc.threshold = log(1.25),verbose = F,max.cells.per.ident = 2000)
write.xlsx(recMarkers,file = "2.fig4B_Human_cor.xlsx")
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
write.xlsx(recMarkers,file = "2.fig4B_Monkey_cor.xlsx")
recMarkers <- subset(recMarkers[grep("^RP[L|S]",recMarkers$gene, ignore.case = FALSE,invert=TRUE),],subset=p_val_adj < 0.05)
region_av <- AverageExpression(gexpr6,group.by = "group",features = unique(recMarkers$gene),assays = "RNA") 
region_av <- region_av$RNA
gene_cell_exp <- t(scale(t(region_av),scale = T,center = T))
dat2_cor <- cor(as.matrix(gene_cell_exp))
p2=corrplot(dat2_cor,order = "original", type = "full", addrect = 4,tl.pos = "lt", tl.cex=.8,col = rev(brewer.pal(n=8, name="RdYlBu")),
            tl.col = "black",cl.cex = 0.7,cl.pos = "r",cl.length = 5,title="Monkey",
            cl.ratio=0.1,mar=c(3, 0, 2, 0)) 


dev.off()


############################################################################################
##fig4c
mat=read.xlsx("4.Shnnon_Entropy.CTX_cgeLin.xlsx")
mydata <- melt(mat,id=c("transAge","species_region"))
colnames(mydata) <- c("transAge","species_region","region","entropy")
mydata <- mydata[!is.na(mydata$entropy),]
library(tidyr)
mydata <- separate(mydata,col = species_region,into=c("species","Cortex"),remove = F)

mydata$transAge <- factor(mydata$transAge,levels = ageRef)
mydata$species <- factor(mydata$species,levels = c("Human","Monkey"))
mydata <- mydata[order(mydata$species,mydata$transAge),]


p_m <- ggplot(data=mydata, aes(transAge, entropy,fill=species,color="black"),show.legend = F) +
  labs(x="Age", y="Entropy") +
  geom_point(shape=21) +
  # ylim(0,3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title=element_text(hjust = 0.5))+
  geom_smooth(level=0.3,aes(group=species,color=species,fill=species,alpha=0.5),linewidth=1.5,span=1)+
  scale_color_manual(name = "Species",values=c("#ea218e","#00aeee"),limits=c("Human","Monkey"))+
  scale_fill_manual(name = "Species",values=c("#ea218e","#00aeee"),limits=c("Human","Monkey"))+
  guides(alpha='none')



############################################################################################
##fig4d
library(dplyr)
library(Seurat)
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
meta = gexpr@meta.data
meta.use <- NULL
rec <- NULL
res <- NULL
resPercent <- NULL
meta.use = subset(meta, lineage == 'cgeLin'& species == 'Human')
rec = table(meta.use$lobe, meta.use$subtype)
res = 100 * rec/apply(rec, 1, sum)
celltypeRef.h = colnames(res)
##
xdata <- NULL
ydata <- NULL
xdata = res['GE',]
ydata = (res['FC', ] + res['OC', ] + res['MSC', ] + res['TC', ] + res['Insula', ])/5
##
xdata.h = log2(xdata + 1)
ydata.h = log2(ydata + 1)


###Monkey
meta.use <- NULL
rec <- NULL
res <- NULL
meta.use = subset(meta, lineage == 'cgeLin'& species == 'Monkey')
rec = table(meta.use$lobe, meta.use$subtype)
res = 100 * rec/apply(rec, 1, sum)
celltypeRef.m = colnames(res)
##
xdata <- NULL
ydata <- NULL
xdata = res['GE',]
ydata = (res['FC', ] + res['OC', ] + res['MSC', ] + res['TC', ] + res['Insula', ])/5
##
xdata.m = log2(xdata + 1)
ydata.m = log2(ydata + 1)


pdf('7.fig4D_cgeLin.ge_vs_ctx.test.pdf')
par(omi = c(0.1, 0.1, 0.1, 0.1))

##Human
plot(xdata.h, ydata.h, type = 'p', pch = 16, xlim = c(0,8), ylim = c(0,8), col = '#ea218e', xlab = 'GE', ylab = 'CTX',)
text(xdata.h, ydata.h, labels = celltypeRef.h, col = '#ea218e',cex = 0.5,pos = 4)

##Monkey
points(xdata.m, ydata.m, type = 'p', pch = 15, col = '#00aeee')
text(xdata.m, ydata.m, labels = celltypeRef.m, col = '#00aeee',cex = 0.5,pos = 4)

dev.off()



############################################################################################
#fig4e
library(dplyr)
library(Seurat)
###feature plot umap colored by CRH expression 
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
pdf("6.fig4E_CRH_featureplot.pdf",10,5)
par(omi=c(0.1,0.1,0.1,0.1))
FeaturePlot(gexpr,features = "CRH",pt.size = 0.1,order = T,cols = c("lightgrey" ,"#DE1F1F"),raster = T,reduction = "scVI",split.by = "species")&
  theme(plot.title = element_text(hjust = 0.5),legend.position = "right",panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid"),aspect.ratio = 1/1)
dev.off()




############################################################################################
##fig4f
load("CRH_famaily_addLister.RData")
set.seed(1234)
gexpr=subset(gexpr,transAge %in% c("GW2trimester","GW3trimester"),invert=T)
gexpr=gexpr %>%  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) 


data.use=DotPlot(gexpr,features = rownames(gexpr),group.by = "age_species")
data.use=data.use$data
data.use=separate(data.use,col = id,into = c("transAge","species"),remove = F)


data.use1=data.use[,c("transAge","species","pct.exp","avg.exp.scaled","features.plot")]
data.use1$transAge=factor(data.use1$transAge,levels = c("GW07","GW08","GW09","GW10","GW11","GW12","GW13","GW14","GW15","GW16",
                                                        "GW17","GW18","GW19","GW20","GW21","GW22","GW24","GW25","GW26","GW27","GW29","GW30","GW33","GW34","juvenile","adult"))

data.use1$species=factor(data.use1$species,levels = rev(c("Human","Chimpanzee","Rhesus","Marmoset","Mouse")))
data.use1=data.use1[order(data.use1$transAge,data.use1$species),]

library(viridis)
pp_list=list()
for(i in 1:length(unique(data.use1$features.plot))) {
  gene.use=unique(data.use1$features.plot)[i]
  data.use2=data.use1[data.use1$features.plot==gene.use,]
  p <- ggplot(data.use2,aes(x=transAge,y=species,color=avg.exp.scaled,size=pct.exp,alpha = pct.exp == 0))+
    geom_point()+
    scale_color_gradientn(colours = viridis(20),
                          guide = guide_colorbar(ticks.colour = "black",frame.colour = "black",
                          ),name="Expression")+
    scale_size("% detected",range = c(0,10))+  
    scale_alpha_manual(values = c(1,0))+
    theme_classic()+
    theme(panel.grid = element_blank(),
          axis.text.x.bottom = element_text(hjust = 0, vjust = 1, angle = 90)
    )+
    scale_y_discrete(breaks=c("Human","Chimpanzee","Rhesus","Marmoset","Mouse"),labels=c("Human","Chimpanzee","Macaque","Marmoset","Mouse"))+
    labs(title = gene.use,x="",y="")
  pp_list[[i]]=p
}

pdf("5.fig4F_CRHfamily_expression.pdf",32,8)
wrap_plots(pp_list,ncol = 3)
dev.off()



############################################################################################
##fig4i,j,k
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
gexpr_1 <- subset(gexpr,lobe %in% c("FC","OC","MSC","TC","GE","Insula","NCX")&species=="Human"&lineage %in% c("cgeLin"))
gexpr_1$group <- "11"
gexpr_1$group[gexpr_1$lobe=="GE"]="GE"
gexpr_1$group[gexpr_1$lobe %in% c("FC","OC","MSC","TC","NCX","Insula")]="CTX"
table(gexpr_1$group)

 
Inscore <- AddModuleScore(gexpr_1,
                          features = list(Gpre),
                          ctrl = 100,
                          name = "CD_Features")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[ncol(Inscore@meta.data)] <- 'Gpre_Score' 
Inscore$group <- factor(Inscore$group,levels = c("CTX","GE"))


comparisons <- list(c("CTX","GE"))
p2=ggviolin(Inscore@meta.data, x="group", y="Gpre_Score", width = 0.8, color = "black",fill="group",palette = c("#e4d1a3", "#439595"), add = 'mean_sd',xlab = F, 
            bxp.errorbar=T,bxp.errorbar.width=0.5,size=0.3,          
            outlier.shape=NA, legend = "bottom",title = "Gpre_CTXvsGE")+
  geom_signif(comparisons = comparisons, 
              map_signif_level = T, 
              textsize =3, 
              test = wilcox.test,  #t.test
              step_increase = 0.2)

pdf("CRH_pathway/3.final.pdf",8,4)
p2
dev.off()



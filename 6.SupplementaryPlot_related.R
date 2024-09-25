##  difflabsite_cell&UMI&gene_barplot&violin in different samples  
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(tibble)
library(patchwork)
gexpr <- readRDS("final_data/mammal_InN_migration.final.p8_01.rds")
dim(gexpr)

#########################################################
###figS1a
regionRef=c("GE","MGE+CGE+LGE","CGE+LGE","MGE","CGE","LGE","xLGE","FR","FC","PFC","DFC","VFC","MFC","OFC",
            "MS","M1C","S1C","A1C","IPC","PCC", "Cing","ACC","TP","ITC","STC","OX","V1C",
            "Insula","Telencephalon","Cortex","HIP","AMY","STR","Putamen","Caudate","NAC")
lobeRef=c("GE","FC","MSC","TC","OC","Insula","NCX","HIP","AMY","STR")
cols=c("#f8766d","#db8e00","#aea200","#64b200","#00bd5c","#00c1a7","#00bade","#b385ff","#ef67eb","#ff63b6")
names(cols)=lobeRef

setdiff(regionRef,unique(gexpr$region))
setdiff(unique(gexpr$region),regionRef)

setdiff(lobeRef,unique(gexpr$lobe))
setdiff(unique(gexpr$lobe),lobeRef)

##for  human
gexpr_h=subset(gexpr,species=="Human")
mydata=FetchData(gexpr_h,vars = c("age","region","lobe"))
mydata=mydata %>% group_by(age,region) %>% mutate(counts=n())%>% unique()
mydata$stage=NA
mydata$stage[which(mydata$age %in% c("GW07","GW08","GW09","GW10","GW11","GW12","GW13"))] <- "1st"
mydata$stage[which(mydata$age %in% c("GW14","GW15","GW16","GW18","GW19","GW20","GW21","GW22","GW25","GW26","GW27","GW2trimester"))] <- "2nt"
mydata$stage[which(mydata$age %in% c("GW29","GW30","GW33","GW3trimester"))] <- "3rt"

mydata$age <- factor(mydata$age,levels = rev(c("GW07","GW08","GW09","GW10","GW12","GW13","GW14","GW16","GW18","GW19","GW20","GW22","GW25","GW27","GW2trimester","GW30","GW3trimester")))
mydata$lobe <- factor(mydata$lobe,levels = lobeRef)
mydata$region <- factor(mydata$region,levels = regionRef[sort(match(unique(mydata$region),regionRef))])

mydata=mydata[order(mydata$age,mydata$lobe,mydata$region),]




p1=ggplot(mydata,aes(x=region,y=age,size=counts,color=lobe,fill=lobe))+
  geom_point()+
  theme_bw()+
  scale_size_continuous(name="Cells",breaks=c(100,500,1000,2500,5000,10000,20000,30000),limits = c(5,30110),
                        labels = c(100,500,'1,000','2,500','5,000','10,000','20,000','30,000'),
                        range = c (0.5, 8))+
  scale_color_manual(name="Region",values = cols)+
  scale_fill_manual(values = cols)+
  # scale_shape_manual(name="Labsite",limits=c("KriegsteinLab","Kriegstein2022","WangLab","Linnarsson2022"),labels=c("Kriegstein 2021","Kriegstein 2022 snRNA-seq","Wang 2021","Linnarsson 2022"),values = c(0,1,2,5))+
  facet_grid(stage~lobe,scales ="free",space = "free")+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5),
        panel.border = element_blank())+
  labs(x=NULL,y=NULL,title = NULL)+guides(fill="none")



##for  monkey
gexpr_m=subset(gexpr,species=="Monkey")
mydata1=FetchData(gexpr_m,vars = c("region","lobe","transAge"))
mydata1=mydata1 %>% group_by(transAge,region) %>% mutate(counts=n())%>% unique()
mydata1$stage=NA
mydata1$stage[which(mydata1$transAge %in% c("GW07","GW08","GW09","GW10","GW11","GW12","GW13"))] <- "1st"
mydata1$stage[which(mydata1$transAge %in% c("GW14","GW15","GW16","GW18","GW19","GW20","GW21","GW22","GW25","GW26","GW27","GW2trimester"))] <- "2nt"
mydata1$stage[which(mydata1$transAge %in% c("GW29","GW30","GW33","GW3trimester"))] <- "3rt"

# mydata1$age <- factor(mydata1$age,levels = rev(c("E37","E40","E42","E43","E50","E54","E62","E64","E65","E77","E78","E80","E90","E93","E100","E110")))
mydata1$transAge <- factor(mydata1$transAge,levels = rev(sort(unique(mydata1$transAge))))
mydata1$lobe <- factor(mydata1$lobe,levels = lobeRef)
mydata1$region <- factor(mydata1$region,levels = regionRef[sort(match(unique(mydata1$region),regionRef))])

mydata1=mydata1[order(mydata1$transAge,mydata1$lobe,mydata1$region),]


p2=ggplot(mydata1,aes(x=region,y=transAge,size=counts,color=lobe,fill=lobe))+
  geom_point()+
  theme_bw()+
  scale_size_continuous(name="Cells",breaks=c(100,500,1000,2500,5000,10000,20000,30000),limits = c(5,30110),
                        labels = c(100,500,'1,000','2,500','5,000','10,000','20,000','30,000'),
                        range = c (0.5, 8))+
  scale_color_manual(name="Region",values = cols)+
  scale_fill_manual(values = cols)+
  facet_grid(stage~lobe,scales ="free",space = "free" )+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5),
        panel.border = element_blank())+
  labs(x=NULL,y=NULL,title = NULL)+guides(fill="none")



#############################################
###figS1b

##for human
mydata <- FetchData(gexpr_h,vars = c("labsite","transAge"))
mydata <- mydata %>% group_by(transAge,labsite) %>% mutate(cell_number=n())
mydata <- unique(mydata)

mydata$transAge <- factor(mydata$transAge,levels = sort(unique(mydata$transAge)))
mydata <- mydata[order(mydata$transAge),]
mydata <- mydata %>% group_by(transAge) %>% mutate(sum_cellnumber=sum(cell_number))

p3 <- ggplot(mydata, aes(x=transAge, y=cell_number, fill = labsite)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(name="Labsite",values = c("#ffc4cd","#85a6d8","#f9a255","#A9D9AB"),limits=c("KriegsteinLab","Kriegstein2022","WangLab","Linnarsson2022"),labels=c("Kriegstein 2021","Kriegstein 2022 snRNA-seq","Wang 2021","Linnarsson 2022"))+
  geom_text(aes(label=sum_cellnumber,y=sum_cellnumber),size=3,vjust=-0.5)+
  theme_classic() +
  # theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
  #       axis.ticks.x = element_blank(),axis.line.x = element_blank(),legend.position = "none")+
  ylab("# Cells")+xlab("Age")




### for monkey
mydata2 <- FetchData(gexpr_m,vars = c("labsite","transAge"))
mydata2 <- mydata2 %>% group_by(transAge,labsite) %>% mutate(cell_number=n())
mydata2 <- unique(mydata2)
mydata2$transAge <- factor(mydata2$transAge,levels = sort(unique(mydata2$transAge)))
mydata2 <- mydata2[order(mydata2$transAge),]
mydata2 <- mydata2 %>% group_by(transAge) %>% mutate(sum_cellnumber=sum(cell_number))

p4 <- ggplot(mydata2, aes(x=transAge, y=cell_number, fill = labsite)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(name="Labsite",values = c("#B6A1DE","#ffdab9"),limits=c("RakicLab","PollenLab"),labels=c("This study","Pollen 2022"))+
  geom_text(aes(label=sum_cellnumber,y=sum_cellnumber),size=3,vjust=-0.5)+
  theme_classic() +
  # theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
  #       axis.ticks.x = element_blank(),axis.line.x = element_blank(),legend.position = "none",
  #       axis.title.y = element_blank())+
  ylab("# Cells")+xlab("Age")



pdf("FigS1ab.pdf",20,8)
par(omi=c(0.1,0.1,0.1,0.1))
p1+p2+p3+p4+plot_layout(ncol = 2,heights = c(3,1))
dev.off()




gexpr <- readRDS("final_data/mammal_InN_migration.final.p8_01.rds")
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(dplyr)
library(patchwork)


set.seed(1234)
# gexpr@reductions$scVI_umap@key <- "scVI_UMAP_"
# colnames(gexpr@reductions$scVI_umap@cell.embeddings) <- c("scVI_UMAP_1","scVI_UMAP_2")

p1=DimPlot(gexpr,group.by = "lineage",reduction ="scANVI_umap",split.by = "species",raster = T,cols = c(alpha(c("#984ea3","#377eb8","#ff7f00","#e41a1c","#00A087FF"),0.8)),order = rev(c("RGC","inIPC","mgeLin","cgeLin","lgeLin")),
           shuffle = T,pt.size = 0.01) &
  theme(legend.position = "right",panel.border = element_blank(),aspect.ratio = 1/1,
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),axis.line = element_blank(),plot.title = element_blank())

p2=DimPlot(gexpr,group.by = "lobe",reduction ="scANVI_umap",split.by = "species",raster = T,cols = c("#f8766d","#db8e00","#aea200","#64b200","#00bd5c","#00c1a7","#00bade","#b385ff","#ef67eb","#ff63b6"),order = rev(c("GE","FC","MSC","TC","OC","Insula","NCX","HIP","AMY","STR")),
           shuffle = T,pt.size = 0.01) &
  theme(legend.position = "right",panel.border = element_blank(),aspect.ratio = 1/1,
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),axis.line = element_blank(),plot.title = element_blank())
p3=DimPlot(gexpr,group.by = "transAge",reduction ="scANVI_umap",split.by = "species",raster = T,cols = rev(colorRampPalette(brewer.pal(n = 11,"RdYlBu"))(23)),
           shuffle = T,pt.size = 0.01) &
  theme(legend.position = "right",panel.border = element_blank(),aspect.ratio = 1/1,
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),axis.line = element_blank(),plot.title = element_blank())
p4=DimPlot(gexpr,group.by = "labsite",reduction ="scANVI_umap",split.by = "species",raster = T,cols = c("#ffc4cd","#85a6d8","#f9a255","#A9D9AB","#B6A1DE","#ffdab9"),
           order = rev(c("KriegsteinLab", "Kriegstein2022","WangLab","Linnarsson2022","RakicLab","PollenLab")),
           shuffle = T,pt.size = 0.01) &
  theme(legend.position = "right",panel.border = element_blank(),aspect.ratio = 1/1,
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),axis.line = element_blank(),plot.title = element_blank())

pdf("figS1c.umap.pdf",12,8)
par(omi=c(0.1,0.1,0.1,0.1))
p4+p2+p3+p1+plot_layout(ncol = 2)
dev.off()





##---------------------------------------------------------
##for tranAge legend
transAgeRef=c("GW07","GW08","GW09","GW10","GW11","GW12","GW13","GW14","GW15","GW16","GW18","GW19","GW20","GW21","GW22","GW25","GW26","GW27","GW2trimester","GW29","GW30","GW33","GW3trimester")
setdiff(unique(gexpr$transAge),transAgeRef)
setdiff(transAgeRef,unique(gexpr$transAge))
gexpr$transAge=factor(gexpr$transAge,levels = transAgeRef)
gexpr$group=as.numeric(factor(gexpr$transAge))
scVI_umap = gexpr@reductions$scANVI_umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(group = gexpr@meta.data$group)



pdf("UMAPtransAge_legend.pdf",6,6)
ggplot(scVI_umap,aes(scANVIUMAP_1,scANVIUMAP_2,color=group))+
  geom_point(size = 0.01 , alpha =1)+
  scale_color_gradientn(colours = rev(colorRampPalette(brewer.pal(n = 11,"RdYlBu"))(23)),limits=c(1,23),breaks=c(1,5,10,15,18,22),labels=c("GW07","GW11","GW16","GW22","GW27","GW33"),
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name="transAge")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))
dev.off()




##---------------------------------------------------------
##for metaNeighbor
library(Seurat)
library(SingleCellExperiment)
# devtools::install_git('https://github.com/gillislab/MetaNeighbor')
library(MetaNeighbor)


#####metaNeighbor######
set.seed(1234)
gexpr.org <- readRDS("mammal_InN_migration.final.p8_01.rds")
dim(gexpr.org)


labsiteRef=c("PollenLab","RakicLab","KriegsteinLab")

for (i in 1:length(labsiteRef)) {
  labsite.use=labsiteRef[i]
  gexpr_Pollen=subset(gexpr.org,labsite==labsite.use)
  dim(gexpr_Pollen)
  
  ######---------------------------Pollen
  ###
  ####pollen
  load("orig_Pollen.RData")
  length(intersect(colnames(gexpr_Pollen),colnames(gexpr)))
  gexpr=gexpr[,intersect(colnames(gexpr_Pollen),colnames(gexpr))]
  dim(gexpr)
  gexpr_Pollen=gexpr_Pollen[,intersect(colnames(gexpr_Pollen),colnames(gexpr))]
  dim(gexpr_Pollen)
  gexpr$labsite="originaldata"
  gexpr$barcode=colnames(gexpr)
  
  gexpr$lineage="mgeLin"
  gexpr$lineage[gexpr$barcode %in% gexpr_Pollen$barcode[gexpr_Pollen$lineage=="mgeLin"]]="mgeLin"
  gexpr$lineage[gexpr$barcode %in% gexpr_Pollen$barcode[gexpr_Pollen$lineage=="cgeLin"]]="cgeLin"
  gexpr$lineage[gexpr$barcode %in% gexpr_Pollen$barcode[gexpr_Pollen$lineage=="lgeLin"]]="lgeLin"
  gexpr$lineage[gexpr$barcode %in% gexpr_Pollen$barcode[gexpr_Pollen$lineage=="inIPC"]]="inIPC"
  gexpr$lineage[gexpr$barcode %in% gexpr_Pollen$barcode[gexpr_Pollen$lineage=="RGC"]]="RGC"
  unique(gexpr$lineage)
  
  aa=intersect(rownames(gexpr),rownames(gexpr_Pollen))
  
  gexpr=gexpr[aa,]
  gexpr_Pollen=gexpr_Pollen[aa,]
  
  
  ##intersect gene
  Idents(gexpr_Pollen) <- "subtype"
  recMarkers <- FindAllMarkers(gexpr_Pollen, only.pos= T, min.pct = 0.1, logfc.threshold = log(1.25)
                               ,max.cells.per.ident = 2000,verbose = T)
  
  write.xlsx(recMarkers,file=paste0(labsite.use,"_subtypeFindAllMarkers.xlsx"))
  Idents(gexpr_Pollen) <- "lineage"
  recMarkers <- FindAllMarkers(gexpr_Pollen, only.pos= T, min.pct = 0.1, logfc.threshold = log(1.25)
                               ,max.cells.per.ident = 5000,verbose = T)
  write.xlsx(recMarkers,file=paste0(labsite.use,"_lineageFindAllMarkers.xlsx"))
  
  idx.pk = which(recMarkers$p_val_adj < 0.01)
  recMarkers = recMarkers[idx.pk, ]
  
  recMarkers1=recMarkers %>% group_by(cluster) %>% top_n(n = 30,wt=avg_log2FC)
  
  
  refdata=gexpr_Pollen
  refsce=as.SingleCellExperiment(refdata)
  refsce$study_id <- "ref"
  
  sce=as.SingleCellExperiment(gexpr)
  sce$study_id="test"
  hvgs = variableGenes(refsce, exp_labels = refsce$study_id)[1:1000]
  #In our experience, we obtained best performance for gene sets ranging from 200 to 1,000 variable genes
  pretrained_model <- MetaNeighbor::trainModel(var_genes = unique(c(hvgs,recMarkers1$gene)),
                                               dat = refsce,
                                               study_id = refsce$study_id, 
                                               cell_type = refsce$subtype 
  )
  
  aurocs <- MetaNeighborUS(trained_model=pretrained_model,
                           dat = sce,
                           study_id = sce$study_id,
                           cell_type = sce$subtype,
                           fast_version=TRUE)
  colnames(aurocs)=sapply(colnames(aurocs),function(x) strsplit(x,fixed = T,split = "|")[[1]][2])
  rownames(aurocs)=sapply(rownames(aurocs),function(x) strsplit(x,fixed = T,split = "|")[[1]][2])
  aurocs_1=t(aurocs)
  rownames(aurocs_1)=subtypeRef[sort(match(rownames(aurocs_1),subtypeRef))]
  
  pdf(paste0("1.",labsite.use,"_metaNeghbor.pdf"),12,12)
  plotHeatmapPretrained(aurocs_1,margins = c(15,15))
  
  
  top_hits = topHits(aurocs_1,
                     
                     dat = sce,
                     
                     study_id = sce$study_id,
                     
                     cell_type = sce$subtype,threshold = 0.8)
  
  best_hits <- MetaNeighborUS(trained_model=pretrained_model,
                              dat = sce,
                              study_id = sce$study_id,
                              cell_type = sce$subtype,
                              fast_version=TRUE,one_vs_best=TRUE,symmetric_output=FALSE)
  colnames(best_hits)=sapply(colnames(best_hits),function(x) strsplit(x,fixed = T,split = "|")[[1]][2])
  rownames(best_hits)=sapply(rownames(best_hits),function(x) strsplit(x,fixed = T,split = "|")[[1]][2])
  best_hits_1=t(best_hits)
  rownames(best_hits_1)=subtypeRef[sort(match(rownames(best_hits_1),subtypeRef))]
  plotHeatmapPretrained(best_hits_1,cex=0.5)
  
  dev.off()
}




##---------------------------------------------------------
##for venn plot 
#####data pre#####
sample_recMarkers <- read.xlsx("2.lineage_allmarkers.xlsx")

sample_recMarkers_1=sample_recMarkers %>% group_by(cluster,gene,species) %>% mutate(counts=length(unique(order)))
sample_recMarkers_1=sample_recMarkers_1[sample_recMarkers_1$counts==10,]
table(sample_recMarkers_1$counts)
sample_recMarkers_1=sample_recMarkers_1 %>% distinct(cluster,gene,species,counts,.keep_all = T)

human_cells_markers <- sample_recMarkers_1[sample_recMarkers_1$species=="Human",]
monkey_cells_markers <- sample_recMarkers_1[sample_recMarkers_1$species=="Monkey",]
mouse_cells_markers <- sample_recMarkers_1[sample_recMarkers_1$species=="Mouse",]

subclasses <- c("mgeLin","cgeLin","lgeLin")
tmp <- 1 #change this value to generate different subclass Venn diagrams
human_genes <- human_cells_markers[grep(subclasses[tmp], human_cells_markers$cluster), ]
monkey_genes <- monkey_cells_markers[grep(subclasses[tmp], monkey_cells_markers$cluster), ]
mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
# human_genes <- human_genes[which(human_genes$avg_diff > 0), ]
# monkey_genes <- monkey_genes[which(monkey_genes$avg_diff > 0), ]
# mouse_genes <- mouse_genes[which(mouse_genes$avg_diff > 0), ]


#venn Diagrams
all.genes <- data.frame(genes = unique(c(human_genes$gene, monkey_genes$gene,mouse_genes$gene)))
all.genes$Human <- as.character(match(all.genes$genes, human_genes$gene))
all.genes$Monkey <- as.character(match(all.genes$genes, monkey_genes$gene))
all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
# all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
all.genes$Human[which(is.na(all.genes$Human))] <- FALSE 
all.genes$Human[which(all.genes$Human != FALSE)] <- TRUE
all.genes$Monkey[which(is.na(all.genes$Monkey))] <- FALSE 
all.genes$Monkey[which(all.genes$Monkey != FALSE)] <- TRUE
all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE
all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
all.genes$Human <- as.logical(all.genes$Human)
all.genes$Monkey <- as.logical(all.genes$Monkey)
all.genes$Mouse <- as.logical(all.genes$Mouse)

library(eulerr)
p1=plot(euler(
  all.genes[ ,2:4]),
  quantities = list(cex = 3),
  labels = NULL,
  main = paste0(subclasses[tmp], " vs. All Inh"),
  fills = c("#ea218e", "#00aeee","#ea7b3d")
)


tmp <- 2 #change this value to generate different subclass Venn diagrams
human_genes <- human_cells_markers[grep(subclasses[tmp], human_cells_markers$cluster), ]
monkey_genes <- monkey_cells_markers[grep(subclasses[tmp], monkey_cells_markers$cluster), ]
mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
# human_genes <- human_genes[which(human_genes$avg_diff > 0), ]
# monkey_genes <- monkey_genes[which(monkey_genes$avg_diff > 0), ]
# mouse_genes <- mouse_genes[which(mouse_genes$avg_diff > 0), ]


#venn Diagrams
all.genes <- data.frame(genes = unique(c(human_genes$gene, monkey_genes$gene,mouse_genes$gene)))
all.genes$Human <- as.character(match(all.genes$genes, human_genes$gene))
all.genes$Monkey <- as.character(match(all.genes$genes, monkey_genes$gene))
all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
# all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
all.genes$Human[which(is.na(all.genes$Human))] <- FALSE 
all.genes$Human[which(all.genes$Human != FALSE)] <- TRUE
all.genes$Monkey[which(is.na(all.genes$Monkey))] <- FALSE 
all.genes$Monkey[which(all.genes$Monkey != FALSE)] <- TRUE
all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE
all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
all.genes$Human <- as.logical(all.genes$Human)
all.genes$Monkey <- as.logical(all.genes$Monkey)
all.genes$Mouse <- as.logical(all.genes$Mouse)

library(eulerr)
p2=plot(euler(
  all.genes[ ,2:4]),
  quantities = list(cex = 3),
  labels = NULL,
  main = paste0(subclasses[tmp], " vs. All Inh"),
  fills = c("#ea218e", "#00aeee","#ea7b3d")
)


tmp <- 3 #change this value to generate different subclass Venn diagrams
human_genes <- human_cells_markers[grep(subclasses[tmp], human_cells_markers$cluster), ]
monkey_genes <- monkey_cells_markers[grep(subclasses[tmp], monkey_cells_markers$cluster), ]
mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
# human_genes <- human_genes[which(human_genes$avg_diff > 0), ]
# monkey_genes <- monkey_genes[which(monkey_genes$avg_diff > 0), ]
# mouse_genes <- mouse_genes[which(mouse_genes$avg_diff > 0), ]


#venn Diagrams
all.genes <- data.frame(genes = unique(c(human_genes$gene, monkey_genes$gene,mouse_genes$gene)))
all.genes$Human <- as.character(match(all.genes$genes, human_genes$gene))
all.genes$Monkey <- as.character(match(all.genes$genes, monkey_genes$gene))
all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
# all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
all.genes$Human[which(is.na(all.genes$Human))] <- FALSE 
all.genes$Human[which(all.genes$Human != FALSE)] <- TRUE
all.genes$Monkey[which(is.na(all.genes$Monkey))] <- FALSE 
all.genes$Monkey[which(all.genes$Monkey != FALSE)] <- TRUE
all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE
all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
all.genes$Human <- as.logical(all.genes$Human)
all.genes$Monkey <- as.logical(all.genes$Monkey)
all.genes$Mouse <- as.logical(all.genes$Mouse)

library(eulerr)
p3=plot(euler(
  all.genes[ ,2:4]),
  quantities = list(cex = 3),
  labels = NULL,
  main = paste0(subclasses[tmp], " vs. All Inh"),
  fills = c("#ea218e", "#00aeee","#ea7b3d")
)

pdf("2.venn_lineage.pdf",9,3)
as.ggplot(p1)+as.ggplot(p2)+as.ggplot(p3)
dev.off()




#######heatmap - conserved genes#####
my_list <- readRDS("2.barcode_downsample.rds")
aa=Reduce(c,my_list) %>% unique()
newdata_2=newdata_1[,aa]
dim(newdata_2)
sample_recMarkers <- read.xlsx("2.lineage_allmarkers.xlsx")

sample_recMarkers_1=sample_recMarkers %>% group_by(cluster,gene,species) %>% mutate(counts=length(unique(order)))
sample_recMarkers_1=sample_recMarkers_1[sample_recMarkers_1$counts==10,]
table(sample_recMarkers_1$counts)
sample_recMarkers_1=sample_recMarkers_1 %>% distinct(cluster,gene,species,counts,.keep_all = T)

human_cells_markers <- sample_recMarkers_1[sample_recMarkers_1$species=="Human",]
monkey_cells_markers <- sample_recMarkers_1[sample_recMarkers_1$species=="Monkey",]
mouse_cells_markers <- sample_recMarkers_1[sample_recMarkers_1$species=="Mouse",]

human_data <- subset(newdata_2,species=="Human")
table(human_data$new_lineage)
monkey_data <- subset(newdata_2,species=="Monkey")
table(monkey_data$new_lineage)
mouse_data <- subset(newdata_2,species=="Mouse")
table(mouse_data$new_lineage)

subclasses <- c("mgeLin","cgeLin","lgeLin")
genes_to_plot <- NA
for(i in 1:length(subclasses)){
  tmp <- i
  human_genes <- human_cells_markers[grep(subclasses[tmp], human_cells_markers$cluster), ]
  human_genes <- human_genes[!duplicated(human_genes$gene),]
  monkey_genes <- monkey_cells_markers[grep(subclasses[tmp], monkey_cells_markers$cluster), ]
  monkey_genes <- monkey_genes[!duplicated(monkey_genes$gene),]
  mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
  mouse_genes <- mouse_genes[!duplicated(mouse_genes$gene),]
  all.genes <- data.frame(genes = unique(c(human_genes$gene, monkey_genes$gene,mouse_genes$gene)))
  all.genes$Human <- as.character(match(all.genes$genes, human_genes$gene))
  all.genes$Monkey <- as.character(match(all.genes$genes, monkey_genes$gene))
  all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
  all.genes$Human[which(is.na(all.genes$Human))] <- FALSE 
  all.genes$Human[which(all.genes$Human != FALSE)] <- TRUE
  all.genes$Monkey[which(is.na(all.genes$Monkey))] <- FALSE 
  all.genes$Monkey[which(all.genes$Monkey != FALSE)] <- TRUE
  all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE
  all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
  all.genes$Human <- as.logical(all.genes$Human)
  all.genes$Monkey <- as.logical(all.genes$Monkey)
  all.genes$Mouse <- as.logical(all.genes$Mouse)
  all.genes <- all.genes[which(all.genes[,2] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,3] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,4] == TRUE), ]
  gc()
  genes_to_plot <- c(genes_to_plot, as.character(all.genes$genes))
}
genes_to_plot1 <- genes_to_plot[-1]

genes_to_plot2 <- human_cells_markers[human_cells_markers$gene %in% genes_to_plot1,]
genes_to_plot2 <- genes_to_plot2[!duplicated(genes_to_plot2$gene),]
genes_to_plot2=genes_to_plot2 %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
genes_to_plot2$cluster=factor(genes_to_plot2$cluster,levels = c("mgeLin","cgeLin","lgeLin"))
genes_to_plot2=genes_to_plot2 %>% arrange(cluster)

gene_cell_exp <- AverageExpression(human_data,
                                   features = unique(genes_to_plot2$gene),
                                   group.by = 'new_lineage',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene_cell_exp <- gene_cell_exp[c("mgeLin","cgeLin","lgeLin")]
# gene_cell_exp1 <- gene_cell_exp %>% rownames_to_column(var = "gene")
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
marker_exp[marker_exp > 2]=2
marker_exp[marker_exp < -2]=-2
library(ComplexHeatmap)
p1=Heatmap(marker_exp,
           show_column_names = T,
           show_row_names = T,
           cluster_rows = F,
           cluster_columns = F,
           #column_order = new_group_lst1,
           name = "Z score",
           row_names_side =  'left',
           col = colorRampPalette(c("#5e126e","white","#fcaf13"))(100),
           column_title = NULL,
           row_title = NULL,
           border = 'black',
           rect_gp = gpar(col = "black", lwd = 1),
           row_names_gp = gpar(fontsize = 5),
           column_names_gp = gpar(fontsize = 5),
           show_heatmap_legend=T,
           heatmap_width = unit(3, "cm"),
           heatmap_height = unit(10, "cm")
           
)

gene_cell_exp <- AverageExpression(monkey_data,
                                   features = unique(genes_to_plot2$gene),
                                   group.by = 'new_lineage',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene_cell_exp <- gene_cell_exp[c("mgeLin","cgeLin","lgeLin")]
gene_cell_exp1 <- gene_cell_exp %>% rownames_to_column(var = "gene")
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
marker_exp[marker_exp > 2]=2
marker_exp[marker_exp < -2]=-2
library(ComplexHeatmap)
p2=Heatmap(marker_exp,
           show_column_names = T,
           show_row_names = T,
           cluster_rows = F,
           cluster_columns = F,
           name = "Z score",
           row_names_side =  'left',
           col = colorRampPalette(c("#5e126e","white","#fcaf13"))(100),
           column_title = NULL,
           row_title = NULL,
           border = 'black',
           rect_gp = gpar(col = "black", lwd = 1),
           row_names_gp = gpar(fontsize = 5),
           column_names_gp = gpar(fontsize = 5),
           show_heatmap_legend=T,
           heatmap_width = unit(3, "cm"),
           heatmap_height = unit(10, "cm")
           
)

gene_cell_exp <- AverageExpression(mouse_data,
                                   features = unique(genes_to_plot2$gene),
                                   group.by = 'new_lineage',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene_cell_exp <- gene_cell_exp[c("mgeLin","cgeLin","lgeLin")]
gene_cell_exp1 <- gene_cell_exp %>% rownames_to_column(var = "gene")
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
marker_exp[marker_exp > 2]=2
marker_exp[marker_exp < -2]=-2
# library(ComplexHeatmap)
p3=Heatmap(marker_exp,
           show_column_names = T,
           show_row_names = T,
           cluster_rows = F,
           cluster_columns = F,
           #column_order = new_group_lst1,
           name = "Z score",
           row_names_side =  'left',
           col = colorRampPalette(c("#5e126e","white","#fcaf13"))(100),
           column_title = NULL,
           row_title = NULL,
           border = 'black',
           rect_gp = gpar(col = "black", lwd = 1),
           row_names_gp = gpar(fontsize = 5),
           column_names_gp = gpar(fontsize = 5),
           show_heatmap_legend=T,
           heatmap_width = unit(3, "cm"),
           heatmap_height = unit(10, "cm")
           
)




#heatmap - human_specific genes
genes_to_plot <- NA
for(i in 1:length(subclasses)){
  tmp <- i
  human_genes <- human_cells_markers[grep(subclasses[tmp], human_cells_markers$cluster), ]
  monkey_genes <- monkey_cells_markers[grep(subclasses[tmp], monkey_cells_markers$cluster), ]
  mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
  all.genes <- data.frame(genes = unique(c(human_genes$gene, monkey_genes$gene,mouse_genes$gene)))
  all.genes$Human <- as.character(match(all.genes$genes, human_genes$gene))
  all.genes$Monkey <- as.character(match(all.genes$genes, monkey_genes$gene))
  all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
  all.genes$Human[which(is.na(all.genes$Human))] <- FALSE 
  all.genes$Human[which(all.genes$Human != FALSE)] <- TRUE
  all.genes$Monkey[which(is.na(all.genes$Monkey))] <- FALSE 
  all.genes$Monkey[which(all.genes$Monkey != FALSE)] <- TRUE
  all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE
  all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
  all.genes$Human <- as.logical(all.genes$Human)
  all.genes$Monkey <- as.logical(all.genes$Monkey)
  all.genes$Mouse <- as.logical(all.genes$Mouse)
  all.genes <- all.genes[which(all.genes[,2] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,3] == FALSE), ]
  all.genes <- all.genes[which(all.genes[,4] == FALSE), ]
  gc()
  genes_to_plot <- c(genes_to_plot, as.character(all.genes$genes))
}
genes_to_plot2 <- genes_to_plot[-1]


genes_to_plot2 <- human_cells_markers[human_cells_markers$gene %in% genes_to_plot2,]
genes_to_plot2 <- genes_to_plot2[!duplicated(genes_to_plot2$gene),]
genes_to_plot2=genes_to_plot2 %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)
genes_to_plot2$cluster=factor(genes_to_plot2$cluster,levels = c("mgeLin","cgeLin","lgeLin"))
genes_to_plot2=genes_to_plot2 %>% arrange(cluster)

gene_cell_exp <- AverageExpression(human_data,
                                   features = unique(genes_to_plot2$gene),
                                   group.by = 'new_lineage',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene_cell_exp <- gene_cell_exp[c("mgeLin","cgeLin","lgeLin")]
gene_cell_exp1 <- gene_cell_exp %>% rownames_to_column(var = "gene")
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
marker_exp[marker_exp > 2]=2
marker_exp[marker_exp < -2]=-2
library(ComplexHeatmap)
p4=Heatmap(marker_exp,
           show_column_names = T,
           show_row_names = T,
           cluster_rows = F,
           cluster_columns = F,
           #column_order = new_group_lst1,
           name = "Z score",
           row_names_side =  'left',
           col = colorRampPalette(c("#5e126e","white","#fcaf13"))(100),
           column_title = NULL,
           row_title = NULL,
           border = 'black',
           rect_gp = gpar(col = "black", lwd = 1),
           row_names_gp = gpar(fontsize = 5),
           column_names_gp = gpar(fontsize = 5),
           show_heatmap_legend=T,
           heatmap_width = unit(3, "cm"),
           heatmap_height = unit(20, "cm")
           
)




pdf("2.pheatmap_conservedGenes.pdf",8,8)
as.ggplot(p1)+as.ggplot(p2)+as.ggplot(p3)
dev.off()



########human specific gene heatmap######
genes_to_plot <- NA
for(i in 1:length(subclasses)){
  tmp <- i
  human_genes <- human_cells_markers[grep(subclasses[tmp], human_cells_markers$cluster), ]
  monkey_genes <- monkey_cells_markers[grep(subclasses[tmp], monkey_cells_markers$cluster), ]
  mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
  all.genes <- data.frame(genes = unique(c(human_genes$gene, monkey_genes$gene,mouse_genes$gene)))
  all.genes$Human <- as.character(match(all.genes$genes, human_genes$gene))
  all.genes$Monkey <- as.character(match(all.genes$genes, monkey_genes$gene))
  all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
  all.genes$Human[which(is.na(all.genes$Human))] <- FALSE 
  all.genes$Human[which(all.genes$Human != FALSE)] <- TRUE
  all.genes$Monkey[which(is.na(all.genes$Monkey))] <- FALSE 
  all.genes$Monkey[which(all.genes$Monkey != FALSE)] <- TRUE
  all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE
  all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
  all.genes$Human <- as.logical(all.genes$Human)
  all.genes$Monkey <- as.logical(all.genes$Monkey)
  all.genes$Mouse <- as.logical(all.genes$Mouse)
  all.genes <- all.genes[which(all.genes[,2] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,3] == FALSE), ]
  all.genes <- all.genes[which(all.genes[,4] == FALSE), ]
  gc()
  genes_to_plot <- c(genes_to_plot, as.character(all.genes$genes))
}
genes_to_plot2 <- genes_to_plot[-1]


genes_to_plot2 <- human_cells_markers[human_cells_markers$gene %in% genes_to_plot2,]
genes_to_plot2 <- genes_to_plot2[!duplicated(genes_to_plot2$gene),]
genes_to_plot2=genes_to_plot2 %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)
genes_to_plot2$cluster <- factor(genes_to_plot2$cluster,levels=c("mgeLin","cgeLin","lgeLin"))
genes_to_plot2 <- genes_to_plot2 %>% arrange(cluster)





gene_cell_exp <- AverageExpression(newdata_2,
                                   features = unique(genes_to_plot2$gene),
                                   group.by = 'species_lineage',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene_cell_exp <- gene_cell_exp[c("Human_mgeLin","Monkey_mgeLin","Mouse_mgeLin",
                                 "Human_cgeLin","Monkey_cgeLin","Mouse_cgeLin",
                                 "Human_lgeLin","Monkey_lgeLin","Mouse_lgeLin")]
gene_cell_exp$max=colnames(gene_cell_exp)[max.col(gene_cell_exp, ties.method = "first")]
gene_cell_exp1 <- gene_cell_exp %>% rownames_to_column(var = "gene")

gene_cell_exp1 <- gene_cell_exp1[
  (gene_cell_exp1$max == "Human_mgeLin" & gene_cell_exp1$gene %in% genes_to_plot2$gene[genes_to_plot2$cluster == "mgeLin"]) |
    (gene_cell_exp1$max == "Human_cgeLin" & gene_cell_exp1$gene %in% genes_to_plot2$gene[genes_to_plot2$cluster == "cgeLin"]) |
    (gene_cell_exp1$max == "Human_lgeLin" & gene_cell_exp1$gene %in% genes_to_plot2$gene[genes_to_plot2$cluster == "lgeLin"]),
]
gene_cell_exp1=gene_cell_exp1 %>% select(-c("max"))
rownames(gene_cell_exp1) <- NULL
gene_cell_exp1=column_to_rownames(gene_cell_exp1,var = "gene")
marker_exp <- t(scale(t(gene_cell_exp1),scale = T,center = T))
marker_exp[marker_exp > 2]=2
marker_exp[marker_exp < -2]=-2
library(ComplexHeatmap)
p7=Heatmap(marker_exp,
           show_column_names = T,
           show_row_names = T,
           cluster_rows = F,
           cluster_columns = F,
           #column_order = new_group_lst1,
           name = "Z score",
           row_names_side =  'left',
           col = colorRampPalette(c("#5e126e","white","#fcaf13"))(100),
           column_title = NULL,
           row_title = NULL,
           border = 'black',
           rect_gp = gpar(col = "black", lwd = 1),
           row_names_gp = gpar(fontsize = 5),
           column_names_gp = gpar(fontsize = 5),
           show_heatmap_legend=T,
           heatmap_width = unit(10, "cm"),
           heatmap_height = unit(21.5, "cm")
           
           
)


####primate specific gene
genes_to_plot <- NA
for(i in 1:length(subclasses)){
  tmp <- i
  human_genes <- human_cells_markers[grep(subclasses[tmp], human_cells_markers$cluster), ]
  monkey_genes <- monkey_cells_markers[grep(subclasses[tmp], monkey_cells_markers$cluster), ]
  mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
  all.genes <- data.frame(genes = unique(c(human_genes$gene, monkey_genes$gene,mouse_genes$gene)))
  all.genes$Human <- as.character(match(all.genes$genes, human_genes$gene))
  all.genes$Monkey <- as.character(match(all.genes$genes, monkey_genes$gene))
  all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
  all.genes$Human[which(is.na(all.genes$Human))] <- FALSE 
  all.genes$Human[which(all.genes$Human != FALSE)] <- TRUE
  all.genes$Monkey[which(is.na(all.genes$Monkey))] <- FALSE 
  all.genes$Monkey[which(all.genes$Monkey != FALSE)] <- TRUE
  all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE
  all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
  all.genes$Human <- as.logical(all.genes$Human)
  all.genes$Monkey <- as.logical(all.genes$Monkey)
  all.genes$Mouse <- as.logical(all.genes$Mouse)
  all.genes <- all.genes[which(all.genes[,2] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,3] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,4] == FALSE), ]
  gc()
  genes_to_plot <- c(genes_to_plot, as.character(all.genes$genes))
}
genes_to_plot2 <- genes_to_plot[-1]


genes_to_plot2 <- human_cells_markers[human_cells_markers$gene %in% genes_to_plot2,]
genes_to_plot2 <- genes_to_plot2[!duplicated(genes_to_plot2$gene),]
genes_to_plot2=genes_to_plot2 %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)
genes_to_plot2$cluster <- factor(genes_to_plot2$cluster,levels=c("mgeLin","cgeLin","lgeLin"))
genes_to_plot2 <- genes_to_plot2 %>% arrange(cluster)


gene_cell_exp <- AverageExpression(newdata_2,
                                   features = unique(genes_to_plot2$gene),
                                   group.by = 'species_lineage',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene_cell_exp <- gene_cell_exp[c("Human_mgeLin","Monkey_mgeLin","Mouse_mgeLin",
                                 "Human_cgeLin","Monkey_cgeLin","Mouse_cgeLin",
                                 "Human_lgeLin","Monkey_lgeLin","Mouse_lgeLin")]
gene_cell_exp$max=colnames(gene_cell_exp)[max.col(gene_cell_exp, ties.method = "first")]
gene_cell_exp1 <- gene_cell_exp %>% rownames_to_column(var = "gene")

gene_cell_exp1 <- gene_cell_exp1[!(gene_cell_exp1$max %in% c("Mouse_mgeLin","Mouse_cgeLin","Mouse_lgeLin")),
]
gene_cell_exp1=gene_cell_exp1 %>% select(-c("max"))
rownames(gene_cell_exp1) <- NULL
gene_cell_exp1=column_to_rownames(gene_cell_exp1,var = "gene")
marker_exp <- t(scale(t(gene_cell_exp1),scale = T,center = T))
marker_exp[marker_exp > 2]=2
marker_exp[marker_exp < -2]=-2
library(ComplexHeatmap)
p8=Heatmap(marker_exp,
           show_column_names = T,
           show_row_names = T,
           cluster_rows = F,
           cluster_columns = F,
           #column_order = new_group_lst1,
           name = "Z score",
           row_names_side =  'left',
           col = colorRampPalette(c("#5e126e","white","#fcaf13"))(100),
           column_title = NULL,
           row_title = NULL,
           border = 'black',
           rect_gp = gpar(col = "black", lwd = 1),
           row_names_gp = gpar(fontsize = 5),
           column_names_gp = gpar(fontsize = 5),
           show_heatmap_legend=T,
           heatmap_width = unit(10, "cm"),
           heatmap_height = unit(26.5, "cm")
           
           
)


pdf("2.pheatmap_human_primate_specific.pdf",15,15)
as.ggplot(p7)+as.ggplot(p8)
dev.off()






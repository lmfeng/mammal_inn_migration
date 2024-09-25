##fig5a
library(SCP)
library(ComplexHeatmap)
library(patchwork)
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")

gexpr1 <- subset(gexpr,lineage=="lgeLin")
gexpr1 = subset(gexpr1,subset = lobe!= "NCX")
gexpr1$lobe <- factor(gexpr1$lobe,levels=c("FC","MSC","OC","TC","Insula","AMY","HIP","STR","GE"))
saveRDS(gexpr1,"lgeLin.Rds")
gexpr1$subtype <- factor(gexpr1$subtype,levels=c("InN MEIS2 PAX6 SCGN","InN MEIS2 FOXP2 TSHZ1","InN MEIS2 ZIC1 ZIC2","InN MEIS2 ISL1 PAX6",
                                                 "InN MEIS2 SIX3 ZFHX3","InN MEIS2 ISL1 ZFHX3","InN MEIS2 ISL1 FOXP1","InN MEIS2 FOXP1 PENK"))
Idents(gexpr1) = "lobe"

p1 <- CellStatPlot(gexpr1, stat.by = c("subtype", "lobe"), plot_type = "sankey")


pdf("Fig5A.pdf", width = 18, height =10)
par(omi=c(0.1,0.1,0.1,0.1))
p1+wrap_elements()
dev.off()

###########################################################################################
##fig5b
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
dim(gexpr)



gexpr$subtype[which(gexpr$subtype %in% c("InN MEIS2 PAX6 SCGN","InN MEIS2 FOXP2 TSHZ1"))] <- "InN MEIS2 SP8 FOXP2"
gexpr.orig <- gexpr
gexpr=NULL

gexpr <- subset(gexpr.orig,subset=subtype=="InN MEIS2 SP8 FOXP2" & transAge!="GW2trimester")
gexpr <- subset(gexpr,subset=transAge!="GW3trimester")
gexpr$transAge <- factor(gexpr$transAge,levels=c(sort(unique(gexpr$transAge))))
gexpr1<-subset(gexpr,subset=lobe %in% c("FC","MSC","GE"))
table(gexpr1$labsite,gexpr1$lobe)


gexpr1$transAge=as.character(gexpr1$transAge)
gexpr1$transAge=as.factor(gexpr1$transAge)
# gexpr1$new_group=paste0(gexpr1$lobe,"_",gexpr1$species)
library(ggpubr)
library(ggbreak)
FC_data=subset(gexpr1,lobe=="FC")
FC_stat = as.data.frame(table(FC_data$species,FC_data$transAge))
colnames(FC_stat) <- c("species","transAge","num")
FC_stat$species <- factor(FC_stat$species,levels=c("Human","Monkey"))
FC_stat = FC_stat[order(FC_stat$species),]
FC_stat$num[FC_stat$num<10]=0
#FC_stat$num[which(FC_stat$num!=0)]=log2(FC_stat$num[which(FC_stat$num!=0)])
p1=ggplot()+geom_bar(data=FC_stat,
                     aes(x=transAge,y=num,fill=species),
                     stat="identity",
                     colour = "black")+
  scale_y_continuous(breaks = c(0,500,1000,3000,5000))+
  scale_fill_manual(breaks  = c("Human","Monkey"),values=c( "#ed0000","#00468b"))+
  theme_bw()+theme(text = element_text(size = 10),panel.grid = element_blank())


# p1=ggbarplot(FC_stat, x="transAge", y="num", fill = "species", color = "white", label =F,lab.pos="out",palette =  "lancet",ylab = "FC")+ 
# 
#     theme(axis.title.x=element_blank(),legend.position = "none",
#         axis.text = element_text(size = 10))+
#   scale_y_break(c(500,700),space = 0.1,scales = "free")+  
#   scale_fill_manual(breaks  = c("Human","Monkey"),values=c( "#ed0000","#00468b"))



MSC_data=subset(gexpr1,lobe=="MSC")
MSC_stat = as.data.frame(table(MSC_data$species,MSC_data$transAge))
colnames(MSC_stat) <- c("species","transAge","num")
MSC_stat$species <- factor(MSC_stat$species,levels=c("Human","Monkey"))
MSC_stat = MSC_stat[order(MSC_stat$species),]
MSC_stat$num[MSC_stat$num<10]=0
#FC_stat$num[which(FC_stat$num!=0)]=log2(FC_stat$num[which(FC_stat$num!=0)])
p2=ggplot()+geom_bar(data=MSC_stat,
                     aes(x=transAge,y=num,fill=species),
                     stat="identity",
                     colour = "black")+
  scale_y_continuous(breaks = c(0,500,1000,2000,3000))+
  scale_fill_manual(breaks  = c("Human","Monkey"),values=c( "#ed0000","#00468b"))+
  theme_bw()+theme(text = element_text(size = 10),panel.grid = element_blank())


#MSC_stat$num[which(MSC_stat$num!=0)]=log2(MSC_stat$num[which(MSC_stat$num!=0)])
# p2=ggbarplot(MSC_stat, x="transAge", y="num", fill = "species", color = "white", label =F,lab.pos="out",palette =  "lancet",ylab = "MSC")+ 
#   theme(axis.title.x=element_blank(),legend.position = "none",
#         axis.text = element_text(size = 10))+
#   scale_y_break(c(200,800),space = 0.1,scales = "free")+  
#   scale_fill_manual(breaks  = c("Human","Monkey"),values=c( "#ed0000","#00468b"))


GE_data=subset(gexpr1,lobe=="GE")
GE_stat = as.data.frame(table(GE_data$species,GE_data$transAge))
colnames(GE_stat) <- c("species","transAge","num")
GE_stat$species <- factor(GE_stat$species,levels=c("Human","Monkey"))
GE_stat = GE_stat[order(GE_stat$species),]
#GE_stat$num[which(GE_stat$num!=0)]=log2(GE_stat$num[which(GE_stat$num!=0)])
GE_stat$num[GE_stat$num<10]=0
#FC_stat$num[which(FC_stat$num!=0)]=log2(FC_stat$num[which(FC_stat$num!=0)])
p3=ggplot()+geom_bar(data=GE_stat,
                     aes(x=transAge,y=num,fill=species),
                     stat="identity",
                     colour = "black")+
  scale_y_continuous(breaks = c(0,500,1000,2000,5000,8000))+
  scale_fill_manual(breaks  = c("Human","Monkey"),values=c( "#ed0000","#00468b"))+
  theme_bw()+theme(text = element_text(size = 10),panel.grid = element_blank())

# p4=ggbarplot(GE_stat, x="transAge", y="num", fill = "species", color = "white", label =F,lab.pos="out",palette =  "lancet",ylab = "GE")+ 
#   theme(axis.title.x=element_blank(),legend.position = "none",
#         axis.text = element_text(size = 10))+
#   scale_y_break(c(600,1200),space = 0.1,scales = "free",ticklabels = c(2000,4000,6000,8000))+  #space截断处的宽度，scales截断上下的比例
#   scale_fill_manual(breaks  = c("Human","Monkey"),values=c( "#ed0000","#00468b"))



pdf("Fig5B.pdf",8,6,onefile = FALSE)
p1+p2+p3+patchwork::plot_layout(ncol = 1,guides = "collect")
dev.off()





###########################################################################################
##fig5c
library(monocle)
library(Seurat)
library(ggplot2)
library(Matrix)
library(SeuratObject)
library(patchwork)
library(dplyr)
library(viridis)
library(openxlsx)

age.use="E93"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  age.use <- args[1]
}

####
set.seed(100)


####loading data
gexpr <- NULL
load('mammal_interneuron_migration.final.p6.1.RData')
gexpr.orig = gexpr




###@@@@@@@@@
gexpr <- NULL
gexpr = subset(gexpr.orig,  subtype %in% c("InN MEIS2 PAX6 SCGN") & age %in% age.use & 
                 # region %in% c("LGE", "DFC", "MFC", "OFC", "M1C", "S1C", "A1C","IPC"))
                 region %in% c("LGE", "DFC", "MFC", "OFC", "M1C", "S1C", "A1C","IPC","VFC","PCC"))



##downsample
Idents(gexpr) <- "lobe"
gexpr6 <- subset(x = gexpr, downsample = 1000)
#gexpr6 = gexpr

#
gexprMat = GetAssayData(gexpr6, slot='counts')
meta = gexpr6@meta.data

#####compute hvg
gexpr6 <- FindVariableFeatures(gexpr6, nfeatures = 500)
hvg <- VariableFeatures(gexpr6)



#####@@@@@@@@@@@@@@@@------run  monocle---			
####pheno type
pd = new('AnnotatedDataFrame', data = meta)

###feature data
geneSymbol = rownames(gexprMat)
geneAnnot = data.frame(gene_short_name = geneSymbol)
rownames(geneAnnot) = geneSymbol
fd = new('AnnotatedDataFrame', data = geneAnnot)

#####-----expression matrix
cds = newCellDataSet(gexprMat, phenoData = pd, featureData = fd,
                     lowerDetectionLimit = 0.1, expressionFamily = negbinomial.size())

####----Estimate size factors and dispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)



cds_subset <- setOrderingFilter(cds, ordering_genes = hvg)
cds_subset <- reduceDimension(cds_subset, max_components =15, reduction_method = 'DDRTree',	norm_method = 'vstExprs'
                              #,residualModelFormulaStr='~orig.ident'
)
suppressWarnings(cds_subset <- orderCells(cds_subset)) 
saveRDS(cds_subset,file = paste0(age.use,"monocle.Rds"))




cds_subset <- orderCells(cds_subset,reverse = T)
##based on lobe
pdf(paste0(age.use,"_lobe_Pseudotime.pdf"),6,6)
par(omi = c(0.1, 0.1, 0.1, 0.1))

##
plot_cell_trajectory(cds_subset, color_by = "Pseudotime",theta=0,cell_size = 1,show_branch_points = F)+ theme(legend.position = "right")+
  theme(aspect.ratio = 1/1)
cols <- c("#69b1d7","#f3b74d","#e38e8b") %>% setNames(c("GE","FC","MSC"))
p1 <- plot_cell_trajectory(cds_subset, color_by = 'lobe', theta=0,cell_size = 1,show_branch_points = F)+
  scale_color_manual(breaks = c("GE", "FC", "MSC"), values=cols) + theme(legend.position = "right")+
  theme(aspect.ratio = 1/1)
print(p1)
dev.off()


pdf(paste0(age.use,"_facetState.pdf"),20,20)
plot_cell_trajectory(cds_subset, color_by = "State",show_branch_points = F)+facet_wrap(~State)
plot_cell_trajectory(cds_subset, color_by = "State",show_branch_points = T)
dev.off()



###############################################pie
mydata=pData(cds_subset)
mydata$lineage_123=case_when(mydata$State %in% c(18:20) ~ "FClineage",
                             mydata$State %in% c(1:4,21) ~ "GElineage",
                             mydata$State %in% c(5,7,8,10,11,12,14:17) ~ "MSClineage",
)
table(mydata$lineage_123,mydata$State)

mydata <- mydata %>% group_by(lineage_123,lobe) %>% summarise(counts=n())
mydata <- mydata %>% group_by(lineage_123) %>% mutate(sums=sum(counts))
mydata$Freq <- round(mydata$counts/mydata$sums,4)
p_list=list()
for (i in 1:3) {
  lineage_123.use=unique(mydata$lineage_123)[i]
  mydata1=mydata[mydata$lineage_123==lineage_123.use,]
  p1=ggplot(mydata1,aes(x="",y=Freq,fill=lobe))+geom_bar(stat = "identity")+labs(title = lineage_123.use,x = '', y = '')+
    scale_fill_manual(values =cols[mydata1$lobe])+coord_polar(theta = "y")+
    theme(axis.text.x = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5,face = "bold"),panel.background = element_blank()
    )
  p_list[[i]]=p1
}

pdf(paste0(age.use,"_pie.pdf"),9,3)
wrap_plots(p_list,ncol = 3,guides = "collect")
dev.off()


##############################BEAM 
branch_point.use=2
BEAM_res <- BEAM(cds_subset[hvg,rownames(mydata1)], branch_point = branch_point.use, cores = 8)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
table(BEAM_res$qval < 1e-4)



tmp1=plot_genes_branched_heatmap(cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),],
                                 branch_point = branch_point.use,
                                 num_clusters = 4,
                                 cores = 8,
                                 use_gene_short_name = TRUE,
                                 show_rownames = TRUE,
                                 return_heatmap = T)
library(ggplotify)
pdf(paste0(age.use,"_BEAM_filter_genes.pdf"),16,8)
as.ggplot(tmp$ph_res)+as.ggplot(tmp1$ph_res)
dev.off()

#########genes plot
pdf(paste0(age.use,"_plotgenes_point.pdf"),12,12)
par(omi=c(0.1,0.1,0.1,0.1))
plot_genes_branched_pseudotime(cds = cds_subset[c("CCND2","ETV1","CCK","PCP4","NR2F1","CADM1","MEIS2","ERBB4","MDK","NR2F2","PROX1"),rownames(mydata1)],
                               branch_point = 2,
                               color_by = "lobe",
                               cell_size=0.3,
                               ncol = 3)+scale_color_manual(breaks = c("GE", "FC", "MSC"), values=cols)+theme(aspect.ratio = 1/1)
dev.off()


###########################################################################################
##fig5f
library(ggpubr)
library(patchwork)
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
dim(gexpr)


gexpr1 <- gexpr
gexpr=NULL
gexpr <- subset(gexpr1,subset=lobe=="FC")
table(gexpr$labsite,gexpr$region)

gexpr <- subset(gexpr,labsite!="Linnarsson2022")
table(gexpr$labsite,gexpr$region)




gexpr$new_region = gexpr$region
gexpr$new_region[which(gexpr$species=="Human")] ="DFC"
gexpr$new_region[which(gexpr$species=="Monkey" & gexpr$new_region %in% c("FR","PFC"))] ="DFC"

gexpr$new_group = gexpr$labsite
gexpr$new_group[which(gexpr$orig.ident=="E93")] ="E93"
gexpr$new_group[which(gexpr$orig.ident=="RMB683")] ="RMB683"
gexpr$new_group[which(gexpr$orig.ident=="RMB691")] ="RMB691"

Fig5E_stat = as.data.frame(table(gexpr$new_group,gexpr$new_region))
gexpr_SCGN=subset(gexpr,subtype=="InN MEIS2 PAX6 SCGN")
Fig5E_stat_SCGN = as.data.frame(table(gexpr_SCGN$new_group,gexpr_SCGN$new_region))
Fig5E_stat=merge(Fig5E_stat,Fig5E_stat_SCGN,by = c("Var1","Var2"),all.x = T)
colnames(Fig5E_stat) <- c("Sample","Region","num","num_SCGN")
Fig5E_stat$Sample <- factor(Fig5E_stat$Sample,levels=c("KriegsteinLab","Kriegstein2022","RakicLab","PollenLab","E93","RMB683","RMB691"))
Fig5E_stat = Fig5E_stat[order(Fig5E_stat$Sample),]
Fig5E_stat$num_SCGN[Fig5E_stat$num_SCGN==0]=NA
Fig5E_stat$num[Fig5E_stat$num==0]=NA
Fig5E_stat$freq=Fig5E_stat$num_SCGN/Fig5E_stat$num

# Fig5E_stat$num[which(Fig5E_stat$num!=0)]=log2(Fig5E_stat$num[which(Fig5E_stat$num!=0)])
# Fig5E_stat1$num[which(Fig5E_stat1$num!=0)]=log2(Fig5E_stat1$num[which(Fig5E_stat1$num!=0)])


p1=ggbarplot(Fig5E_stat, x="Region", y="num", fill = "Sample", size=1,color = "white", label =T,palette = "npg",ylab ="Cell number",
             x.text.angle=0,position = position_dodge(),title = "all subtypes in FC")+scale_y_continuous(limits = c(0,10428),breaks = c(0,2500,5000,7500,10000))

p2=ggbarplot(Fig5E_stat, x="Region", y="num_SCGN", fill = "Sample", size=1,color = "white", label =T,palette = "npg",ylab ="Cell number",
             x.text.angle=0,position = position_dodge(),title = "InN MEIS2 PAX6 SCGN in FC")+scale_y_continuous(limits = c(0,2490),breaks = c(0,500,1000,1500,2000,2500))





pdf("Fig5F.pdf", width = 6, height =10)
p1+p2+plot_layout(nrow = 2,guides = "collect")&theme(legend.position = "bottom")
dev.off()



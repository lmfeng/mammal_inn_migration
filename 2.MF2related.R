#####fig2a
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
gexpr$lobe[gexpr$lobe %in% c("FC","MSC","TC","OC","Insula","NCX")]="CTX"
gexpr$group="other"
gexpr$group[gexpr$lineage=="inIPC"&gexpr$lobe=="GE"&gexpr$species=="Human"]="Human"
gexpr$group[gexpr$lineage=="inIPC"&gexpr$lobe=="GE"&gexpr$species=="Monkey"]="Macaque"


gexpr$group1="other"
gexpr$group1[gexpr$lineage=="inIPC"&gexpr$lobe=="CTX"&gexpr$species=="Human"]="Human"
gexpr$group1[gexpr$lineage=="inIPC"&gexpr$lobe=="CTX"&gexpr$species=="Monkey"]="Macaque"

library(ggrastr)
scviumap_tx = gexpr@reductions$scVI@cell.embeddings %>% 
  as.data.frame() %>% cbind(group = gexpr$group,group1=gexpr$group1)

set.seed(1234)
df=scviumap_tx[sample(rownames(scviumap_tx)[scviumap_tx$group!="other"]),]
set.seed(1234)
df1=scviumap_tx[sample(rownames(scviumap_tx)[scviumap_tx$group1!="other"]),]

p3=ggplot(scviumap_tx, aes(x=scVI_1, y=scVI_2, color=group)) + geom_point_rast(size=0.1) + 
  scale_color_manual(values=c("Human" = alpha("#ea218e",0.5), 
                              "Macaque" = alpha("#00aeee",0.5),
                              "other"="lightgrey"),
                     labels=c("Human inIPC","Macaque inIPC","Non-inIPC"))+
  geom_point_rast(data=df, mapping=aes(scVI_1, scVI_2,color=group), size=0.1)+theme_bw()+
  theme(aspect.ratio=1/1,panel.grid = element_blank(),plot.title = element_text(hjust = 0.5),axis.title = element_blank(),
        axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),panel.border = element_blank())+
  labs(title = "GE")+
  guides(color = guide_legend(override.aes = list(size = c(2, 2, 2))))
p4=ggplot(scviumap_tx, aes(x=scVI_1, y=scVI_2, color=group1)) + geom_point_rast(size=0.1) + 
  scale_color_manual(values=c("Human" = alpha("#ea218e",0.5), 
                              "Macaque" = alpha("#00aeee",0.5),
                              "other"="lightgrey"),
                     labels=c("Human inIPC","Macaque inIPC","Non-inIPC"))+
  geom_point_rast(data=df1, mapping=aes(scVI_1, scVI_2,color=group1), size=0.1)+theme_bw()+
  theme(aspect.ratio=1/1,panel.grid = element_blank(),plot.title = element_text(hjust = 0.5),axis.title = element_blank(),
        axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),panel.border = element_blank())+
  labs(title = "CTX")+
  guides(color = guide_legend(override.aes = list(size = c(2, 2, 2))))



pdf("Fig2A_UMAP.pdf",10,10)
library(patchwork)
p1+p2+p3+p4+plot_layout(ncol=2)
dev.off()




####fig2b
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
gexpr1 <- subset(gexpr,subset = lineage=="inIPC" &species=="Human" & lobe %in% c("FC","Insula","MSC","OC","NCX","TC","GE"))
gexpr1@meta.data$new_lobe <- "-"
gexpr1@meta.data$new_lobe[which(gexpr1@meta.data$lobe == "GE")] <- "GE"
gexpr1@meta.data$new_lobe[which(gexpr1@meta.data$lobe != "GE")] <- "CTX"
human_stat <- as.data.frame(table(gexpr1$new_lobe,gexpr1$transAge))

gexpr2 <- subset(gexpr,subset = lineage=="inIPC" &species=="Monkey" & lobe %in% c("FC","Insula","MSC","OC","NCX","TC","GE"))
gexpr2@meta.data$new_lobe <- "-"
gexpr2@meta.data$new_lobe[which(gexpr2@meta.data$lobe == "GE")] <- "GE"
gexpr2@meta.data$new_lobe[which(gexpr2@meta.data$lobe != "GE")] <- "CTX"
monkey_stat <- as.data.frame(table(gexpr2$new_lobe,gexpr2$transAge))

colnames(human_stat) <- c("region","age","num")
human_stat$age <- as.character(human_stat$age)
#ge_stat <- aggregate(num ~ celltype+ age, data = ge_stat, sum)
human_stat$species <- "Human"
colnames(monkey_stat) <- c("region","age","num")
monkey_stat$age <- as.character(monkey_stat$age)
#ctx_stat <- aggregate(num ~ celltype+ age, data = ctx_stat, sum)
monkey_stat$species <- "Monkey"


result <- rbind(human_stat,monkey_stat)
result$num  <- as.numeric(result$num)
result  <-result[which(result$num!=0),]
result$num  <-log(result$num+1,10)
result$num[which(result$species == "Monkey")] <- result$num[which(result$species == "Monkey")]*-1
result  <-result[which(result$num!=0),]

library(ggplot2)
library(forcats)

pdf("Fig2B_cellnum_barplot.pdf",8,4)
par(omi=c(0.1,0.1,0.1,0.1))
p1 <- ggplot()+geom_bar(data=result,
                        aes(x=num,y=age,fill=region),
                        stat="identity",
                        colour = "grey50")+
  scale_fill_manual(values = c("#e4d1a3", "#439595"))+
  #scale_x_continuous(breaks = seq(-10,10,5),labels = as.character(abs(seq(-10,10,5))))+
  theme_bw()+
  labs(x="Number of cells (log1p)",y="transAge")+
  coord_flip()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#X轴
  geom_segment(aes(y = 1, yend = 1,x = -5, xend = -8),arrow = arrow(length = unit(0.2, "cm")),colour="grey36",linewidth=1)+
  geom_segment(aes(y = 1, yend = 1,x = 5, xend = 8),arrow = arrow(length = unit(0.2, "cm")),colour="grey36",linewidth=1)+
  annotate("text", x = 9 , y = 1,label = "Human",colour="grey36",size=4)+ 
  annotate("text", x = -9 , y = 1,label = "Monkey",colour="grey36",size=4)+
  theme(text = element_text(size = 10),panel.grid = element_blank()) 
print(p1)
dev.off()



###fig2d
library(openxlsx)
mydata <- read.xlsx("fig2D.organized_data.xlsx")

mydata1=mydata[mydata$freq >= 5,]  
mydata1=mydata1[mydata1$avg_log2FC>log2(1.5) &mydata1$p_val_adj<0.05,]
mydata1 %>% group_by(species,cluster) %>% summarise(count=length(unique(gene)))


mydata1$avg_log2FC <- round(mydata1$avg_log2FC,2)
mydata1$log10_p_val_adj <- -round(log10(mydata1$p_val_adj1),digits = 2)  
mydata1$avg_log2FC[mydata1$cluster=="GE"] <- -mydata1$avg_log2FC[mydata1$cluster=="GE"]
mydata1$log10_p_val_adj[mydata1$log10_p_val_adj>=30] <- 30


#Draw background column
mydata1 %>% group_by(species) %>% summarise(max=max(avg_log2FC),min=min(avg_log2FC))
dfbar<-data.frame(x=c("Human","Monkey"),
                  y=c(4,2.5))
dfbar$x <- factor(dfbar$x,levels = c("Human","Monkey"))

dfbar1<-data.frame(x=c("Human","Monkey"),
                   y=c(-2,-2))
dfbar1$x <- factor(dfbar1$x,levels = c("Human","Monkey"))

## set label
aa= mydata1 %>% group_by(gene) %>% mutate(reversed=length(unique(cluster)))
aa=aa %>% group_by(cluster,gene) %>% mutate(shared=n())
##select shared gene
bb=aa[aa$shared>1&aa$reversed==1,]
bb %>% group_by(species,cluster) %>% summarise(count=length(unique(gene)))
mm=bb %>% group_by(species,cluster) %>% top_n(3,wt = abs(avg_log2FC))
bb=bb[bb$gene %in% unique(mm$gene),]


##select top gene
cc=aa[aa$shared==1&aa$reversed==1,]
cc=cc %>% group_by(species,cluster) %>% top_n(8,wt = abs(avg_log2FC))
cc %>% group_by(species,cluster) %>% summarise(count=length(unique(gene)))

top10sig0 <- rbind(bb,cc) %>% unique()
##remove TMPO, TMEM123,ANKS1B,CADPS,ANP32E
top10sig0=top10sig0[!(top10sig0$gene %in% c("TMPO","TMEM123","ANKS1B","CADPS","ANP32E")),]
write.xlsx(top10sig0,file = "fig2D.label_genes.xlsx")


#Add the cluster color block label for the X axis
dfcol<-data.frame(x=c("Human","Monkey"),y=0,labels=c("Human","Monkey"))
dfcol$x <- factor(dfcol$x,levels = c("Human","Monkey"))
mycol <- c("#ea218e","#00aeee")
# mycol <- c("#ff4545","#4DBBD5FF","#f49a9b","#87e0ff")
# mycol <- c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F","#F39B7F7F")


library(ggplot2)
library(ggrepel)
set.seed(1234)
p2=ggplot()+
  #Draw background column
  geom_col(data = dfbar,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)+
  #volcano plot
  geom_jitter(data = mydata1,
              aes(x = species, y =avg_log2FC,color=log10_p_val_adj,size=log10_p_val_adj),width =0.4)+
  scale_color_gradientn(values = seq(0,1,0.25),
                        colors = c("#4575b4","#91bfdb","#e0f3f8","#fc8d59","#d73027"),breaks=c(10,20,30),labels=c(10,20,">=30"))+
  scale_size_continuous(range = c(0.1,3),breaks=c(10,20,30),labels=c(10,20,">= 30"))+
  scale_y_continuous(limits = c(-3.5,4),breaks = c(-3,-2,-1,0,1,2,3,4))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.line.y = element_line(color = "black",
                                   linewidth = 1.2),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title =element_blank(),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6))+
  theme(aspect.ratio = 1/1)+
  xlab("Species")+
  ylab("Log2FC)")+
  geom_tile(data = dfcol,aes(x=x,y=y),height=0.4,color = "black",alpha = 0.6,fill=mycol,show.legend = F)+
  geom_text(data = dfcol,aes(x=x,y=y,label=labels),size =4,color ="white",fontface="bold")+
  geom_text_repel(data=top10sig0, aes(x=species,y=avg_log2FC,label=gene),color="black",size = 1.5,fontface="bold",max.overlaps = 50,box.padding=unit(0.3, "lines"), point.padding=unit(1.6, "lines"),seed=1234,min.segment.length = Inf)+
  annotate('text',x=1,y=3.5,label='CTX',size=5,color='black')+
  annotate('text',x=1,y=-3.5,label='GE',size=5,color='black')

pdf("fig2D.volcano_M&H.pdf",4,6)
par(mar=c(0.1,0.1,0.1,0.1))
p2+theme(legend.position = "bottom")
dev.off()



##fig2e
library(clusterProfiler)
library(org.Hs.eg.db)
mydata <- read.xlsx("fig2D.organized_data.xlsx")
ids=bitr(mydata1$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') 
setdiff(unique(mydata1$gene),unique(ids$SYMBOL))
markers=merge(mydata1,ids,by.x='gene',by.y='SYMBOL')

mydata1$avg_log2FC[mydata1$cluster=="GE"] <- -mydata1$avg_log2FC[mydata1$cluster=="GE"]
gene_up <- markers$ENTREZID[markers$cluster=="CTX"]
gene_down <- markers$ENTREZID[markers$cluster=="GE"]
gene_diff <- c(gene_up,gene_down)




##CTX
go.up <- enrichGO(gene = gene_up,
                  OrgDb = org.Hs.eg.db,
                  ont = "ALL" ,
                  keyType = "ENTREZID",  
                  pAdjustMethod = "fdr",  
                  readabl = TRUE       
)



#GE
go.down <- enrichGO(gene = gene_down,
                    OrgDb = org.Hs.eg.db,
                    ont = "ALL" ,
                    keyType = "ENTREZID",  
                    pAdjustMethod = "fdr",
                    readabl = TRUE       
)



up.data<-go.up %>% as.data.frame() %>% subset(p.adjust<0.000001)
up.data$group=-1
up.data$region="CTX"
down.data<-go.down %>% as.data.frame() %>% subset(p.adjust<0.000001)
down.data$group=1
down.data$region="GE"
all_data <- rbind(up.data,down.data)


library(ggthemes)
go.GO_plot <- function(up.data,down.data){
  dat=rbind(up.data,down.data)
  colnames(dat)
  dat$p.adjust = -log10(dat$p.adjust)
  dat$p.adjust=dat$p.adjust*dat$group 

  aa=dat[dat$Description %in% c("mitotic cell cycle phase transition"),]
  dat=dat[order(dat$p.adjust,decreasing = T),]  
  dat=dat %>% group_by(group,ONTOLOGY) %>% top_n(6,wt = abs(p.adjust))
  
  dat=rbind(dat,aa)
  
  gk_plot <- ggplot(dat,aes(reorder(Description, p.adjust), y=p.adjust)) +
    geom_bar(aes(fill=factor(group)),stat="identity", width=0.7, position=position_dodge(0.7)) +
    coord_flip() +
    scale_fill_manual(values=c("#e4d1a3", "#439595"),breaks = c(-1,1),labels=c("CTX","GE"),name="Region") +
    scale_y_continuous(breaks = c(-10,0,10,20,30,40,50,60))+
    labs(x="Description", y="P.adjust" ) +
    theme_pander()  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks.length.x = unit(-0.10, "cm"),
          axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),
          axis.ticks.x = element_line(colour = "black",size = 0.3) ,                      
          axis.ticks.y = element_blank(),
          axis.text.y  = element_text(hjust=0),
          panel.background = element_rect(fill=NULL, colour = 'white')
    )
}

p3 <- go.GO_plot(up.data ,down.data)

pdf("fig2E.GOenrichment.pdf",6,6)
par(mar=c(0.1,0.1,0.1,2))
p3
dev.off()


##fig2f
##complexHeatmap
library(dplyr)
library(openxlsx)
gexpr <- readRDS("mammal_InN_migration.final.p8_01.rds")
gexpr.m_h <- subset(gexpr,lineage =="inIPC" & lobe %in% c("FC","OC","MSC","TC"))
# gexpr.m_h <- subset(gexpr.m_h,subtype != "inIPC ASCL1 GPX1")
gexpr.m_h=subset(gexpr.m_h,labsite!="Kriegstein2022")

table(gexpr.m_h$labsite,gexpr.m_h$lobe)




gexpr.m_h$new_group_2 <- paste(gexpr.m_h$lobe,gexpr.m_h$species,gexpr.m_h$subtype,sep="_")
aa=as.data.frame(table(gexpr.m_h$new_group_2)>10)
aa=aa[aa[,1]==TRUE,,drop=FALSE]
gexpr.m_h_1 <- subset(gexpr.m_h,new_group_2 %in% rownames(aa))


dim(gexpr.m_h_1)
# [1] 14198  3054

dim(gexpr.m_h)
# [1] 14198  3108


gexpr.m_h_1$new_group <- paste0(gexpr.m_h_1$species,"_",gexpr.m_h_1$subtype)
markers_result<- data.frame(matrix(ncol = 8,nrow=0))
colnames(markers_result) <- c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","gene","region")
for (i in c("FC","MSC","TC","OC")){
  gexpr_t <- subset(gexpr.m_h_1,subset=lobe==i)
  Idents(gexpr_t) <- "new_group"
  markers <- FindAllMarkers(gexpr_t, min.pct = 0.1, logfc.threshold = log(1.25),only.pos = T)
  markers1 <- markers[which(markers$p_val_adj < 0.01),]
  gene_1st <- grep('^RPL|^RPS', rownames(markers1), value=T ,invert=T)
  markers1 <- markers1[gene_1st,]
  markers1$region <- i
  markers_result <- rbind(markers_result,markers1)
}

markers_final<- data.frame(matrix(ncol = 8,nrow=0))
colnames(markers_final) <- c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","gene","region")
for (n in c("FC","MSC","TC","OC")){
  markers_result1 <- subset(markers_result,subset=region==n)
  for (m in unique(markers_result1$cluster)){
    gene_lst1 <- unique(markers_result1$gene[which(markers_result1$cluster==m)])
    gene_lst2 <- unique(markers_result$gene[which(markers_result$cluster!=m)])
    gene_lst3 <- setdiff(gene_lst1, gene_lst2)
    markers_result2 <- markers_result1[which(markers_result1$gene %in% gene_lst3),]
    markers_result2 <- markers_result2[order(markers_result2$avg_log2FC,decreasing = TRUE),]
    markers_result2 <- head(markers_result2,10)
    markers_final <- rbind(markers_final,markers_result2)
  }
}



unique(gexpr.m_h_1$subtype[gexpr.m_h_1$species=="Human"])
unique(gexpr.m_h_1$subtype[gexpr.m_h_1$species=="Monkey"])

lobe1 <- c("FC","MSC","TC","OC")
celltype1 <- c("inIPC ASCL1 DLX1","inIPC ASCL1 GADD45G","inIPC ASCL1 NKX2-1","inIPC ASCL1 SP8 PAX6")
celltype2 <- c("inIPC ASCL1 DLX1","inIPC ASCL1 NKX2-1","inIPC ASCL1 SP8 PAX6")
new_group_lst <- c()
for (n in (lobe1)){
  for (m in celltype1){
    p = paste0(n,"_Human_",m)
    new_group_lst <- c(new_group_lst,p)
  }
  for (m in celltype2){
    p = paste0(n,"_Monkey_",m)
    new_group_lst <- c(new_group_lst,p)
  }
}


gexpr.m_h_1@meta.data$new_group_2 <- factor(gexpr.m_h_1@meta.data$new_group_2,levels = new_group_lst[sort(match(unique(gexpr.m_h_1$new_group_2),new_group_lst))])







markers_final$group=paste(markers_final$region,markers_final$cluster,sep = "_")
markers_final$group <- factor(markers_final$group,levels = new_group_lst[sort(match(unique(markers_final$group),new_group_lst))])
markers_final=markers_final %>% arrange(group)
markers_final_1=markers_final %>% distinct(gene,.keep_all = T)




uni_gene <- markers_final_1$gene
gene_cell_exp <- AverageExpression(gexpr.m_h_1,
                                   features = unique(uni_gene),
                                   group.by = 'new_group_2',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene_cell_exp1 <- gene_cell_exp %>% rownames_to_column(var = "gene")
write.xlsx(gene_cell_exp1,file = "fig2F.Expression.xlsx")
library(openxlsx)
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))


marker_exp[marker_exp > 2]=2
marker_exp[marker_exp < -2]=-2

pdf("fig2C.whole_PerLobe_heatmap_forGO.pdf",6,12)
p1 <- Heatmap(marker_exp,
              show_column_names = T,
              show_row_names = T,
              cluster_rows = F,
              cluster_columns = F,
              #column_order = new_group_lst1,
              name = "Z score",
              row_names_side =  'left',
              col = colorRampPalette(c("#a0d8ef","white","#d9333f"))(100),
              column_title = NULL,
              row_title = NULL,
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 5),
              show_heatmap_legend=T
              
)
p1
dev.off()


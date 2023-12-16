# Expression Analysis CUP/EAE
# 06.01.2020

# Libraries and pheno data ####
library(beadarray)
library(lumi)
library(illuminaMousev2.db)
library(RColorBrewer)
library(openxlsx)
library(limma)

anno <- read.xlsx("./Cup_EAE_SampleSheet.xlsx")
anno <- anno[order(anno$SentrixID),]
anno <- anno[anno$Sample.ID!="33",] # Remove outlier Sample ID 33

# Read files, normalize and filter ####
idatFiles <- list.files(path = "./idat", pattern = ".idat", full.names=TRUE)
BSData <- readIdatFiles(idatFiles)

# QC
boxplot(as.data.frame(log2(exprs(BSData))), las = 2, outline = FALSE, ylab = "log2(intensity)")
boxplot(as.data.frame(nObservations(BSData)), las = 2, outline = FALSE, ylab = "number of beads")

# Norm
BSData.quantile = normaliseIllumina(BSData, method = "quantile", transform = "log2")
#BSData.vsn = normaliseIllumina(BSData, method = "vsn", transform = "log2")
#BSData.neqc = normaliseIllumina(BSData, method = "neqc")
# QC2
par(mfrow = c(1, 2))
boxplot(as.data.frame(log2(exprs(BSData))), las = 2, outline = FALSE, ylab = "log2(intensity)")
boxplot(as.data.frame((exprs(BSData.quantile))), las = 2, outline = FALSE, ylab = "log2(intensity)")
# Filtering
ids <- as.character(featureNames(BSData.quantile))
qual <- unlist(mget(ids, illuminaMousev2PROBEQUALITY, ifnotfound=NA))
table(qual)
rem <- qual == "No match" | qual == "Bad" | is.na(qual)
BSData.quantile.flt <- BSData.quantile[!rem,]

# Expression of genes of interest -----------------------------------------
library(dplyr)
BSData.quantile.flt <- addFeatureData(BSData.quantile.flt, toAdd = c("SYMBOL","GENOMICLOCATION"))

#genes<-read.xlsx("./Genliste_Array.xlsx")

# GO Tert "myelin sheath"
# GO:0043209 (mouse)
genes <- read.xlsx("./GO_term_summary_20220518_070401.xlsx")
genes <- genes[genes$Annotated.Term == "myelin sheath",]
genes <- genes[genes$Evidence=="IDA",]
genes <- unique(genes[,2])
length(genes)

#all_genes<-as.character(genes$all)

df<-data.frame(exprs(BSData.quantile.flt))
df$gene<-fData(BSData.quantile.flt)[,"SYMBOL"]

df<-df[is.na(df$gene)==FALSE,]
df1 = as.data.frame(df %>% dplyr::group_by(gene) %>% dplyr::summarize_all(mean))

all<-df1[df1$gene%in%genes,] # or: genes$Gene.Symbol

mat<-as.matrix(all[,2:length(colnames(all))])
rownames(mat)<-all$gene

annotation <- data.frame(sample=as.character(anno$group),
                         model=as.character(anno$model)) 
#weeks=as.character(anno1$weeks), 
#batch=as.character(anno1$SentrixID))
row.names(annotation) <- colnames(mat)
mat_colors <- list(sample = brewer.pal(6, "Set1"),
                   model = brewer.pal(3,"Dark2")[1:2])

names(mat_colors$sample) <- unique(anno$group)
names(mat_colors$model) <- unique(anno$model)

my_breaks <- c(seq(7, 9, by=0.05)) 
my_colors <- c(colorRampPalette(colors = c("blue", "yellow"))(length(my_breaks)))

sd=apply(mat,1,sd)
m.sd <- names(sort(sd, decreasing = TRUE))
mat1=mat[rev(order(sd))[1:25],]

pheatmap(mat,
         clustering_method = "average",
         annotation = annotation,
         annotation_colors = mat_colors,
         breaks = my_breaks,
         color = my_colors, 
         fontsize_row = 8,
         show_rownames = T, show_colnames = F)

# w/o clustering
annotation <- annotation[order(annotation$model),]
annotation <- annotation[order(annotation$sample),]
mat_ordered <- mat[,rownames(annotation)]

my_breaks <- c(seq(4, 10, by=0.05)) 
my_colors <- c(colorRampPalette(colors = c("blue", "yellow"))(length(my_breaks)))

pheatmap(mat_ordered,
         cluster_cols = FALSE,
         annotation = annotation,
         annotation_colors = mat_colors,
         legend_breaks = c(4, 6, 8, 10),
         legend_labels = c("4", "6", "8", "Expression"),
         legend = TRUE,
         main = "Normalized Expression of\nMyelin-associated Genes (GO:0043209)",
         fontsize = 11,
         breaks = my_breaks,
         color = my_colors, 
         fontsize_row = 8,
         show_rownames = T, show_colnames = F)


# Differential Expression Analysis ----------------------------------------
anno <- read.xlsx("./Cup_EAE_SampleSheet.xlsx")
anno <- anno[order(anno$SentrixID),]
anno <- anno[anno$Sample.ID!="33",] # Remove outlier Sample ID 33

# Read files, normalize and filter ####
idatFiles <- list.files(path = "./idat", pattern = ".idat", full.names=TRUE)
BSData <- readIdatFiles(idatFiles)

# Norm
BSData.quantile = normaliseIllumina(BSData, method = "quantile", transform = "log2")

# Filtering
ids <- as.character(featureNames(BSData.quantile))
qual <- unlist(mget(ids, illuminaMousev2PROBEQUALITY, ifnotfound=NA))
table(qual)
rem <- qual == "No match" | qual == "Bad" | is.na(qual)
BSData.quantile.flt <- BSData.quantile[!rem,]

design  <- model.matrix(~0 + anno$model)
colnames(design) <- c("CUP","CUP_EAE")
aw <- arrayWeights(BSData.quantile.flt, design)
fit  <- lmFit(BSData.quantile.flt, design , weights = aw)
contrasts  <- makeContrasts(CUP_EAE - CUP, levels = design)
contr.fit  <- eBayes(contrasts.fit(fit , contrasts ))
topTable(contr.fit , coef = 1)
par(mfrow = c(1, 1))
volcanoplot(contr.fit , main = "CUP_EAE vs. CUP")

library(illuminaMousev2.db) # load mouse array annotation
ids  <- as.character(contr.fit$genes$ArrayAddressID)
ids2  <- unlist(mget(ids , revmap(illuminaMousev2ARRAYADDRESS),ifnotfound = NA))

chr <- mget(ids2 , illuminaMousev2CHR , ifnotfound = NA)
refseq  <- mget(ids2 , illuminaMousev2REFSEQ , ifnotfound = NA)
entrezid  <- mget(ids2 , illuminaMousev2ENTREZID , ifnotfound = NA)
symbol  <- mget(ids2 , illuminaMousev2SYMBOL , ifnotfound = NA)
genename  <- mget(ids2 , illuminaMousev2GENENAME , ifnotfound = NA)
anno_genes  <- data.frame(Ill_ID = ids2 , Chr = as.character(chr),
                    EntrezID = as.numeric(entrezid), 
                    #RefSeq = as.character(refseq),
                    Symbol = as.character(symbol),
                    Name = as.character(genename ))
contr.fit$genes  <- anno_genes
top_results <- topTable(contr.fit,adjust.method="BH",p.value = 0.05,number=Inf)
write.xlsx(as.data.frame(top_results), file = "CUP_EAE_vs_CUP.xlsx")

# plot top results
top_genes <- unique(top_results$Symbol)
top_genes <- top_genes[-2] # remove NA outlier (cannot be ommited with na.omit)
length(top_genes)

BSData.quantile.flt <- addFeatureData(BSData.quantile.flt, toAdd = c("SYMBOL","GENOMICLOCATION"))

CUP <- list(0)
CUP_EAE <- list(0)

for (i in 1:100) { #1:length(top_genes)
  x <- which(fData(BSData.quantile.flt)[, "SYMBOL"] == as.character(top_genes[i]))
  CUP[i] <- mean(as.numeric(colMeans(exprs(BSData.quantile.flt[x,anno$model=="Cup"]))))
  CUP_EAE[i] <- mean(as.numeric(colMeans(exprs(BSData.quantile.flt[x,anno$model=="Cup_EAE"]))))
}

df <- data.frame(as.numeric(CUP), as.numeric(CUP_EAE))
colnames(df) <- c("CUP","CUP_EAE")
rownames(df) <- as.character(top_genes[1:100])
#df <- scale(df, center = F)

mat <- as.matrix(df)
mat <- mat[1:25,]

pheatmap(mat, 
         color = colorRampPalette(c("blue","yellow"))(50),
         scale = "none",
         fontsize = 8,
         cellwidth = 14,
         cellheight = 10,
         treeheight_row = 0,
         cluster_rows = T, 
         cluster_cols = F,
         show_rownames = T, 
         show_colnames = T)


# Enrichment Analysis -----------------------------------------------------
library(clusterProfiler); library(enrichplot); library(ggplot2)
geneList <- as.numeric(top_results[,6])
names(geneList) <- as.character(top_results[,3])
unique <- unique(names(geneList))
geneList_unique <- geneList[na.omit(unique)]
geneList_unique <- sort(geneList_unique, decreasing = T)
gene.df <- bitr(names(geneList_unique), fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Mm.eg.db)
head(gene.df)

ego <- enrichGO(gene = gene.df$ENTREZID,
                keyType = "ENTREZID",
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",   
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
#ego@result$Description %<>% stringr::str_sub(., 1, 50)
dotplot(ego, showCategory=10, font.size=10)
#emapplot(ego, showCategory = 10, orderBy="p.adjust", font.size=8)
head(ego@result$Description %>% stringr::str_sub(., 1, 75))
ego@result$Description %<>% stringr::str_sub(., 1, 75)
ego@result$p.adjust %<>% signif(.,digits=2)
ego <- pairwise_termsim(ego)
treeplot(ego, fontsize=4, showCategory = 20, label_format = 20, cex_category=0.6, offset=1.3) +
  theme(legend.position="right") + ggtitle("Hierarchical Clustering of enriched GO Terms") 
x <- treeplot(ego, fontsize=4, showCategory = 20, label_format = 20, cex_category=0.6, offset=1.3) +
  theme(legend.position="right") + ggtitle("Hierarchical Clustering of enriched GO Terms")

# How many of the n=393 genes are related to inflammatory terms?
plotted_terms <- x$data$label[1:20]
head(plotted_terms)
tmp <- as.data.frame(ego)
head(tmp$Description)
tmp <- tmp[match(plotted_terms,tmp$Description),] # select terms from tree plot
tmp <- tmp$geneID
tmp <- strsplit(tmp,split="/")
tmp <- unlist(tmp)
tmp <- unique(tmp)
tmp %in% top_results$Symbol
length(tmp)

ego1 <- gseGO(gene = geneList_unique,
              keyType = "ENTREZID",
              OrgDb = org.Mm.eg.db,
              ont = "BP",
              pvalueCutoff = 0.05)
dotplot(ego1, showCategory=10)



egox <- setReadable(ego, 'org.Mm.eg.db', 'ENTREZID')
cnetplot(egox, foldChange=geneList_unique)
cnetplot(egox, foldChange=geneList_unique, circular = TRUE, colorEdge = TRUE)
heatplot(egox, foldChange=geneList_unique)




kegg <- enrichKEGG(gene = gene.df$ENTREZID,
                   organism = 'mmu',
                   pvalueCutoff = 0.05)
dotplot(kegg, showCategory=20)


# Cluster differentially expressed genes ----------------------------------
df<-data.frame(exprs(BSData.quantile.flt))
df$gene<-fData(BSData.quantile.flt)[,"SYMBOL"]

df<-df[is.na(df$gene)==FALSE,]
df1 = as.data.frame(df %>% dplyr::group_by(gene) %>% dplyr::summarize_all(mean))

all<-df1[df1$gene%in%tmp,] # or: genes$Gene.Symbol

mat<-as.matrix(all[,2:length(colnames(all))])
rownames(mat)<-all$gene

annotation <- data.frame(sample=as.character(anno$group),
                         model=as.character(anno$model)) 
#weeks=as.character(anno1$weeks), 
#batch=as.character(anno1$SentrixID))
row.names(annotation) <- colnames(mat)
mat_colors <- list(sample = brewer.pal(6, "Set1"),
                   model = brewer.pal(3,"Dark2")[1:2])

names(mat_colors$sample) <- unique(anno$group)
names(mat_colors$model) <- unique(anno$model)

my_breaks <- c(seq(7, 9, by=0.05)) 
my_colors <- c(colorRampPalette(colors = c("blue", "yellow"))(length(my_breaks)))

sd=apply(mat,1,sd)
m.sd <- names(sort(sd, decreasing = TRUE))
mat1=mat[rev(order(sd))[1:25],]

pheatmap(mat,
         main = "Hierarchical Clustering of immune-related genes",
         clustering_method = "average",
         annotation = annotation,
         annotation_colors = mat_colors,
         breaks = my_breaks,
         color = my_colors, 
         fontsize_row = 8,
         show_rownames = FALSE, 
         show_colnames = FALSE)


# Cluster Genes with Inflammatory GO Term ---------------------------------
df<-data.frame(exprs(BSData.quantile.flt))
df$gene<-fData(BSData.quantile.flt)[,"SYMBOL"]

df<-df[is.na(df$gene)==FALSE,]
df1 = as.data.frame(df %>% dplyr::group_by(gene) %>% dplyr::summarize_all(mean))

# inflammatory response (GO:0006954)
tmp <- read.xlsx("./GO_term_summary_20220519_045802.xlsx")
tmp <- unique(tmp$Symbol)
all<-df1[df1$gene%in%tmp,] # or: genes$Gene.Symbol

mat<-as.matrix(all[,2:length(colnames(all))])
rownames(mat)<-all$gene

annotation <- data.frame(sample=as.character(anno$group),
                         model=as.character(anno$model)) 
#weeks=as.character(anno1$weeks), 
#batch=as.character(anno1$SentrixID))
row.names(annotation) <- colnames(mat)
mat_colors <- list(sample = brewer.pal(6, "Set1"),
                   model = brewer.pal(3,"Dark2")[1:2])

names(mat_colors$sample) <- unique(anno$group)
names(mat_colors$model) <- unique(anno$model)

my_breaks <- c(seq(3, 9, by=0.05)) 
my_colors <- c(colorRampPalette(colors = c("blue", "yellow"))(length(my_breaks)))

#sd=apply(mat,1,sd)
#m.sd <- names(sort(sd, decreasing = TRUE))
#mat1=mat[rev(order(sd))[1:25],]

pheatmap(mat,
         main = "Hierarchical Clustering of inflammatory response (GO:0006954)",
         clustering_method = "average",
         annotation = annotation,
         annotation_colors = mat_colors,
         breaks = my_breaks,
         color = my_colors, 
         fontsize_row = 8,
         show_rownames = FALSE, 
         show_colnames = FALSE)


# differential expression analysis week 0 ---------------------------------
anno <- read.xlsx("./Cup_EAE_SampleSheet.xlsx")
anno <- anno[anno$weeks==0,]

# Read files, normalize and filter ####
idatFiles <- list.files(path = "./idat", pattern = ".idat", full.names=TRUE)
BSData <- readIdatFiles(idatFiles)
BSData <- BSData[,grepl(anno$SentrixID[1], BSData$sectionNames)]

# Norm
BSData.quantile = normaliseIllumina(BSData, method = "quantile", transform = "log2")

# Filtering
ids <- as.character(featureNames(BSData.quantile))
qual <- unlist(mget(ids, illuminaMousev2PROBEQUALITY, ifnotfound=NA))
table(qual)
rem <- qual == "No match" | qual == "Bad" | is.na(qual)
BSData.quantile.flt <- BSData.quantile[!rem,]

design  <- model.matrix(~0 + anno$model)
colnames(design) <- c("CUP","CUP_EAE")
aw <- arrayWeights(BSData.quantile.flt, design)
fit  <- lmFit(BSData.quantile.flt, design , weights = aw)
contrasts  <- makeContrasts(CUP_EAE - CUP, levels = design)
contr.fit  <- eBayes(contrasts.fit(fit , contrasts ))
topTable(contr.fit , coef = 1)
par(mfrow = c(1, 1))
volcanoplot(contr.fit , main = "CUP_EAE vs. CUP")

library(illuminaMousev2.db) # load mouse array annotation
ids  <- as.character(contr.fit$genes$ArrayAddressID)
ids2<-unlist(mget(ids,revmap(illuminaMousev2ARRAYADDRESS),ifnotfound = NA))

chr <- mget(ids2 , illuminaMousev2CHR , ifnotfound = NA)
refseq  <- mget(ids2 , illuminaMousev2REFSEQ , ifnotfound = NA)
entrezid  <- mget(ids2 , illuminaMousev2ENTREZID , ifnotfound = NA)
symbol  <- mget(ids2 , illuminaMousev2SYMBOL , ifnotfound = NA)
genename  <- mget(ids2 , illuminaMousev2GENENAME , ifnotfound = NA)
anno_genes  <- data.frame(Ill_ID = ids2 , Chr = as.character(chr),
                          EntrezID = as.numeric(entrezid), 
                          #RefSeq = as.character(refseq),
                          Symbol = as.character(symbol),
                          Name = as.character(genename ))
contr.fit$genes  <- anno_genes
top_results <- topTable(contr.fit,adjust.method="BH",p.value = 0.05,number=Inf)
write.xlsx(as.data.frame(top_results), file = "CUP_EAE_vs_CUP_week0.xlsx")

# plot top results
top_genes <- unique(top_results$Symbol)
top_genes <- top_genes[-2] # remove NA outlier (cannot be ommited with na.omit)

BSData.quantile.flt <- addFeatureData(BSData.quantile.flt, toAdd = c("SYMBOL","GENOMICLOCATION"))

CUP <- list(0)
CUP_EAE <- list(0)

for (i in 1:100) { #1:length(top_genes)
  x <- which(fData(BSData.quantile.flt)[, "SYMBOL"] == as.character(top_genes[i]))
  CUP[i] <- mean(as.numeric(colMeans(exprs(BSData.quantile.flt[x,anno$model=="Cup"]))))
  CUP_EAE[i] <- mean(as.numeric(colMeans(exprs(BSData.quantile.flt[x,anno$model=="Cup_EAE"]))))
}

df <- data.frame(as.numeric(CUP), as.numeric(CUP_EAE))
colnames(df) <- c("CUP","CUP_EAE")
rownames(df) <- as.character(top_genes[1:100])
df <- scale(df, center = F)

mat <- as.matrix(df)
mat <- mat[1:20,]
pheatmap(mat, 
         color = colorRampPalette(c("blue", "yellow"))(50),
         fontsize = 8,
         cellwidth = 14,
         cellheight = 10,
         treeheight_row = 0,
         cluster_rows = T, 
         cluster_cols = F,
         show_rownames = T, 
         show_colnames = T)
#####

#####

# Explorative analysis ####
library(pheatmap)
d <- exprs(BSData.quantile.flt[,anno$Sample.ID!="33"]) # Remove outlier Sample ID 33
anno1 <- anno[anno$Sample.ID!="33",]
sampleDists <- dist(t(d))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(anno1$Sample_Name)
colnames(sampleDistMatrix) <- as.character(anno1$SampleID)
#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
annotation <- data.frame(sample=as.character(anno1$group),
                         model=as.character(anno1$model), 
                         weeks=as.character(anno1$weeks), 
                         batch=as.character(anno1$SentrixID))
row.names(annotation) <- colnames(sampleDistMatrix)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         clustering_method = "average",
         annotation = annotation,
         #col=colors, 
         fontsize_col = 7,
         show_rownames = F, show_colnames = T)

# tSNE analysis
library(Rtsne)
library(ggplot2)
dist <- dist(t(d))
res <- Rtsne(dist, 
             pca = FALSE, 
             verbose = TRUE, 
             is_distance = TRUE, 
             theta = 0, 
             max_iter = 2000, 
             perplexity = 3)$Y
df <- data.frame(res)

ggplot(df, aes(x=df[,1], y=df[,2])) +
  geom_point(aes(color=anno1$group),alpha= 0.7, size=3, show.legend = T) +
  #stat_ellipse(aes(fill=time), geom = "polygon", alpha=0.2, show.legend = F) +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.justification=c(0,0),
        legend.title = element_text(face = "bold"),
        legend.background = element_blank(),
        legend.position = "right",
        panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.key.size = unit(1, 'mm'))
#####


# Plot single genes -------------------------------------------------------
library(reshape2)
library(ggpubr)
df<-data.frame(exprs(BSData.quantile.flt))
df$gene<-fData(BSData.quantile.flt)[,"SYMBOL"]
df<-df[is.na(df$gene)==FALSE,]
all <- as.data.frame(df %>% dplyr::group_by(gene) %>% dplyr::summarize_all(mean))
df2<-all
head(df2)
#df2<-df2[df2$gene=="Tnf"|df2$gene=="Ifng",]
df2<-df2[df2$gene=="Il17a"|df2$gene=="Csf2",]
df2<-df2[df2$gene=="Plp1"|df2$gene=="Plp1",]

#colnames(df2)<-c("gene",paste0(anno$Sample.ID,"_",anno$model))
colnames(df2)<-c("gene",anno$Sample.ID)
dfm <- melt(df2, id=c("gene"))
dfm$group <- anno$group[match(dfm$variable, anno$Sample.ID)]
#dfm$group<-ifelse(grepl("*Cup_EAE",dfm$variable),"Cup_EAE","Cup")
head(dfm)

ggplot(dfm, aes(x=group, y=value, fill=group)) +
  geom_boxplot(alpha=0.8) +
  geom_point() +
  #scale_fill_manual(values=c("blue","orange")) +
  facet_wrap(~gene, scales = "free") +
  stat_compare_means(size=3, method="t.test",
                     comparisons = list(c("Cup", "Cup_EAE"))) +
  #ylim(0,10) +
  xlab("") + ylab("") +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank(),
  #       axis.title.y=element_blank(),
  #       axis.text.y=element_blank(),
  #       axis.ticks.y=element_blank())

m1<-df1[df1$gene%in%genes$M1,]
m2a<-df1[df1$gene%in%genes$M2a,]
m2b<-df1[df1$gene%in%genes$M2b,]
m2c<-df1[df1$gene%in%genes$M2c,]
m2d<-df1[df1$gene%in%genes$M2d,]
act<-df1[df1$gene%in%genes$activation,]


colnames(mat) <- colnames(as.matrix(exprs(BSData.quantile.flt)))
annotation <- data.frame(age=as.character(anno$Age), condition=as.character(anno$Condition), batch=as.character(anno$Sentrix_ID))
row.names(annotation) <- colnames(mat)
pheatmap(mat,
         annotation = annotation,
         #col=colors, 
         fontsize_col = 7,
         show_rownames = T, show_colnames = F)

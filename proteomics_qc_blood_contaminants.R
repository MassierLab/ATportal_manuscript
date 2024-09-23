##### Libraries #####

library(tidyverse)
library(correlation)
library(DataExplorer)
library(data.table)
library(NbClust)
library(tidymodels)
library(ggrepel)
library(GGally)
library(ggcorrplot)
library(tidyselect)
library(embed)
library(mclust)
library(pheatmap)
library(stats)
library(factoextra)
library(FactoMineR)
library(tidyr)
library(rstatix)
library(ggpubr)
library(edgeR)
library(readxl)


##### Functions #####

stderror <- function(x) sd(x, na.rm = T)/sqrt(length(x))
maxmin <- function(x, na.rm=TRUE){
  if(is.vector(x)==TRUE){
    maxs <- max(x, na.rm = T)
    mins <- min(x, na.rm = T)
    scale(x,center=mins,scale=maxs-mins)
  } else {
    maxs <- apply(x, 2, max)
    mins <- apply(x, 2, min)
    scale(x, center = mins, scale = maxs - mins)
  }
}

portalcol <- colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(51)
portalcol2 <- c("#1a4659", "#e2c744", "#ba5c28", "#4f736a", "#f0cbbf", "#b2dfda", "#f59b7c", "#7dbfb3", "#939a64", "#8ca5a5", "#adbec4", "#c3b6b6")
portalcol3 <- colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(25)

Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))}}
Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))}}
Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)}

##### Raw data #####


# Clinical data

clin.portal <- read.csv("sample_list_portal.csv")
clin.portal$ExpID <- as.character(clin.portal$ExpID)

# Proteins
raw <- read.csv("WAT-panel-full_Report_PG.csv", sep = ",")
annot <- raw[,1:7]
colnames(annot) <-  sub("PG.", "", colnames(annot))

# Remove MaxQuant contaminants

raw <- raw %>% filter(! grepl("MaxQuant Contaminants", annot$FastaFiles))

# Clean protein file

pr.int <- raw[,c(1,8:107)]

colnames(pr.int) <-  sub(".htrms.PG.Quantity", "", colnames(pr.int))
colnames(pr.int) <-  sub("*K.LM_WAT.panel_", "", colnames(pr.int))
colnames(pr.int) <-  sub(".*F", "", colnames(pr.int))
colnames(pr.int) <-  sub("\\_.*", "", colnames(pr.int))
colnames(pr.int) <-  sub("PG.", "", colnames(pr.int))

pr.bc <- raw[,c(1,8:107)] # this is to check for batch effects between different days of measurement

colnames(pr.bc) <-  sub(".htrms.PG.Quantity", "", colnames(pr.bc))
colnames(pr.bc) <-  sub("*K.LM_WAT.panel_", "", colnames(pr.bc))
colnames(pr.bc) <-  sub("\\_.*", "", colnames(pr.bc))
colnames(pr.bc) <-  sub("PG.", "", colnames(pr.bc))
colnames(pr.bc) <-  sub("..20", "_", colnames(pr.bc))
colnames(pr.bc) <-  sub(".*?_", "", colnames(pr.bc))
colnames(pr.bc)[21] <-  sub("..20", "", colnames(pr.bc)[21])


pr.ln <- pr.int %>% pivot_longer(cols = 2:101, names_to = "ExpID", values_to = "Intensity") 
pr.ln$ExpID <- gsub("^0", "", pr.ln$ExpID)

pr.ln$Intensity[is.nan(pr.ln$Intensity)]<-NA
pr.ln$Intensity[is.infinite(pr.ln$Intensity)]<-NA

plot_intro(pr.ln)

##### QC #####


# Batch effect


pr.bc.ln <- pr.bc %>% pivot_longer(cols = 2:101, names_to = "Date", values_to = "Intensity") 

pr.bc.ln <- pr.bc %>% pivot_longer(cols = 2:ncol(pr.bc), names_to = "Date", values_to = "Intensity")
pr.bc.ln$Intensity[is.nan(pr.bc.ln$Intensity)]<-NA
pr.bc.ln$Intensity[is.infinite(pr.bc.ln$Intensity)]<-NA

ggplot(pr.bc.ln, aes(x=as.character(Date), y=Intensity)) + geom_violin()


stderror <- function(x) sd(x, na.rm = T)/sqrt(length(x))

pr.bc.stats <- pr.bc.ln %>% group_by(Date) %>%
  summarise(Min = min(Intensity,na.rm = TRUE),
            Q1 = quantile(Intensity,probs = .25,na.rm = TRUE),
            Median = median(Intensity, na.rm = TRUE),
            Q3 = quantile(Intensity,probs = .75,na.rm = TRUE),
            Max = max(Intensity,na.rm = TRUE),
            Mean = mean(Intensity, na.rm = TRUE),
            SD = sd(Intensity, na.rm = TRUE),
            CV = abs((sd(Intensity, na.rm = TRUE)/mean(Intensity, na.rm = TRUE)*100)),
            MAD = mad(Intensity, na.rm = TRUE),
            SEM = stderror(Intensity),
            iqr = IQR(Intensity, na.rm = TRUE),
            Var = var(Intensity, na.rm = TRUE),
            Total = sum(Intensity, na.rm = T),
            n=n(),
            samples = length(unique(Date))) 

ggplot(pr.bc.stats, aes(x=as.character(Date), y=Median)) + geom_col()
ggplot(pr.bc.stats, aes(x=as.character(Date), y=Total)) + geom_col()
ggplot(pr.bc.stats, aes(x=as.character(Date), y=Total/samples)) + geom_col()

## Count per sample

prt_count <- pr.ln %>% group_by(ExpID) %>%
  summarise(count = sum(!is.na(Intensity)))

ggplot(prt_count, aes(x=reorder(ExpID, count), y=count)) + geom_col(fill = "#21918c") +
  geom_hline(yintercept = 3000, linetype = "dotted") +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1))

low.count.out <- prt_count %>% filter(count < 3000)
low.count.out <- low.count.out$ExpID

## Count per protein


prt_count <- pr.ln %>% group_by(ExpID) %>%
  summarise(count = sum(!is.na(Intensity)))

ggplot(prt_count, aes(x=reorder(ExpID, count), y=count)) + geom_col(fill = "#21918c") +
  geom_hline(yintercept = 3000, linetype = "dotted") +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1))

low.count.out <- prt_count %>% filter(count < 3000)
low.count.out <- low.count.out$ExpID


prt_count_prt <- pr.ln %>% group_by(ProteinGroups) %>%
  summarise(count = sum(!is.na(Intensity)))

low_prt <- prt_count_prt %>% filter(count <= 25)
low_prt <- low_prt$ProteinGroups

ggplot(prt_count_prt, aes(x=reorder(ProteinGroups, count), y=count)) + geom_col(fill = "#21918c") +
  geom_hline(yintercept = 25, linetype = "dotted") +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_blank())

## Clean and log tranfrom


pr.ln <- left_join(pr.ln, clin.portal[,c(1,4)]) 

pr.ln.cl <- pr.ln[! pr.ln$ExpID %in% low.count.out,]
pr.ln.cl <- pr.ln.cl[! pr.ln.cl$ProteinGroups %in% low_prt,]

pr.ln.cl <- pr.ln.cl %>% mutate(log2Int = log2(Intensity),
                                BMI_group = case_when(BMI < 25 ~ "BMI < 25",
                                                      BMI >= 30 ~ "BMI >= 30",
                                                      TRUE ~ "BMI 25 - 29.9"))


# Corelation Matrix 

pr.cor <- pr.ln.cl 
pr.cor$ID <- paste(pr.cor$ExpID, pr.cor$BMI_group, sep="_")
colnames(pr.cor)
pr.cor <- pr.cor[,c(1,5,7)] %>% pivot_wider(names_from = ID, values_from = log2Int) %>% as.data.frame()
rownames(pr.cor) <- pr.cor$ID

corr_mat <- round(cor(pr.cor[,-1], use="pairwise.complete.obs"),2)  
pheatmap(corr_mat)

## PCA

pca.imp <- pr.ln.cl %>% group_by(ProteinGroups, BMI_group) %>% 
  mutate_at('log2Int', ~replace_na(., median(., na.rm=T))) %>%
  ungroup()

colnames(pca.imp)

pr.pca <- pca.imp[,c(1,2,5,6)] %>% 
  pivot_wider(names_from = ProteinGroups, values_from = log2Int) %>% as.data.frame()
pr.pca$ID <- paste(pr.pca$ExpID, pr.pca$BMI_group)
pr.pca <- pr.pca %>% select(ID, everything())
rownames(pr.pca) <- pr.pca$ID
pr.pca$obese <- as.factor(pr.pca$BMI_group)
pr.pca$ExpID <- as.factor(pr.pca$ExpID)

pca_rec <- recipe(~., data = pr.pca) %>%
  update_role(ExpID, BMI_group, new_role = "BMI_group") %>%
  step_pca(all_numeric())

pca_prep <- prep(pca_rec)
juice(pca_prep) %>%
  ggplot(aes(PC1, PC2, label = ExpID)) +
  geom_point(aes(color = as.factor(BMI_group)), alpha = 0.7, size = 2) + 
  stat_conf_ellipse(aes(color = BMI_group)) +
  geom_text(check_overlap = TRUE, hjust = "inward") +
  labs(color = NULL) +
  theme_classic()

##### Blood contaminants ####

at.markers <- c("ADIPOQ", "PLIN1", "PLIN4", "LEP", "ABHD1", "CIDEA", "PDLIM1", "SLC2A4", "ABHD5", "LIPE", "PNPLA2", "FABP4")
at.markers <- annot[annot$Genes %in% at.markers, c(1,2)]
at.markers$CellType <- "AT"

bl.markers <- read.csv("mann_blood.csv") %>% drop_columns(2) %>% left_join(annot[,1:2])

meta.markers <- read.csv("meta_celltypes.csv") %>% left_join(annot[,1:2])

at.markers <- at.markers[,c(2,3,1)]

markers <- rbind(at.markers, meta.markers, bl.markers)
markers <- markers[!duplicated(markers$ProteinGroups),]

pr.cor <- pr.ln.cl[pr.ln.cl$ProteinGroups %in% bl.markers$ProteinGroups,]
pr.cor$ID <- paste(pr.cor$ExpID, pr.cor$Comments, sep="_")
colnames(pr.cor)
pr.cor <- pr.cor[,c(1,10,9)] %>% pivot_wider(names_from = ID, values_from = log2Int) %>% as.data.frame()
rownames(pr.cor) <- pr.cor$ID

corr_mat <- round(cor(pr.cor[,-1], use="pairwise.complete.obs"),2)  
pheatmap(corr_mat)

pr.at <- pr.ln.cl %>% left_join(annot[,1:2]) %>% filter(Genes %in% at.markers$Genes) 

pr.at %>% ggplot(aes(x=reorder(ExpID, log2(Intensity)), y=Genes, fill=log2(Intensity))) + geom_tile() + 
  scale_fill_viridis() + theme_classic()
pr.at %>% filter(Genes == "FABP4") %>% 
  ggplot(aes(x=reorder(ExpID, log2(Intensity)), y=log2(Intensity), fill=Comments)) + geom_col() + theme_classic()

pr.bl <- pr.ln.cl %>% left_join(annot[,1:2]) %>% filter(Genes %in% bl.markers$Genes) 

pr.bl %>% ggplot(aes(x=reorder(ExpID, log2(Intensity)), y=Genes, fill=log2(Intensity))) + geom_tile() + 
  scale_fill_viridis() + theme_classic()
pr.bl %>% filter(Genes == "HBB") %>% 
  ggplot(aes(x=reorder(ExpID, log2(Intensity)), y=log2(Intensity), fill=Comments)) + geom_col() + theme_classic()

pr.mt <- pr.ln.cl %>% left_join(annot[,1:2]) %>% filter(Genes %in% meta.markers$Genes) 

pr.mt %>% ggplot(aes(x=reorder(ExpID, log2(Intensity)), y=Genes, fill=log2(Intensity))) + geom_tile() + 
  scale_fill_viridis() + theme_classic()
pr.mt %>% filter(Genes == "FBN1") %>% 
  ggplot(aes(x=reorder(ExpID, log2(Intensity)), y=log2(Intensity), fill=Comments)) + geom_col() + theme_classic()


# Protein correlations - Markers only

pr.cor.markers <- pr.ln.cl[pr.ln.cl$ProteinGroups %in% markers$ProteinGroups, ]
pr.cor.markers$ID <- paste(pr.cor.markers$ExpID, pr.cor.markers$Comments, sep="_")
colnames(pr.cor.markers)
pr.cor.markers <- pr.cor.markers[,c(1,10,9)] %>% pivot_wider(names_from = ProteinGroups, values_from = log2Int) %>% as.data.frame()
rownames(pr.cor.markers) <- pr.cor.markers$ID

corr.mat.markers <- round(cor(pr.cor.markers[,-1], use="pairwise.complete.obs", method = "spearman"),2) %>% as.data.frame()
pheatmap(as.matrix(corr.mat.markers), 
         color = colorRampPalette(c(portalcol))(20), 
         clustering_method = "ward.D2",
         show_rownames = F, show_colnames = F, breaks = breaksList)

cor.markers <- correlation(pr.cor.markers[-1], method = "spearman") %>% as.data.frame()

# Clustering

dist.markers <- dist(corr.mat.markers, method = 'euclidean')

set.seed(123)
hclust.markers <- hclust(dist.markers, method = 'ward.D2')
plot(hclust.markers, labels = F)
rect.hclust(hclust.markers , h = 11, border = 2:6)

hcl.markers.grp <- cutree(hclust.markers, h = 11)
table(hcl.markers.grp)
hcl.markers.grp <- cutree(hclust.markers, h = 10) %>% as.data.frame() %>% rownames_to_column("ProteinGroups") %>%
  left_join(markers)
colnames(hcl.markers.grp)[2] <- c("Cluster")
rownames(hcl.markers.grp) <- hcl.markers.grp$ProteinGroups

breaksList = seq(-1, 1, by = 0.1)

cl1 <- corr.mat.markers[rownames(corr.mat.markers) %in% rownames(hcl.markers.grp)[hcl.markers.grp$Cluster == 1],]
pheatmap(as.matrix(cl1), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl2 <- corr.mat.markers[rownames(corr.mat.markers) %in% rownames(hcl.markers.grp)[hcl.markers.grp$Cluster == 2],]
pheatmap(as.matrix(cl2), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl3 <- corr.mat.markers[rownames(corr.mat.markers) %in% rownames(hcl.markers.grp)[hcl.markers.grp$Cluster == 3],]
pheatmap(as.matrix(cl3), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl4 <- corr.mat.markers[rownames(corr.mat.markers) %in% rownames(hcl.markers.grp)[hcl.markers.grp$Cluster == 4],]
pheatmap(as.matrix(cl4), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl5 <- corr.mat.markers[rownames(corr.mat.markers) %in% rownames(hcl.markers.grp)[hcl.markers.grp$Cluster == 5],]
pheatmap(as.matrix(cl5), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl6 <- corr.mat.markers[rownames(corr.mat.markers) %in% rownames(hcl.markers.grp)[hcl.markers.grp$Cluster == 6],]
pheatmap(as.matrix(cl6), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))


hcl.markers.grp <- left_join(hcl.markers.grp, annot[annot$ProteinGroups %in% hcl.markers.grp$ProteinGroups, ], multiple = "all")

final.markers <- data.frame(ProteinGroups = c(rownames(cl2), rownames(cl4), rownames(cl3))) %>% left_join(annot)


# Global marker correlation profiling

pr.markers.final <- pr.cor.markers[,names(pr.cor.markers) %in% c(rownames(cl2), rownames(cl4), rownames(cl3))] 
pr.markers.final$ID <- rownames(pr.markers.final)

pr.cor.global <- pr.ln.cl %>% left_join(annot[,1:2])
pr.cor.global$ID <- paste(pr.cor.global$ExpID, pr.cor.global$Comments, sep="_")
colnames(pr.cor.global)
pr.cor.global <- pr.cor.global[,c(1,11,9)] %>% filter(! pr.cor.global$ProteinGroups %in% colnames(pr.markers.final)) 

pr.cor.global <- pr.cor.global %>% nest(data = -ProteinGroups)

chunk_size <- 100
num_chunks <- ceiling(nrow(pr.cor.global) / chunk_size)
chunk_indices <- split(1:nrow(pr.cor.global), cut(1:nrow(pr.cor.global), num_chunks, labels = FALSE))

bl.cor.res.s <- data.frame()

for (i in 1:num_chunks) { #1:nrow(pr.cor.global)y
  start_row <- (i - 1) * chunk_size + 1
  end_row <- min(i * chunk_size, nrow(pr.cor.global))
  
  start_time <- Sys.time()
  
  dt <- pr.cor.global[start_row:end_row,]
  dt <- dt %>% unnest(data) %>% pivot_wider(names_from = ProteinGroups, values_from = log2Int) %>% left_join(pr.markers.final)
  corel <- correlation(dt[,-c(1)], method = "spearman")
  corel <- as.data.frame(corel)
  
  bl.cor.res.s <- rbind(bl.cor.res.s, corel)
  
  end_time <- Sys.time()
  processing_time <- end_time - start_time
  
  print(paste("Processed chunk", i, "of", num_chunks, "- Rows left:", nrow(pr.cor.global) - max(unlist(chunk_indices[1:i]))))
  
  print(paste("Total processing time:", processing_time))
  
  gc()
}

bl.cor.res.s.f <- distinct(bl.cor.res.s, Parameter1, Parameter2, .keep_all = T)

par1 <- bl.cor.res.s.f$Parameter1
par2 <- bl.cor.res.s.f$Parameter2

bl.cor.res.s.f2 <- bl.cor.res.s.f
bl.cor.res.s.f2$Parameter1 <- par2
bl.cor.res.s.f2$Parameter2 <- par1

rbind(bl.cor.res.s.f[,1:3], bl.cor.res.s.f2[,1:3]) %>%
  dplyr::summarise(n = dplyr::n(), .by = c(Parameter1, Parameter2)) %>%
  dplyr::filter(n > 1L) 

bl.cor.res.s.f.wd <- rbind(bl.cor.res.s.f[,1:3], bl.cor.res.s.f2[,1:3]) 
bl.cor.res.s.f.wd <- distinct(bl.cor.res.s.f.wd, Parameter1, Parameter2, .keep_all = T)
bl.cor.res.s.f.wd <- bl.cor.res.s.f.wd %>% pivot_wider(names_from = Parameter2, values_from = rho) %>% as.data.frame() #, bl.cor.res.s.f2[,1:3]
bl.cor.res.s.f.wd[is.na(bl.cor.res.s.f.wd)] <- 0
rownames(bl.cor.res.s.f.wd) <- bl.cor.res.s.f.wd$Parameter1
bl.cor.res.s.f.wd <- bl.cor.res.s.f.wd[bl.cor.res.s.f.wd$Parameter1 %in% c(rownames(cl2), rownames(cl4), rownames(cl5)),]

# Clustering

dist.markers.bl <- dist(t(bl.cor.res.s.f.wd[,-1]), method = 'euclidean')
#dist.markers.bl <- distances(bl.cor.res.s.f.wd[,-1])

set.seed(123)
hclust.markers.bl <- hclust(dist.markers.bl, method = 'ward.D2')
plot(hclust.markers.bl, labels = F)
rect.hclust(hclust.markers.bl, h = 30, border = 2:6)

hcl.markers.grp.bl <- cutree(hclust.markers.bl, h = 30)
table(hcl.markers.grp.bl)
hcl.markers.grp.bl <- cutree(hclust.markers.bl, h = 30) %>% as.data.frame() %>% rownames_to_column("ProteinGroups") %>%
  left_join(markers)
colnames(hcl.markers.grp.bl)[2] <- c("Cluster")
rownames(hcl.markers.grp.bl) <- hcl.markers.grp.bl$ProteinGroups

breaksList = seq(-1, 1, by = 0.1)

pheatmap(as.matrix(bl.cor.res.s.f.wd[,-1]), 
         #color = colorRampPalette(rev(brewer.pal(name = "PRGn")))(100),
         color = colorRampPalette(c(portalcol))(20), 
         clustering_method = "ward.D2",
         show_rownames = F, show_colnames = F, breaks = breaksList)

cl1.bl <- bl.cor.res.s.f.wd[,names(bl.cor.res.s.f.wd) %in% rownames(hcl.markers.grp.bl)[hcl.markers.grp.bl$Cluster == 1]]
pheatmap(as.matrix(cl1.bl), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl2.bl <- bl.cor.res.s.f.wd[,names(bl.cor.res.s.f.wd) %in% rownames(hcl.markers.grp.bl)[hcl.markers.grp.bl$Cluster == 2]]
pheatmap(as.matrix(cl2.bl), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl3.bl <- bl.cor.res.s.f.wd[,names(bl.cor.res.s.f.wd) %in% rownames(hcl.markers.grp.bl)[hcl.markers.grp.bl$Cluster == 3]]
pheatmap(as.matrix(cl3.bl), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl4.bl <- bl.cor.res.s.f.wd[,names(bl.cor.res.s.f.wd) %in% rownames(hcl.markers.grp.bl)[hcl.markers.grp.bl$Cluster == 4]]
pheatmap(as.matrix(cl4.bl), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl5.bl <- bl.cor.res.s.f.wd[,names(bl.cor.res.s.f.wd) %in% rownames(hcl.markers.grp.bl)[hcl.markers.grp.bl$Cluster == 5]]
pheatmap(as.matrix(cl5.bl), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))


hcl.markers.grp.bl <- left_join(hcl.markers.grp.bl, annot[annot$ProteinGroups %in% hcl.markers.grp.bl$ProteinGroups, ], multiple = "all") %>% na.omit()


blood.contaminants <- rownames(hcl.markers.grp.bl)[hcl.markers.grp.bl$Cluster == 5] 

at.proteins <- rownames(hcl.markers.grp.bl)[hcl.markers.grp.bl$Cluster == 3]

blood.contaminants <- annot[annot$ProteinGroups %in% blood.contaminants,]
at.proteins <- annot[annot$ENTREZID %in% at.proteins,]

blood.contaminants.list <- data.frame(ProteinGroups = blood.contaminants) %>% left_join(annot)


##### Clean-up and normalization #####

## Remove blood contaminants 

pr.ln.cl.f <- pr.ln.cl[! pr.ln.cl$ProteinGroups %in% c(blood.contaminants),] 
plot_intro(pr.ln.cl.f)
colnames(pr.ln.cl.f)

## Scale & Normalize

#vsn

norm.vsn <- pr.ln.cl.f[! pr.ln.cl.f$ExpID %in% samples.to.remove,1:3] %>% pivot_wider(names_from = ProteinGroups, values_from = Intensity) %>% as.data.frame()
rownames(norm.vsn) <- norm.vsn[,1]
norm.vsn <- norm.vsn[,-c(1)]
norm.vsn <- t(norm.vsn)
norm.vsn <- normalizeVSN(norm.vsn)
meanSdPlot(norm.vsn)

norm.vsn <- as.data.frame(norm.vsn)
norm.vsn$ProteinGroups <- rownames(norm.vsn)
norm.vsn <- norm.vsn %>% pivot_longer(cols = 1:92, names_to = "ExpID", values_to = "norm.vsn") 

sample_count <- pr.ln.norm[,c(1,2,3,7,8)] %>% pivot_wider(names_from = "ProteinGroups", values_from = "Intensity") 
table(sample_count$obese)
table(sample_count$obese, sample_count$sex)

contingency_table <- as.data.frame(table(sample_count$obese, sample_count$sex))
colnames(contingency_table) <- c("Obese", "Sex", "Count")

PieDonut(contingency_table, aes(Obese, Sex, count=Count), title = "Sex by Class")

##  DEP 

library(DEP)

LFQ <- pr.ln.cl.f[! pr.ln.cl.f$ExpID %in% samples.to.remove,-c(4:9)] %>% pivot_wider(names_from = ExpID, values_from = Intensity)

data_unique <- make_unique(left_join(LFQ, annot[,1:2]), "Genes", "ProteinGroups", delim = ";")
data_cols <- c(2:93)
data_annot <- data_unique[,c(1,95)]

des.ma <- clin.portal[clin.portal$ExpID %in% colnames(LFQ),c(1,24)] %>% group_by(obese) %>% 
  mutate(replicate = row_number(obese)) %>% as.data.frame()
colnames(des.ma) <- c("label", "condition", "replicate")
des.ma$replicate <- as.character(des.ma$replicate)

data.se <- make_se(data_unique, data_cols, des.ma)

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data.se)
# Filter for proteins that are identified in  replicates of at least one condition
data_filt <- filter_missval(data.se, thr = 1)
# Plot a barplot of the protein identification overlap between samples
plot_coverage(data.se)
# Plot a barplot of the number of identified proteins per samples
plot_numbers(data.se)
# Plot a heatmap of proteins with missing values
plot_missval(data_filt)
# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)
# Normalize the data
data_norm <- normalize_vsn(data_filt)
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

## Corelation Matrix 

pr.cor.f <- pr.ln.norm 
pr.cor.f$ID <- paste(pr.cor.f$ExpID, pr.cor.f$obese, sep="_")
colnames(pr.cor.f)
pr.cor.f <- pr.cor.f[,c(1,12,11)] %>% pivot_wider(names_from = ID, values_from = norm.vsn) %>% as.data.frame()
rownames(pr.cor.f) <- pr.cor.f$ID

corr_mat <- round(cor(pr.cor.f[,-1], use="pairwise.complete.obs"),2)  
pheatmap(corr_mat, color = colorRampPalette(c(portalcol))(20)) #, breaks = breaksList

## PCA

colnames(pr.ln.norm)

pca.imp <- pr.ln.norm[! pr.ln.norm$ExpID %in% samples.to.remove, ] %>% 
  group_by(ProteinGroups, obese) %>% 
  mutate_at('Intensity', ~replace_na(., median(., na.rm=T))) %>% 
  mutate_at('median.norm', ~replace_na(., median(., na.rm=T))) %>%
  mutate_at('log2Int', ~replace_na(., median(., na.rm=T))) %>%
  mutate_at('norm.vsn', ~replace_na(., median(., na.rm=T))) %>%
  ungroup()

colnames(pca.imp)

pr.pca <- pca.imp[, c(1,2,4,8,11)] %>% 
  pivot_wider(names_from = ProteinGroups, values_from = norm.vsn) %>% as.data.frame()
pr.pca$ID <- paste(pr.pca$ExpID, pr.pca$BMI_group)
pr.pca <- pr.pca %>% select(ID, everything())
rownames(pr.pca) <- pr.pca$ID
pr.pca$obese <- as.factor(pr.pca$BMI_group)
pr.pca$ExpID <- as.factor(pr.pca$ExpID)

pca_rec <- recipe(~., data = pr.pca) %>%
  update_role(ExpID, BMI_group, new_role = "BMI_group") %>%
  step_pca(all_numeric())

pca_prep <- prep(pca_rec)
summary(pca_prep$steps[[1]]$res)

juice(pca_prep) %>%
  ggplot(aes(PC1, PC2, label = ExpID)) +
  geom_point(aes(color = as.factor(BMI_group)), alpha = 0.7, size = 2) + 
  stat_conf_ellipse(aes(color = BMI_group)) +
  geom_text(check_overlap = TRUE, hjust = "inward") +
  labs(color = NULL) +
  theme_classic()

tidied_pca <- tidy(pca_prep, 1)
tidied_pca %>%
  filter(component %in% paste0("PC", 1:4)) %>%
  group_by(component) %>%
  top_n(8, abs(value)) %>%
  ungroup() %>%
  mutate(terms = reorder_within(terms, abs(value), component)) %>%
  ggplot(aes(abs(value), terms, fill = value > 0)) +
  geom_col() +
  facet_wrap(~component, scales = "free_y") +
  scale_y_reordered() +
  labs(
    x = "Absolute value of contribution",
    y = NULL, fill = "Positive?")


samples.to.remove <- c("94", "70", "45", "49")

### ICA

library(fastICA)

data_matrix <- as.matrix(pr.pca[ , -c(1:4)])

ica_result <- fastICA(data_matrix, n.comp = 2, method = "C")
ica_components <- ica_result$S
colnames(ica_components) <- c("IC1", "IC2")
total_signal <- sum(abs(ica_result$S))
contribution_IC1 <- sum(abs(ica_result$S[, 1])) / total_signal * 100
contribution_IC2 <- sum(abs(ica_result$S[, 2])) / total_signal * 100

ica_data <- cbind(pr.pca, ica_components) %>% left_join(clin[,c(1, 6)]) %>% left_join(homa_data) %>% left_join(clin_age)

ggplot(ica_data, aes(x = IC1, y = IC2, color = factor(AgeGroup))) +
  geom_point(aes(shape=sex), size = 2) +
  stat_conf_ellipse(aes(color = AgeGroup)) +
  labs(x = paste0("Independent Component 1 (", round(contribution_IC1, 1), "%)"),
       y = paste0("Independent Component 2 (", round(contribution_IC2, 1), "%)"),
       color = "Age Group") +
  theme_classic()

### CLEAN DATA

pr.ln.norm <- pr.ln.norm[! pr.ln.norm$ExpID %in% samples.to.remove,] 

write_rds(pr.ln.norm, "normalized_data.RDS")

sample.list <- pr.ln.norm[,c(2,7,8)]
sample.list <- sample.list[!duplicated(sample.list),]
colnames()


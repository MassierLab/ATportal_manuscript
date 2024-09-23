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

raw <- read.csv("WAT-panel-full_Report_PG.csv", sep = ",")
annot <- raw[,1:7]
pr.ln.norm <- read_rds("normalized_data.RDS")
clin.portal <- read.csv("sample_list_portal.csv")
clin.portal$ExpID <- as.character(clin.portal$ExpID)

##### Data descritpion #####

# Clinical

clin.desc <- clin[clin$ExpID %in% pr.ln.norm$ExpID,1:24] #clin.portal[,1:24]

table(clin.desc$sex)
table(clin.desc$obese)

clin.desc.gen <- clin.desc[,c(1,6:8,15, 24)] %>% 
  mutate(obese2 = case_when(BMI < 25 ~ "NonObese",
                            BMI >= 25 & BMI < 30 ~ "Overweight",
                            BMI >= 30 & BMI < 35 ~ "ObeseI",
                            BMI >= 35 & BMI < 40 ~ "ObeseII",
                            BMI >= 40 ~ "ObeseIII"),
         obese3 = case_when(BMI < 30 ~ "NonObese",
                            BMI >= 30 ~ "Obese"),
         obese = case_when(obese == "lean" ~ "NonObese",
                           obese == "over weight" ~ "Overweight",
                           obese == "obese" ~ "Obese")) %>%
  filter(ExpID %in% pr.ln.cl.f$ExpID)

clin.desc.gen$obese <- factor(clin.desc.gen$obese, levels = c("NonObese", "Overweight", "Obese"))

ggplot(clin.desc.gen %>% pivot_longer(cols = 3:5, names_to = "ClinParam", values_to = "Value"), aes(x=ClinParam, y=Value)) + 
  geom_violin() + facet_wrap(~ ClinParam, scales = "free") + theme_minimal()

ggplot(clin.desc.gen %>% pivot_longer(cols = 3:5, names_to = "ClinParam", values_to = "Value"), aes(x=obese, y=Value, fill=obese)) + 
  geom_violin() + facet_wrap(~ ClinParam, scales = "free") + theme_light()


sex <- table(clin.desc.gen$sex) %>% as.data.frame() %>% mutate(Perc = Freq/sum(Freq)*100, ypos = cumsum(Perc)- 0.5*Perc)
obese3 <- table(clin.desc.gen$obese3) %>% as.data.frame() %>% mutate(Perc = Freq/sum(Freq)*100, ypos = cumsum(Perc)- 0.5*Perc)

ggplot(obese3, aes(x="", y=Perc, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = ypos, label = Perc), color = "white", size=6)+
  theme_void()


##### Clinical correlations #####

colnames(pr.ln.norm)

pr.clin.cor <- pr.ln.norm[! pr.ln.norm$ExpID %in% samples.to.remove, c(1,2,11)] %>% nest(data = -ProteinGroups) 

clin.cor.dt <- clin.portal[clin.portal$ExpID %in% pr.ln.norm$ExpID, c(1,7:23)] %>% 
  mutate(basal_TG = basal/TG, iso_basal = iso/basal, iso_TG = iso/TG) %>%
  pivot_longer(cols = 2:21, names_to = "ClinParam", values_to = "Value") %>%
  #group_by(ClinParam) %>% mutate(Value = (Value - median(Value, na.rm = T))/sd(Value, na.rm = T)) %>% ungroup() %>%
  pivot_wider(names_from = ClinParam, values_from = Value)

#clin.cor.dt <- clin.portal[clin.portal$ExpID %in% pr.ln.norm$ExpID, c(1,7:23)]

clin.cor.res.s <- data.frame()
clin.cor.res.p <- data.frame()

chunk_size <- 100
num_chunks <- ceiling(nrow(pr.clin.cor) / chunk_size)
chunk_indices <- split(1:nrow(pr.clin.cor), cut(1:nrow(pr.clin.cor), num_chunks, labels = FALSE))

for (i in 1:num_chunks) { #nrow(pr.clin.cor)
  
  start_row <- (i - 1) * chunk_size + 1
  end_row <- min(i * chunk_size, nrow(pr.clin.cor))
  
  start_time <- Sys.time()
  
  dt <- pr.clin.cor[start_row:end_row,] #pr.clin.cor[i,2]
  dt <- dt %>% unnest(data) %>% left_join(clin.cor.dt) %>% pivot_wider(names_from = ProteinGroups, values_from = norm.vsn)
  
  corel <- correlation(dt[,c(2:21)], dt[,22:ncol(dt)], method = "spearman") # 18, 19
  corel <- as.data.frame(corel)
  #corel$Parameter2 <- pr.clin.cor[i,1]$ProteinGroups
  
  corelp <- correlation(dt[,c(2:21)], dt[,22:ncol(dt)], method = "pearson")
  corelp <- as.data.frame(corelp)
  #corelp$Parameter2 <- pr.clin.cor[i,1]$ProteinGroups
  
  clin.cor.res.s <- rbind(clin.cor.res.s, corel)
  clin.cor.res.p <- rbind(clin.cor.res.p, corelp)
  
  end_time <- Sys.time()
  processing_time <- end_time - start_time
  
  print(paste("Processed chunk", i, "of", num_chunks, "- Rows left:", nrow(pr.clin.cor) - max(unlist(chunk_indices[1:i]))))
  print(paste("Total processing time:", processing_time))
  
  gc()
}

colnames(clin.cor.res.s)[1:2] <- c("ClinParam", "ProteinGroups")
colnames(clin.cor.res.p)[1:2] <- c("ClinParam", "ProteinGroups")

clin.cor.res.s <- clin.cor.res.s %>% left_join(annot[,1:2])
clin.cor.res.p <- clin.cor.res.p %>% left_join(annot[,1:2])

#write_rds(clin.cor.res.s, "clin.cor.res.s.full.vsn.2.outliers.rds")
#write_rds(clin.cor.res.p, "clin.cor.res.p.full.vsn.2.outliers.rds")

# <- readRDS("clin.cor.res.s.rds")
#clin.cor.res.p <- readRDS("clin.cor.res.p.rds")

clin.cor.res.s$stars <- cut(clin.cor.res.s$p, breaks=c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), label=c("****", "***", "**", "*", "")) 
clin.cor.res.p$stars <- cut(clin.cor.res.p$p, breaks=c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), label=c("****", "***", "**", "*", "")) 

## Spearman

clin.cor.res.s %>% filter(Genes == "PCCA") %>% 
  ggplot(aes(x=reorder(ClinParam, desc(rho)), y=rho, fill = -log(p,10))) + geom_col() +
  geom_text(aes(label=stars), color="black", size=5) +
  scale_fill_distiller(palette = "PRGn") +
  theme_linedraw()

clin.cor.res.s.sig <- clin.cor.res.s %>% filter(p <= 0.05)
length(unique(clin.cor.res.s.sig$ProteinGroups))

ggplot(clin.cor.res.s[clin.cor.res.s$ProteinGroups %in% clin.cor.res.s.sig$ProteinGroups,], aes(x=ClinParam, y=Genes, fill=rho)) + geom_tile() + 
  scale_fill_distiller(palette = "PRGn") +
  geom_text(aes(label=stars), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="r") +
  theme_classic() + theme(axis.text.x=element_text(angle = -90, hjust = 0))

clin.cor.spear.bmi <- clin.cor.res.s.sig %>% filter(ClinParam == "BMI_0") %>% filter(abs(rho) >= 0.3)
clin.cor.spear.homa <- clin.cor.res.s.sig %>% filter(ClinParam == "HOMA_IR") %>% filter(abs(rho) >= 0.3)


#Screen

screen$stars <- cut(screen$p, breaks=c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), label=c("****", "***", "**", "*", "")) 
screen.p$stars <- cut(screen.p$p, breaks=c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), label=c("****", "***", "**", "*", "")) 

screen$Omic <- "Transcriptomics"
screen.p$Omic <- "Proteomics"

screen.all <- rbind(screen, screen.p)
screen.all %>% 
  ggplot(aes(x=ClinParam, y=Genes, fill=rho)) + geom_tile() + 
  scale_fill_distiller(palette = "PRGn") +
  geom_text(aes(label=stars), color="black", size=5) +
  labs(y=NULL, x=NULL, fill="rho") +
  theme_classic() + theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  facet_wrap(~ Omic)


screen.cor <- clin.cor.res.s %>% filter(Genes %in% screen$Genes & ClinParam %in% c("BMI", "iso", "HOMA", "basal"))
screen.cor$Omic <- "Proteomics"
ggplot(screen.cor, aes(x=ClinParam, y=Genes, fill=rho)) + geom_tile() + 
  scale_fill_distiller(palette = "PRGn") +
  geom_text(aes(label=stars), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="rho") +
  theme_classic() + theme(axis.text.x=element_text(angle = -90, hjust = 0))

screen$Omic <- "Transcriptomics"

screen %>% pivot_longer(cols = 2:6, names_to = "ClinParam", values_to = "rho") %>% 
  filter(ClinParam != "isobasal") %>% na.omit() %>%
  ggplot(aes(x=ClinParam, y=Genes, fill=rho)) + geom_tile() + 
  scale_fill_distiller(palette = "PRGn") +
  #geom_text(aes(label=stars), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="rho") +
  theme_classic() + theme(axis.text.x=element_text(angle = -90, hjust = 0))

screen.all <- rbind(screen %>% pivot_longer(cols = 2:6, names_to = "ClinParam", values_to = "rho"), screen.cor[,c(11,13,1,3)])
screen.all %>% filter(ClinParam != "isobasal") %>% na.omit() %>%
  ggplot(aes(x=ClinParam, y=Genes, fill=rho)) + geom_tile() + 
  scale_fill_distiller(palette = "PRGn") +
  labs(y=NULL, x=NULL, fill="rho") +
  theme_classic() + theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  facet_wrap(~ Omic)


## Pearson

clin.cor.res.p %>% filter(Genes == "CALR") %>% 
  ggplot(aes(x=reorder(ClinParam, desc(r)), y=r, fill = -log(p,10))) + geom_col() +
  geom_text(aes(label=stars), color="black", size=5) + ylim(-1,1) +
  scale_fill_distiller(palette = "PRGn") +
  theme_linedraw()

clin.cor.res.p.sig <- clin.cor.res.p %>% filter(p <= 0.05)
length(unique(clin.cor.res.p.sig$ProteinGroups))

ggplot(clin.cor.res.p[clin.cor.res.p$ProteinGroups %in% clin.cor.res.p.sig$ProteinGroups,], aes(x=ClinParam, y=Genes, fill=r)) + geom_tile() + 
  scale_fill_distiller(palette = "PRGn") +
  geom_text(aes(label=stars), color="black", size=5) + ylim(-1,1) + 
  labs(y=NULL, x=NULL, fill="r") +
  theme_classic() + theme(axis.text.x=element_text(angle = -90, hjust = 0))


##### NO Vs Ob #####

pr.ln.norm.2 <- pr.ln.norm[! pr.ln.norm.2$ExpID %in% samples.to.remove,] %>% 
  mutate(obese = case_when(obese == "lean" ~ "Non_Obese",
                           obese == "over weight" ~ "Non_Obese",
                           obese == "obese" ~ "Obese")) 
pr.ln.norm.no <- pr.ln.norm.2 %>% 
  filter(obese == "Non_Obese") %>% group_by(ProteinGroups) %>%
  summarise(lean.med = median(norm.vsn, na.rm = T))
pr.ln.norm.fc.2 <- pr.ln.norm.2 %>% left_join(pr.ln.norm.no) %>% mutate(log2fc = norm.vsn - lean.med)

t.test.obesity.2 <- pr.ln.norm.fc.2  %>%
  group_by(ProteinGroups) %>%
  pairwise_t_test(log2fc ~ obese, p.adjust.method = "fdr") 
t.test.obesity.2

t.test.obesity.2 <- t.test.obesity.2 %>% left_join(annot[,1:2]) 
length(unique(t.test.obesity.2$ProteinGroups[t.test.obesity.2$p.adj <= 0.05]))


wilc.obesity <- pr.ln.norm.2 %>%
  group_by(ProteinGroups) %>%
  pairwise_wilcox_test(norm.vsn ~ obese, p.adjust.method = "fdr") 
wilc.obesity

wilc.obesity <- wilc.obesity %>% left_join(annot[,1:2]) 
length(unique(wilc.obesity$ProteinGroups[wilc.obesity$p.adj <= 0.05]))

##### ANOVA ####

pr.ln.norm.sex.sep <- pr.ln.norm.2[! pr.ln.norm.2$ExpID %in% samples.to.remove, c(1,2,7,8,11)]

pr.ln.norm.sex.sep.sum <- pr.ln.norm.2[! pr.ln.norm.2$ExpID %in% samples.to.remove, c(1,2,7,8,11)] %>%
  group_by(ProteinGroups, sex, obese) %>%
  summarise(norm.vsn = median(norm.vsn, na.rm=T))

pr.ln.norm.sex.sep.lean <- pr.ln.norm.sex.sep %>% filter(obese == "Non_Obese") %>%
  group_by(ProteinGroups) %>%
  summarise(median.lean = median(norm.vsn, na.rm = T))
pr.ln.norm.sex.sep <- pr.ln.norm.sex.sep %>% left_join(pr.ln.norm.sex.sep.lean[,-2]) %>%
  mutate(log2fc = norm.vsn - median.lean) %>% drop_columns(c(6)) %>% ungroup()

anova.sex.sep <- pr.ln.norm.sex.sep %>% nest(data = -ProteinGroups)


anova.sex.sep.res.logfc.tukey <- list()

for (i in 1:nrow(anova.sex.sep)) {
  dt <- anova.sex.sep$data[[i]] %>% unnest(cols = c(log2fc, sex, obese))
  
  tryCatch({
    anova <- aov(log2fc ~ obese * sex, data = dt)
    anova_tidy <- tidy(anova)
    
    tukey_hsd <- TukeyHSD(anova)
    tukey_hsd_tidy <- tidy(tukey_hsd)
    
    
    result <- list(anova_tidy = anova_tidy,
                   posthoc_tidy = tukey_hsd_tidy)
    
    anova.sex.sep.res.logfc.tukey[[i]] <- result
  }, error = function(e) {
    # Print an error message if ANOVA or post hoc fails
    message("Error in processing ProteinGroups ", i, ": ", e$message)
    
    # Store NA if an error occurs
    anova.sex.sep.res.logfc.tukey[[i]] <- list(anova_tidy = NA, posthoc_tidy = NA)
  })
}

names(anova.sex.sep.res.logfc.tukey) <- anova.sex.sep$ProteinGroups

anova_results <- anova.sex.sep.res.logfc %>%
  map_dfr(pluck, "anova_tidy", .id = "ProteinGroups") %>%
  group_by(term) %>% mutate(padj = p.adjust(p.value, method = "fdr")) %>%
  left_join(annot[,1:2]) %>%
  filter(term != "Residuals")

anova_results_filtered <- anova_results %>% filter(padj <= 0.05)

posthoc_results_ob <- anova.sex.sep.res.logfc %>% 
  map_dfr(pluck, "posthoc_tidy_obese", .id = "ProteinGroups") %>%
  left_join(annot[,1:2])

posthoc_results_sex <- anova.sex.sep.res.logfc %>% 
  map_dfr(pluck, "posthoc_tidy_sex", .id = "ProteinGroups") %>%
  left_join(annot[,1:2])

posthoc_results_tukey <- anova.sex.sep.res.logfc.tukey %>% 
  map_dfr(pluck, "posthoc_tidy", .id = "ProteinGroups") %>%
  left_join(annot[,1:2])

proteins_obese <- posthoc_results_tukey %>% filter(adj.p.value <= 0.05 & term == "obese") 
length(unique(proteins_obese$ProteinGroups))
proteins_sex <- posthoc_results_tukey %>% filter(adj.p.value <= 0.05 & term == "sex") 
length(unique(proteins_sex$ProteinGroups))
proteins_interaction <- posthoc_results_tukey %>% filter(adj.p.value <= 0.05 & term == "obese:sex") 
length(unique(proteins_interaction$ProteinGroups))

pr.ln.norm.sex.sep.sum.obese <- pr.ln.norm.sex.sep %>% filter(obese != "Non_Obese") %>%
  group_by(ProteinGroups) %>% summarise(log2fc = median(log2fc, na.rm = T))

#anova.to.join <- posthoc_results_ob[,c(1,9,10)] 
anova.to.join <- posthoc_results_tukey[,c(1,2,8,9)] %>% filter(term == "obese")

pr.ln.norm.sex.sep.sum.obese <- left_join(pr.ln.norm.sex.sep.sum.obese, anova.to.join) %>%  
  mutate(levels = case_when(log2fc >= 0.3 & adj.p.value <= 0.05 ~ "Up-Regulated",
                            log2fc <= -0.3 & adj.p.value <= 0.05 ~ "Down-Regulated",
                            TRUE ~ "Unchanged"),
         logp = log(adj.p.value, 10)) %>%
  left_join(annot[,1:2])

pr.ln.norm.sex.sep.sum.obese$logp[pr.ln.norm.sex.sep.sum.obese$logp == -Inf] <- -15

top <- 20
top_genes <- bind_rows(pr.ln.norm.sex.sep.sum.obese  %>% filter(levels == 'Up-Regulated') %>%  arrange(logp, desc(abs(log2fc))) %>% head(top),
                       pr.ln.norm.sex.sep.sum.obese  %>% filter(levels == 'Down-Regulated') %>% arrange(logp, desc(abs(log2fc))) %>% head(top))
top_genes

options(ggrepel.max.overlaps = Inf)
pr.ln.norm.sex.sep.sum.obese %>% 
  ggplot(aes(x=log2fc, y=-logp)) + geom_point(aes(color = levels)) + ylim(0,15) + xlim(-2,2) +
  scale_color_manual(values = c("#004b6e", "#ebebeb", "#E2C109")) + #ebebeb E2C109 004b6e
  geom_label_repel(data = top_genes,
                   mapping = aes(log2fc, -logp, label = Genes)) +
  xlab("log2FC") + ylab("-log10(adj.p-value)") +
  theme_classic()



### Post Hoc

posthoc_results <- anova.sex.sep.res.logfc.tukey %>% 
  map_dfr(pluck, "posthoc_tidy", .id = "ProteinGroups") %>%
  left_join(annot[,1:2])

posthoc_results_filtered <- posthoc_results %>% filter(adj.p.value <= 0.05 & term == "obese:sex")
length(unique(posthoc_results_filtered$ProteinGroups))
unique(posthoc_results_filtered$contrast)
#posthoc_results_filtered <- posthoc_results_filtered %>% filter(contrast %in% c(  "Obese:male-Obese:female",  "Non_Obese:male-Non_Obese:female", "Non_Obese:male-Obese:female", "Obese:male-Non_Obese:female"))
length(unique(posthoc_results_filtered$ProteinGroups))

top.genes.list.sex.sep.anova <- list(
  NOF_NOM = posthoc_results_filtered$ProteinGroups[posthoc_results_filtered$contrast == "Non_Obese:male-Non_Obese:female"],
  NOF_OBM = posthoc_results_filtered$ProteinGroups[posthoc_results_filtered$contrast == "Obese:male-Non_Obese:female"],
  OBF_NOM = posthoc_results_filtered$ProteinGroups[posthoc_results_filtered$contrast == "Non_Obese:male-Obese:female"],
  OBF_OBM = posthoc_results_filtered$ProteinGroups[posthoc_results_filtered$contrast == "Obese:male-Obese:female"],
  NOF_OBF = posthoc_results_filtered$ProteinGroups[posthoc_results_filtered$contrast == "Obese:female-Non_Obese:female"],
  NOM_OBM = posthoc_results_filtered$ProteinGroups[posthoc_results_filtered$contrast == "Obese:male-Non_Obese:male"])


UpSetR::upset(fromList(top.genes.list.sex.sep.anova),
              order.by = "freq", decreasing = T, nintersects = NA, nsets = 6,
              #scale.intersections = "log10", scale.sets = "log10",
              matrix.color="black", point.size=5) 


pr.ln.norm.sex.sep.sum.sex.infl <- pr.ln.norm.sex.sep[pr.ln.norm.sex.sep$ProteinGroups %in% posthoc_results_filtered$ProteinGroups,] %>% 
  left_join(annot[,1:2]) %>%
  group_by(ProteinGroups, sex, obese) %>%
  summarise(log2fc = median(log2fc, na.rm=T),
            norm.vsn = median(norm.vsn, na.rm = T)) %>% ungroup()
pr.ln.norm.sex.sep.sum.sex.infl$merged <- paste(pr.ln.norm.sex.sep.sum.sex.infl$obese, pr.ln.norm.sex.sep.sum.sex.infl$sex, sep = " - ")

pr.ln.norm.sex.sep.sum.sex.infl.wd <- pr.ln.norm.sex.sep.sum.sex.infl[,-c(2:3,5)] %>% filter(merged %in% c("Obese - female", "Obese - male")) %>%
  filter(ProteinGroups %in%  posthoc_results_filtered$ProteinGroups[posthoc_results_filtered$contrast == "Obese:male-Obese:female"]) %>%
  pivot_wider(names_from = merged, values_from = log2fc) %>% left_join(annot[,1:2]) %>% as.data.frame()
rownames(pr.ln.norm.sex.sep.sum.sex.infl.wd) <- pr.ln.norm.sex.sep.sum.sex.infl.wd$ProteinGroups
pheatmap(pr.ln.norm.sex.sep.sum.sex.infl.wd[,-c(1,4)], 
         color = colorRampPalette(c(portalcol))(20), 
         clustering_method = "ward.D2",
         show_rownames = F, show_colnames = T, breaks = breaksList,
         cutree_rows = 5)
length(pr.ln.norm.sex.sep.sum.sex.infl.wd$ProteinGroups)

obm_vs_obf <- pr.ln.norm.sex.sep.sum.sex.infl[,-c(2:3,5)] %>% filter(merged %in% c("Obese - female", "Obese - male")) %>%
  filter(ProteinGroups %in%  posthoc_results_filtered$ProteinGroups[posthoc_results_filtered$contrast == "Obese:male-Obese:female"])

# Clustering

dist.markers.ob <- dist(pr.ln.norm.sex.sep.sum.sex.infl.wd[,-c(1,4)], method = 'euclidean')

set.seed(123)
hclust.markers.ob <- hclust(dist.markers.ob, method = 'ward.D2')
plot(hclust.markers.ob, labels = F)
rect.hclust(hclust.markers.ob, h = 5, border = 2:6)

hcl.markers.grp.ob <- cutree(hclust.markers.ob, h = 5)
table(hcl.markers.grp.ob)
hcl.markers.grp.ob <- cutree(hclust.markers.ob, h = 5) %>% as.data.frame() %>% rownames_to_column("ProteinGroups") %>%
  left_join(pr.ln.norm.sex.sep.sum.sex.infl.wd)
colnames(hcl.markers.grp.ob)[2] <- c("Cluster")
rownames(hcl.markers.grp.ob) <- hcl.markers.grp.ob$ProteinGroups

breaksList = seq(-1, 1, by = 0.1)

cl1.ob <- pr.ln.norm.sex.sep.sum.sex.infl.wd[pr.ln.norm.sex.sep.sum.sex.infl.wd$ProteinGroups %in% hcl.markers.grp.ob$ProteinGroups[hcl.markers.grp.ob$Cluster == 1],]
pheatmap(as.matrix(cl1.ob[,-c(1,4)]), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl2.ob <- pr.ln.norm.sex.sep.sum.sex.infl.wd[pr.ln.norm.sex.sep.sum.sex.infl.wd$ProteinGroups %in% hcl.markers.grp.ob$ProteinGroups[hcl.markers.grp.ob$Cluster == 2],]
pheatmap(as.matrix(cl2.ob[,-c(1,4)]), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl3.ob <- pr.ln.norm.sex.sep.sum.sex.infl.wd[pr.ln.norm.sex.sep.sum.sex.infl.wd$ProteinGroups %in% hcl.markers.grp.ob$ProteinGroups[hcl.markers.grp.ob$Cluster == 3],]
pheatmap(as.matrix(cl3.ob[,-c(1,4)]), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl4.ob <- pr.ln.norm.sex.sep.sum.sex.infl.wd[pr.ln.norm.sex.sep.sum.sex.infl.wd$ProteinGroups %in% hcl.markers.grp.ob$ProteinGroups[hcl.markers.grp.ob$Cluster == 4],]
pheatmap(as.matrix(cl4.ob[,-c(1,4)]), show_rownames = F, show_colnames = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))
cl5.ob <- pr.ln.norm.sex.sep.sum.sex.infl.wd[pr.ln.norm.sex.sep.sum.sex.infl.wd$ProteinGroups %in% hcl.markers.grp.ob$ProteinGroups[hcl.markers.grp.ob$Cluster == 5],]
pheatmap(as.matrix(cl5.ob[,-c(1,4)]), show_rownames = T, show_colnames = F, cluster_rows = F, breaks = breaksList, color = colorRampPalette(c(portalcol))(20))


gene1 <- cl1.ob$Genes
id1 <- mapIds(org.Hs.eg.db, gene1, 'ENTREZID', 'SYMBOL')
gene2 <- cl2.ob$Genes
id2 <- mapIds(org.Hs.eg.db, gene2, 'ENTREZID', 'SYMBOL')
gene3 <- cl3.ob$Genes
id3 <- mapIds(org.Hs.eg.db, gene3, 'ENTREZID', 'SYMBOL')
gene4 <- cl4.ob$Genes
id4 <- mapIds(org.Hs.eg.db, gene4, 'ENTREZID', 'SYMBOL')

ego1 <- enrichGO(gene          = id1,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05)

ego2 <- enrichGO(gene          = id2,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.06)

ego3 <- enrichGO(gene          = id3,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.06)

ego4 <- enrichGO(gene          = id4,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.06)

cowplot:: plot_grid(dotplot(ego1), dotplot(ego2), dotplot(ego3), dotplot(ego4))

ego1res <- ego1@result %>% filter(p.adjust <= 0.06) %>% mutate(cluster = "1") %>% arrange(p.adjust) %>%  slice_head(n = 5)  
ego2res <- ego2@result %>% filter(p.adjust <= 0.06) %>% mutate(cluster = "2") %>% arrange(p.adjust) %>%  slice_head(n = 5)  
ego3res <- ego3@result %>% filter(p.adjust <= 0.06) %>% mutate(cluster = "3") %>% arrange(p.adjust) %>%  slice_head(n = 5)  
ego4res <- ego4@result %>% filter(p.adjust <= 0.06) %>% mutate(cluster = "4") %>% arrange(p.adjust) %>%  slice_head(n = 5)  
egores <- rbind(ego1res, ego2res, ego3res, ego4res)

ggplot(egores, aes(x=cluster, y=Description, size=-log(p.adjust,10), colour=p.adjust)) + geom_point() + 
  theme_light() + scale_colour_gradientn(colours = portalcol) +
  facet_wrap(~cluster, scales = "free")




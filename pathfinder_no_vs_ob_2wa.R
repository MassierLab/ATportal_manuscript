##### Libraries ####

Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_202\\')

library(pathfindR)
library(tidyverse)
library(msigdbr)
library(viridis)

portalcol <- colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(51)
portalcol2 <- c("#1a4659", "#e2c744", "#ba5c28", "#4f736a", "#f0cbbf", "#b2dfda", "#f59b7c", "#7dbfb3", "#939a64", "#8ca5a5", "#adbec4", "#c3b6b6")

##### Data #####

getwd()
setwd("~/")
setwd("WAT Portal/WAT_panel/Results")

#load(("no_vs_ob_pathfinder.RData"))

pr.ln.norm.fc <- readRDS("pr.ln.norm.sex.sep.sum.obese.tukey.RDS")

genes <- pr.ln.norm.fc %>% separate_rows(Genes, sep=";")
genes <- genes$Genes


####

pf.input.ob <- pr.ln.norm.fc[,c(5,2,4)] %>%  ungroup() %>%
  filter(adj.p.value <= 0.05) %>% separate_rows(Genes, sep=";")
colnames(pf.input.ob) <- c("Gene.symbol", "logFC", "adj.P.Val")

gsets_list <- get_gene_sets_list(
  source = "MSigDB",
  collection = "H")
hallmark_gsets <- gsets_list[[1]]
hallmark_descriptions <- gsets_list[[2]]

##### Pathway Analysis Obese #####

###### Manual

processed.ob <- input_processing(
  input = pf.input.ob, # the input: in this case, differential expression results
  p_val_threshold = 0.05, # p value threshold to filter significant genes
  pin_name_path = "IntAct", # the name of the PIN to use for active subnetwork search
  convert2alias = TRUE # boolean indicating whether or not to convert missing symbols to alias symbols in the PIN
)

duplicated_genes <- pf.input.ob$Gene.symbol[duplicated(pf.input.ob$Gene.symbol)]

n_iter <- 10 ## number of iterations
combined_res_ob <- NULL ## to store the result of each iteration
current_res <- NULL

for (i in 1:n_iter) {
  ###### Active Subnetwork Search
  snws_file <- paste0("active_snws_ob_vs_no_2wa", i) # Name of output file
  active_snws <- active_snw_search(
    input_for_search = processed.ob,
    pin_name_path = "IntAct",
    snws_file = snws_file,
    score_quan_thr = 0.8, # you may tweak these arguments for optimal filtering of subnetworks
    sig_gene_thr = 0.05, # you may tweak these arguments for optimal filtering of subnetworks
    search_method = "GR", # we suggest using GR
    seedForRandom = i # setting seed to ensure reproducibility per iteration
  )
  
  ###### Enrichment Analyses
  current_res <- enrichment_analyses(
    snws = active_snws,
    sig_genes_vec = processed.ob$GENE,
    pin_name_path = "IntAct",
    genes_by_term = hallmark_gsets,
    term_descriptions = hallmark_descriptions,
    adj_method = "fdr",
    enrichment_threshold = 0.05,
    list_active_snw_genes = TRUE
  ) # listing the non-input active snw genes in output
  
  ###### Combine results via `rbind`
  combined_res_ob <- rbind(combined_res_ob, current_res)
}

###### Summarize Combined Enrichment Results
summarized_df_ob <- summarize_enrichment_results(combined_res_ob,
                                                 list_active_snw_genes = TRUE)

###### Annotate Affected Genes Involved in Each Enriched Term
final_res_ob <- annotate_term_genes(
  result_df = summarized_df_ob,
  input_processed = processed.ob,
  genes_by_term = hallmark_gsets)

enrichment_chart(result_df = final_res_ob,top_terms = 15) + scale_colour_gradientn(colours = portalcol)

final_res_ob$Term_Description <- factor(final_res_ob$Term_Description )

top15 <- final_res_ob %>% arrange(lowest_p) %>% slice_head(n = 15) %>% select(Term_Description)    

final_res_ob_15 <- final_res_ob[final_res_ob$Term_Description %in% top15$Term_Description,] %>%
  mutate(logp = -log10(lowest_p))

ggplot(final_res_ob[final_res_ob$Term_Description %in% top15$Term_Description,], #%>% filter(lowest_p <= 0.05 & occurrence > 5), 
       aes(x=Fold_Enrichment, y=reorder(Term_Description, desc(lowest_p)), size=occurrence, colour=lowest_p)) + geom_point() +
  scale_colour_gradientn(colours = portalcol) +
  theme_light()

genes_in_pathways_ob <- final_res_ob[final_res_ob$Term_Description %in% top15$Term_Description,c(2,3,6,9,10)] %>% 
  pivot_longer(cols = 4:5, names_to = "Regulation", values_to = "Genes") %>%
  separate_rows(Genes, sep = ", ") %>% mutate(Comparison = "OB_Vs_NO") %>%
  group_by(Term_Description) %>% mutate(NoGenes = nrow(Genes))

count <- genes_in_pathways_ob %>%
  group_by(Term_Description, Regulation) %>%
  summarise(Gene_Count = n(), .groups = 'drop') %>%
  group_by(Term_Description) %>% mutate(Genes_total = sum(Gene_Count)) %>% ungroup() %>%
  mutate(Perc = Gene_Count/Genes_total*100)

final_res_ob_15 <- final_res_ob_15 %>% left_join(count[,c(1,4)])

ggplot(final_res_ob_15, 
       aes(x=Fold_Enrichment, y=reorder(Term_Description, logp), size=Genes_total, colour=logp)) + geom_point() +
  scale_colour_gradientn(colours = c("#F2EBDD", "#E2C744")) +
  theme_light()

count$Term_Description <- factor(count$Term_Description, levels = top15$Term_Description)
ggplot(count, aes(x=Regulation, y=Term_Description, fill = Perc)) + geom_tile() +
  scale_fill_gradientn(colours = c("#FFFFFF", "#1A4659"), limits=c(0,100)) +
  theme_minimal()

#cluster_enriched_terms(final_res_ob, method = "fuzzy", use_description = T, use_active_snw_genes = T)
#term_gene_graph(final_res_ob, use_description = TRUE)
#term_gene_graph(final_res_ob, num_terms = 5, use_description = TRUE, node_size = "p_val")
#gg_list_ob <- visualize_term_interactions(final_res_ob, pin_name_path = 'IntAct', show_legend = TRUE)

#save.image("no_vs_ob_pathfinder_2wa.RData")

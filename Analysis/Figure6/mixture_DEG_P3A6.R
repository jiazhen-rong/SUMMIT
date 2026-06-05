library(spotglm)
library(ggplot2) # for plotting
library(dplyr) # for manipulating data tables
library(patchwork) # for qqplots side-by-side
library(Matrix)
library(Seurat)
library(grid)
library(gridExtra)
library(spacexr)
library(CARD)
library(dplyr)
library(ggplot2)
library(SummarizedExperiment)
library(Matrix)
library(ape)
library(distances)
library(sf)
library(spdep)
library(cowplot)
library(SpatialMT)
library(stringr)

#data =sSeurat#data =spotglm::read_example_visium_colorectal_cancer_data()
setwd("~/nzhanglab/project/jrong/mito_LT/scripts/our_model/")
source("example_data/210215_FunctionsGeneral.R") # from MAESTER paper

slide_name = "P3A6"

# load transcriptomics/seeker data
seu <-readRDS("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/P3A6/seurat_final.rds")
mt.genes <- grep("^MT-", rownames(seu), value = TRUE)
filtered.genes <- rownames(seu)[!rownames(seu) %in% mt.genes]
seu <- subset(seu,features=filtered.genes)
coords = seu@reductions$SPATIAL@cell.embeddings
# RCTD ratio
rctd_ratio_major = readRDS("../../results/Sydney_Bracht/sample_level/P3A6/RCTD/Dale_P03B2_HCC/RCTD_res.rds")
# normalize the celltype weights to sum to 1 in each spot
rctd_ratio_major = as.data.frame(as.matrix(normalize_weights(rctd_ratio_major@results$weights))) 
rownames(rctd_ratio_major) = paste0(rownames(rctd_ratio_major),"-1")
#Ws = rctd_ratio_major[colnames(af.dm),]
#celltypes= colnames(Ws)
deconv=rctd_ratio_major 
#deconv=as.matrix(deconv)
# create niche by VOI:
maegatk.rse = readRDS("../../data/Sydney_Bracht/P3A6/maegatk_final.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))#*100 # all possible mitochodnrial mutations' (4* 16K) x spot' VAF
# prepare coverage N, # spot x each chrM location's (16K) coverage
counts=as.matrix(maegatk.rse @assays@data$coverage)
rownames(counts) = 1:nrow(counts);
colnames(counts) = maegatk.rse @colData@rownames
N=as(as.matrix(counts[sapply(strsplit(rownames(af.dm),"_"),"[[",1),]), "sparseMatrix")
rownames(N) <- rownames(af.dm)
# load variant of interest
voi = read.table(paste0("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/P3A6/P3A6_voi_be25vaf.tsv"),sep="\t")[,1]
voi = c("15723_G>A")
niche = as.matrix(Matrix(0,nrow=dim(seu)[2],ncol=2)) # X covariate matrix
rownames(niche) = colnames(seu);colnames(niche)=c(voi ,"Rest")
niche[colnames(af.dm)[af.dm["15723_G>A",]>0.25],"15723_G>A"]=1
#niche[intersect(colnames(af.dm)[af.dm["11456_G>A",]>0.25],rownames(niche)),"11456_G>A"]=1
niche[(niche[,"15723_G>A"]==0)  ,
        "Rest"] = 1

data <- list()
inter_bc = intersect(colnames(seu),rownames(deconv))
data$coords = coords[inter_bc ,]
data$deconv = deconv[inter_bc ,]
# Subset on the 1000 most variable genes
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
data$counts = t(seu@assays$RNA$counts)[inter_bc ,VariableFeatures(seu)]
data$library_size = rowSums(data$counts)[inter_bc ]
data$niche=niche[inter_bc,]
data$niche <- as.matrix(data$niche)
#data$deconv[is.na(data$deconv)]=0.0001
#rm(af.dm)
rm(maegatk.rse)
gc()

# Fitting SpotGLM
num_genes = ncol(data$counts)
res = vector("list",num_genes) # container for the results

# grab the current grad
g_old <- spotglm:::spot_poisson[["grad"]]

# define patched version
g_new <- function(X, y, beta, lambda, offset) {
  means <- spotglm:::spot_poisson[["predict"]](X, beta, lambda, offset)
  CT_means <- means$individual
  total_means <- means$total
  weights <- -lambda * CT_means + y * (1 / total_means) * lambda * CT_means
  # Key change
  grad_beta <- t(as.matrix(X)) %*% as.matrix(weights)
  
  list(grad = grad_beta, weights = weights)
}

# replace 
assignInNamespace("spot_poisson",
                  within(spotglm:::spot_poisson, { grad <- g_new }),
                  ns = "spotglm")

make_zero_result <- function(X, lambda) {
  list(
    beta_estimate = matrix(
      0,
      nrow = ncol(X),
      ncol = ncol(lambda),
      dimnames = list(colnames(X), colnames(lambda))
    ),
    vcov = matrix(
      0,
      nrow = ncol(X) * ncol(lambda),
      ncol = ncol(X) * ncol(lambda)
    ),
    converged = FALSE,
    niter = NA,
    error = TRUE
  )
}

fit_one_gene <- function(j,data) {
  y <- as.numeric(data$counts[, j])
  
  # Optional: skip all-zero genes early
  if (sum(y, na.rm = TRUE) == 0) {
    return(make_zero_result(data$niche, data$deconv))
  }
  
  tryCatch(
    spotglm::run_model(y = data$counts[,j],
                       X = data$niche,
                       lambda = data$deconv,
                       family = "spot poisson",
                       offset = log(data$library_size),
                       initialization = T,batch_size = 250), 
   error = function(e) {
      message(j,": Gene ", colnames(data$counts)[j], " failed, setting result to 0")
      make_zero_result(data$niche, data$deconv)
    }
  )
}


t1 = Sys.time()
for(j in c(1:num_genes)){
  if(j%%100 == 0){
    cat("Fitting model for gene ",j," out of ", num_genes,"\n")
    print(Sys.time() - t1)
  }
  res[[j]] <- fit_one_gene(j,data)

  # res[[j]] = spotglm::run_model(y = data$counts[,j],
  #                               X = data$niche,
  #                               lambda = data$deconv,
  #                               family = "spot poisson",
  #                               offset = log(data$library_size),
  #                               initialization = T,batch_size = 250)
  
}

names(res) = colnames(data$counts)
saveRDS(res,"~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/P3A6/spotglm_res.rds")


res <- readRDS("~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/P3A6/spotglm_res.rds")


library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(pheatmap)

cell_types <- colnames(deconv)          # c("BE","SQ","GC","TC","IM","FB","VC","NS")
niches <- colnames(niche)               # 4 niches/covariates
pairs <- combn(niches, 2, simplify = FALSE)

all_contrasts <- purrr::map_dfr(cell_types, function(ct) {
  purrr::map_dfr(pairs, function(p) {
    n1 <- p[1]; n2 <- p[2]
    out <- compute_contrast_significance(
      input_list = res,
      cell_type = ct,
      effect_names = c(n1, n2),
      beta_name = "beta_estimate",
      covariance_name = "vcov",
      #sided = 1,direction = "pos",
      sided=2
    )
    
    # If gene names are rownames, promote them to a column
    if (!("gene" %in% names(out))) out <- out %>% tibble::rownames_to_column("gene")
    
    out %>%
      mutate(
        contrast = paste0(n1, " vs ", n2)
      )
  })
})
#saveRDS(all_contrasts,"~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/A200_S6/deg.rds")
saveRDS(all_contrasts,"~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/P3A6/deg_2sided.rds")

all_contrasts <- readRDS("~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/P3A6/deg_2sided.rds")
# remove MT genes
all_contrasts = all_contrasts[!grepl("MT-",all_contrasts$name),]
summary_by_ct_contrast <- all_contrasts %>%
  group_by(cell_type, contrast) %>%
  summarise(
    n_q05 = n_distinct(gene[qval < 0.05]),
    genes_q05 = list(sort(unique(gene[qval < 0.05]))),
    n_q10 = n_distinct(gene[qval < 0.10]),
    genes_q10 = list(sort(unique(gene[qval < 0.10]))),
    .groups = "drop"
  ) %>%
  arrange(cell_type, contrast)

#summary_by_ct_contrast
sig_gene_list <- all_contrasts %>%
  filter(qval < 0.1) %>%
  group_by(cell_type, contrast) %>%
  summarise(sig_genes = paste0(list(unique(gene))), n_sign_genes=length(unique(gene)),
            .groups = "drop")

# Plotting DEGs
source("~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/mixture_DEG_pathway_functions.R")
pdf_path="~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/P3A6/celltype_deg_top20_2sided_MTremoved.pdf"
plot_DEG(all_contrasts, cell_types=cell_types,
         pdf_path=pdf_path, effect_col="test_statistic",
         top_k=20, q_thresh=0.1,qcap_val=10)



# #### Pathway Enrichment Analysis
# library(dplyr)
# ct <- "HCC"
# contrast_name <- "3054_G>C vs 15777_G>C"
# 
# rank_df <- all_contrasts %>%
#   filter(cell_type == ct, contrast == contrast_name) %>%
#   filter(!is.na(pval)) %>%
#   mutate(
#     stat = test_statistic
#   ) %>%
#   select(gene, stat) %>%
#   distinct(gene, .keep_all = TRUE)
# 
# ranks <- rank_df$stat
# names(ranks) <- rank_df$gene
# ranks <- sort(ranks, decreasing = TRUE)
# 
# library(msigdbr)
# 
# hallmark <- msigdbr(
#   species = "Homo sapiens",
#   category = "H"
# ) %>%
#   split(x = .$gene_symbol, f = .$gs_name)
# 
# reactome <- msigdbr(
#   species = "Homo sapiens",
#   category = "C2",
#   subcategory = "CP:REACTOME"
# ) %>%
#   split(x = .$gene_symbol, f = .$gs_name)
# 
# go_bp <- msigdbr(
#   species = "Homo sapiens",
#   category = "C5",
#   subcategory = "GO:BP"
# )
# 
# go_mf <- msigdbr(
#   species = "Homo sapiens",
#   category = "C5",
#   subcategory = "GO:MF"
# )
# 
# go_cc <- msigdbr(
#   species = "Homo sapiens",
#   category = "C5",
#   subcategory = "GO:CC"
# )
# 
# go <- bind_rows(go_bp, go_mf, go_cc) %>%
#   split(x = .$gene_symbol, f = .$gs_name)
# 
# pathways_list <- list(
#   Hallmark = hallmark,
#   Reactome = reactome,
#   GO = go
# )
# 
# rank_df <- all_contrasts %>%
#   filter(cell_type == "BE", contrast == "3054_G>C vs 15777_G>C") %>%
#   filter(!is.na(pval))
# 
# # Preferred ranking statistic
# rank_df <- rank_df %>%
#   mutate(stat = estimate / se)   # or use stat column directly
# 
# ranks <- rank_df$stat
# names(ranks) <- rank_df$gene
# ranks <- sort(ranks, decreasing = TRUE)
# 
# ibrary(fgsea)
# 
# gsea_results <- purrr::imap_dfr(pathways_list, function(pw, pw_name) {
#   
#   fg <- fgsea(
#     pathways = pw,
#     stats = ranks,
#     nperm = 10000,
#     minSize = 15,
#     maxSize = 500
#   )
#   
#   fg %>%
#     mutate(
#       collection = pw_name,
#       contrast = "3054_G>C vs 15777_G>C",
#       cell_type = "BE"
#     )
# })
# 
# sig_pw <- gsea_results %>%
#   filter(padj < 0.05) %>%
#   arrange(collection, padj)
# 
# head(sig_pw, 20)
# 
# library(ggplot2)
# 
# top_pw <- sig_pw %>%
#   group_by(collection) %>%
#   slice_head(n = 10)
# 
# ggplot(top_pw,
#        aes(x = reorder(pathway, NES),
#            y = NES,
#            size = -log10(padj),
#            color = collection)) +
#   geom_point() +
#   coord_flip() +
#   theme_minimal() +
#   labs(x = NULL, title = "GSEA results")

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

# =====================================================================
# A200_s6 mixture-DEG analysis with spotGLM.
# Identical to mixture_DEG_path_A200_s6_sec1.R EXCEPT the GLM family is
# changed from "spot poisson" -> "spot negative binomial".
#
# Notes on the NB change:
#  - spotGLM selects the family by string: family = "spot negative binomial".
#  - NB estimates an overdispersion parameter internally; each gene's fit
#    in the result list gains a $dispersion element.
#  - The Poisson script monkey-patched spot_poisson[["grad"]]. That patch is
#    NOT needed here: (1) it is now identical to the package's built-in
#    Poisson gradient, and (2) the NB family ships its own dispersion-aware
#    gradient spot_negative_binomial[["grad"]]. So no assignInNamespace patch.
#  - Offset (log library size) is supported for NB just like Poisson.
# =====================================================================

setwd("~/nzhanglab/project/jrong/mito_LT/scripts/our_model/")
source("example_data/210215_FunctionsGeneral.R") # from MAESTER paper

slide_name = "A200_s6"

# load transcriptomics/seeker data
seu <- readRDS("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/A200_s6/MAESTER/MAESTER_subsets/a200_s6_final.rds")
coords = seu@reductions$SPATIAL@cell.embeddings
# RCTD ratio
rctd_ratio_major = readRDS("../../results/Sydney_Bracht/sample_level/a200_s6_high_cov/RCTD/Major/RCTD_res.rds")
# normalize the celltype weights to sum to 1 in each spot
rctd_ratio_major = as.data.frame(as.matrix(normalize_weights(rctd_ratio_major@results$weights)))
rownames(rctd_ratio_major) = paste0(rownames(rctd_ratio_major),"-1")
deconv=rctd_ratio_major
# create niche by VOI:
maegatk.rse = readRDS("../../data/Sydney_Bracht/A200_s6/MAESTER/MAESTER_subsets/maegatk_mr1.rds")
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse))# all possible mitochodnrial mutations' (4* 16K) x spot' VAF
# prepare coverage N, # spot x each chrM location's (16K) coverage
counts=as.matrix(maegatk.rse @assays@data$coverage)
rownames(counts) = 1:nrow(counts);
colnames(counts) = maegatk.rse @colData@rownames
N=as(as.matrix(counts[sapply(strsplit(rownames(af.dm),"_"),"[[",1),]), "sparseMatrix")
rownames(N) <- rownames(af.dm)
# load variant of interest
voi = read.table(paste0("~/nzhanglab/project/jrong/mito_LT/data/Sydney_Bracht/CRC076_C/MAESTER/lineage/CRC076_voi_vaf.tsv"),sep="\t")[,1]
voi = c("3054_G>C","15777_G>C","3071_T>C")
niche = as.matrix(Matrix(0,nrow=dim(seu)[2],ncol=4)) # X covariate matrix
rownames(niche) = colnames(seu);colnames(niche)=c("3054_G>C","15777_G>C","3071_T>C","Rest")
niche[colnames(af.dm)[af.dm["3054_G>C",]>0.25],"3054_G>C"]=1
niche[intersect(colnames(af.dm)[af.dm["15777_G>C",]>0.25],rownames(niche)),"15777_G>C"]=1
niche[intersect(colnames(af.dm)[af.dm["3071_T>C",]>0.25],rownames(niche)),"3071_T>C"]=1
niche[(niche[,"3054_G>C"]==0) & (niche[,"15777_G>C"]==0) &(niche[,"3071_T>C"]==0),
        "Rest"] = 1

data <- list()
inter_bc = intersect(colnames(seu),rownames(deconv))
data$coords = coords[inter_bc ,]
data$deconv = deconv[inter_bc ,]
# Subset on the 2000 most variable genes
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
data$counts = t(seu@assays$RNA$counts)[inter_bc ,VariableFeatures(seu)]
data$library_size = rowSums(data$counts)[inter_bc ]
data$niche=niche[inter_bc,]
data$niche <- as.matrix(data$niche)
rm(af.dm)
rm(maegatk.rse)
gc()

# Fitting SpotGLM (Negative Binomial)
num_genes = ncol(data$counts)
res = vector("list",num_genes) # container for the results

# NOTE: no spot_poisson gradient patch here. The NB family uses its own
# dispersion-aware gradient (spot_negative_binomial[["grad"]]).

t1 = Sys.time()
for(j in c(1:num_genes)){
  if(j%%100 == 0){
    cat("Fitting model for gene ",j," out of ", num_genes,"\n")
    print(Sys.time() - t1)
  }
  res[[j]] = spotglm::run_model(y = data$counts[,j],
                                X = data$niche,
                                lambda = data$deconv,
                                family = "spot negative binomial",
                                offset = log(data$library_size),
                                initialization = T,batch_size = 250)
}
names(res) = colnames(data$counts)
saveRDS(res,"~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/A200_S6/spotglm_res_nb.rds")

res <- readRDS("~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/A200_S6/spotglm_res_nb.rds")


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
saveRDS(all_contrasts,"~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/A200_S6/deg_2sided_nb.rds")

all_contrasts <- readRDS("~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/A200_S6/deg_2sided_nb.rds")
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

sig_gene_list <- all_contrasts %>%
  filter(qval < 0.1) %>%
  group_by(cell_type, contrast) %>%
  summarise(sig_genes = paste0(list(unique(gene))), n_sign_genes=length(unique(gene)),
            .groups = "drop")

# Plotting DEGs
source("~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/mixture_DEG_pathway_functions.R")
pdf_path="~/nzhanglab/project/jrong/mito_LT/scripts/SpotGLM/A200_S6/celltype_deg_top20_2sided_MTremoved_nb.pdf"
plot_DEG(all_contrasts, cell_types=cell_types,
         pdf_path=pdf_path, effect_col="test_statistic",
         top_k=20, q_thresh=0.1,qcap_val=10)

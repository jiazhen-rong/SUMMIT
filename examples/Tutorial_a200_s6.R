library(grid)
library(gridExtra)
library(spacexr)
library(CARD)
library(dplyr)
library(ggplot2)
library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(SUMMIT)

# load data (X - alternate allele counts, N - total counts, Ws - celltype)
setwd("~/Documents/GitHub/SpatialMT/examples/")
source("utility_prev_literature/210215_FunctionsGeneral.R") # from MAESTER paper
#source("../diagnostic_plot.R")

# load Spatial Transcriptomics/Curio Seeker data
seu <-readRDS("example_data/a200_s6_final.rds")
# obtain spatial coordinates
spatial_coords = seu@reductions$SPATIAL@cell.embeddings

# load MT data
maegatk.rse = readRDS("example_data//maegatk_mr1.rds")
# all possible mitochodnrial mutations' (4* 16K) x spot' VAF
af.dm <- data.matrix(computeAFMutMatrix(maegatk.rse)) # 0-1 VAF scale
# prepare coverage N, # spot x each chrM location's (16K) coverage
counts=as.matrix(maegatk.rse @assays@data$coverage)
rownames(counts) = 1:nrow(counts);
colnames(counts) = maegatk.rse @colData@rownames
N=as(as.matrix(counts[sapply(strsplit(rownames(af.dm),"_"),"[[",1),]), "sparseMatrix")
rownames(N) <- rownames(af.dm)

# load variant of interest
voi =  c("3054_G>C","15777_G>C","3071_T>C")

# subset to matrix to be much smaller sizes
af.dm=af.dm[voi,]
N = N[voi,]
N = N[,colnames(af.dm)]

# Load RCTD ratio
rctd_ratio_major = readRDS("example_data//RCTD_res.rds")
# normalize the celltype weights to sum to 1 in each spot
rctd_ratio_major = as.data.frame(as.matrix(normalize_weights(rctd_ratio_major@results$weights)))
rownames(rctd_ratio_major) = paste0(rownames(rctd_ratio_major),"-1")
Ws = rctd_ratio_major[colnames(af.dm),]
celltypes= colnames(Ws)

save_path = paste0("example_data/results")
dir.create(save_path)

# global celltype test
res_lg = celltype_test(celltypes = celltypes,voi=voi,N=as.matrix(N),
                       vaf=as.matrix(af.dm), Ws=Ws, spatial_coords = spatial_coords,
                       test_type = "linear",plot=T,
                       save_path=save_path,method="FDR",
                       figure_height = 15, figure_width = 15)

beta_list = power_analysis_all(voi=voi,celltypes = celltypes,Ws=Ws,
                               N=N,vaf=af.dm,X=NULL,
                               sample_num=100,alpha=0.05,n_sim=10,
                               beta_threshold =0.5,plot=T,save_path = save_path)

celltype_colors = c("BE" = "#E00A82","SQ" = "#ea6860","IM" = "#538e8e","FB" = "#8bc3c3","VC" = "#a2b7d7")
celltypes = c("BE" ,"SQ" ,"IM","FB","VC")
p<-plot_lineage_significance(
  res = res_lg,
  celltypes=celltypes,
  fill_values=celltype_colors,
  order_by = "BE",
  title="Global Cell Type Test",
  outfile = file.path(save_path, "lineage_significance_allcelltype.pdf"),
  return_plot=T
)
non_log_p <- plot_lineage_significance(
  res = res_lg,
  celltypes=celltypes,
  fill_values=celltype_colors,
  order_by = "BE",
  title="Global Cell Type Test",
  outfile = file.path(save_path, "lineage_significance_allcelltype.pdf"),
  return_plot=T,
  logy=F
)


png("example_data/results/lineage_significance_allcelltype.png",width=1000,height=800)
print(p)
dev.off()

celltype_pattern = paste(celltypes,collapse = "|")
effect_sizes = seq(0.1, 1, by = 0.2)
plot_power_curves(
  beta_list = beta_list,
  effect_sizes = effect_sizes,
  celltype_pattern=celltype_patterns,
  desired_order = c("3054_G>C","3071_T>C","15777_G>C"),
  save_pdf = TRUE,
  save_path = save_path,
  file_name = "Power_curve.pdf",
  return_plot=T
)
png("example_data/results/power_curve.png",width=1000,height=800)
print(p2)
dev.off()

#' Regression Test for Identifying Significant Mutations and Carrier Celltypes.
#'
#' `celltype_test()` implements the statistical tests for identifying significant celltypes given spatial decomposition and varaint information.
#'
#' @param celltypes A string vector containing all unique celltypes.
#' @param voi A string vector containing names of variant of interest. Usually a subset of variants that went through basic filtering.
#' @param N A MT variant x cell total coverage matrix.
#' @param X A MT variant x cell alternative allele count matrix. Optional if vaf matrix is provided.
#' @param vaf A MT variant x cell variant allele frequnecy (VAF) matrix. Equivalent as af.dm matrix.
#' @param Ws A cell x celltype weight matrix. Representing spatial decomposition results from tools like RCTD.
#' @param spatial_coords A cell x spatial coordinates matrix. Representing the actual spot/cell location in space.
#' @param test_type Type of Celltype test. Options from c("linear","weighted").
#' @param permute_num A numeric number representing permuation number for each variant.
#' @param save_path path for output results
#' @param figure_height output pdf figure height
#' @param figure_width output pdf figure width
#' @param method P-value correction method: Options: c("FDR", "FWER", "Raw").
#'
#' @param ... Other options used to control matching behavior between duplicate
#'   strings. Passed on to [stringi::stri_opts_collator()].
#' @returns A list object containing fitted parameters (intercept, coefficient and p-values).
#' @import Matrix
#'
#' @examples
#' celltype_test(celltypes=c("BE","DYSP","Normal"), voi=voi_list,
#' N=N,vaf=af.dm, test_type="linear",permute_num=1000)
#'
#' @export
#'
celltype_test <- function(celltypes=NULL,voi=NULL,N=NULL,vaf=NULL,X=NULL,Ws=NULL,spatial_coords=NULL,
                          test_type=c("linear","weighted"),permute_num=1000,plot=F,save_path=NULL,
                          figure_height=10,figure_width=10,method="Raw"){
  # check for input and conditions
  if(test_type == "linear"){
    message("Performing Linear Rgression ...")
  }else if(test_type == "weighted"){
    message("Performing Linear Rgression weighted by coverage...")
  }else{
    error("Not a valid test. Please choose from documented tests.")
  }

  if(is.null(vaf)){
    if(!is.null(X) && !is.null(N)){
      vaf = X / N
    }else{
      stop("Missing variant allele frequency matrix 'vaf'; to compute it instead, both the alternative allele count matrix 'X' and coverage matrix 'N' must be supplied.")
    }
  }
  # check if celltypes aligned with weight matrix
  if(is.null(celltypes)){
    message("Missing celltypes...Set as the column names of celltype weight matrix Ws.")
    celltypes = colnames(Ws)
  }else if(any(celltypes != colnames(Ws))){
    error("Celltypes and the column names of celltype weight matrix Ws not matched...Please check the variables celltypes or Ws.")
  }
  if(is.null(save_path)){
    message("Default saving results to 'output' folder...")
    save_path = "output/"
  }else{
    dir.create(save_path)
  }

  # create plotting data frame
  if(plot){
    intersect_bc=intersect(intersect(colnames(N),rownames(Ws)),rownames(spatial_coords))
    plot_df=cbind(spatial_coords[intersect_bc,],
                  Ws[intersect_bc,],
                  as.data.frame(t(vaf[,intersect_bc])))
    # add in coverage of variant of interest
    cov = N[,intersect_bc];rownames(cov) =  paste0(rownames(N),"_cov")
    plot_df = cbind(plot_df,t(cov))
  }

  # for each variant j, fit a Xs/Ns vs Ws Linear regression
  intercept_df = Matrix(NA,nrow=length(voi),ncol=length(celltypes))
  coef_df = Matrix(NA,nrow=length(voi),ncol=length(celltypes))
  pval_df=Matrix(NA,nrow=length(voi),ncol=length(celltypes))
  adjusted_pval_df <- Matrix(NA_real_, nrow = length(voi), ncol = length(celltypes),sparse=F)

  rownames(intercept_df) <- voi
  colnames(intercept_df) <- celltypes
  rownames(coef_df) <- voi
  colnames(coef_df) <- celltypes
  rownames(pval_df) <- voi
  colnames(pval_df) <- celltypes
  rownames(adjusted_pval_df) <- voi
  colnames(adjusted_pval_df) <- celltypes

  all_pvals <- c()  # Store all p-values for multiple testing correction

  for(j in 1:length(voi)){
    var = voi[j]
    print(paste0(j,": ", var))
    N_j = N[j,]
    vaf_j = vaf[j,]

    for(k in 1:length(celltypes)){ # for each celltype

      data = data.frame(vaf_j=vaf_j,N_j,W_sk=Ws[,k])
      cell_ratio_bins = seq(from=0,to=1,by=0.05)

      if(test_type == "linear"){
        fit=lm(vaf_j ~ W_sk,data=data)# weights equals to inverse of variance (1/sigma_i^2)
      }else if(test_type == "weighted"){
        fit=lm(vaf_j ~ W_sk,data=data, weights=sqrt(N_j))# weights equals to inverse of variance (1/sigma_i^2)
      }
      res=summary(fit)
      intercept_df[j,k]=res$coefficients[1,1]
      coef_df[j,k]= res$coefficients[2,1]

      if(test_type == "linear"){
        pval_df[j,k] = res$coefficients[2,4]
      }else if(test_type == "weighted"){
        # doing permutations and see where the data locates
        set.seed(42)
        coef_list <- c()
        #  decouple VAF with Ws
        for(l in 1:permute_num){
          idx=sample(1:length(vaf_j),size=length(vaf_j))
          vaf_shuffled = vaf_j[idx]
          #vaf_shuffled = pmin(X_shuffled/(N_j+1e-6),1) # cap max VAF to be by 1
          fit=lm(vaf_shuffled ~ Ws[,k],weights=sqrt(N_j[idx]))# weights equals to inverse of variance (1/sigma_i^2)
          res_shuf=summary(fit)
          coef_list <- append(coef_list,res_shuf$coefficients[2,1])
        }
        pval_df[j,k] = 1- sum(res$coefficients[2,1] >= coef_list)/length(coef_list)
      }
    }
    # Collect p-values
    all_pvals <- c(all_pvals, pval_df[j, ])
  }

  # Apply Multiple Testing Correction
  if (method == "FDR") {
    all_pvals_adj <- p.adjust(all_pvals, method = "fdr")  # FDR correction
    message("Applying FDR correction to p-values.")
  } else if (method == "FWER") {
    all_pvals_adj <- p.adjust(all_pvals, method = "bonferroni")  # Bonferroni for FWER
    message("Applying FWER correction to p-values.")
  } else {
    all_pvals_adj <- all_pvals  # Keep raw p-values
    message("Using raw p-values.")
  }

  # Restore adjusted p-values into `adjusted_pval_df`
  index <- 1
  for(j in 1:length(voi)){
    adjusted_pval_df[j,1:length(celltypes)] = as.numeric(all_pvals_adj[index:(index + length(celltypes) - 1)])
    index <- index + length(celltypes)
  }


  rownames(adjusted_pval_df) <- voi
  colnames(adjusted_pval_df) <- celltypes

  if(plot){
    message("Plotting raw p-values...")
    pdf(paste0(save_path,"/",test_type,"_regression_VAF_vs_celltype_plot.pdf"),
        height=figure_height,width=figure_width)
    # for variant j, plot out the fitted values on the diagnostic plots
    for(j in 1:length(voi)){
      intercept=intercept_df[j,]
      coef=coef_df[j,]
      pval=pval_df[j,]
      plot_vaf_cellprop(j,vaf,Ws,plot_df,intercept,coef,pval)
    }
    dev.off()
    # plot adjusted p-value as well
    if(method!="Raw"){
      message("Plotting adjusted p-values")
      pdf(paste0(save_path,"/",test_type,"_regression_VAF_vs_celltype_plot_",method,"_adjustedP.pdf"),
          height=figure_height,width=figure_width)
      # for variant j, plot out the fitted values on the diagnostic plots
      for(j in 1:length(voi)){
        intercept=intercept_df[j,]
        coef=coef_df[j,]
        pval=adjusted_pval_df[j,]
        plot_vaf_cellprop(j,vaf,Ws,plot_df,intercept,coef,pval)
      }
      dev.off()
    }
    message("Diagnostic plots generated.")
  }

  return(res=list(intercept=intercept_df,coef=coef_df,pval=pval_df,
                  adjusted_pval = adjusted_pval_df,adjusted_method=method))
}

#' Filter Significant Variants Based on Adjusted P-values and Coefficients
#'
#' This function takes a results object (from a celltype testing analysis) that contains raw and adjusted
#' p-values as well as coefficients for each variant across multiple cell types. It filters for variants
#' with an adjusted p-value below a specified alpha threshold and a positive coefficient, and then saves
#' the significant variants for each cell type into a text file.
#'
#' @param results A list object returned by \code{celltype_test()} that contains matrices:
#'   \code{adjusted_pval} (adjusted p-values) and \code{coef} (regression coefficients) with rows as variants
#'   and columns as cell types.
#' @param alpha_threshold A numeric value specifying the significance threshold for adjusted p-values (default is 0.05).
#' @param output_file A character string specifying the filename where the significant variants will be saved (default is "significant_variants.txt").
#'
#' @return A list where each element corresponds to a cell type and contains the names of the variants that meet the criteria.
#'
#' @examples
#' \dontrun{
#'   # With 'results'object returned from celltype_test():
#'   sig_vars <- filter_significant_variants(results, alpha_threshold = 0.05, output_file = "output/sig_variants.txt")
#' }
#'
#' @export
filter_significant_variants <- function(results, alpha_threshold = 0.05, output_file = "output/significant_variants.txt") {
  # 'results' should be a list with at least:
  #   results$adjusted_pval: a matrix with variants as rows and cell types as columns
  #   results$coef: a matrix with variants as rows and cell types as columns

  # Get the adjusted p-values and coefficient matrices
  adjusted_pval <- results$adjusted_pval
  coef_mat <- results$coef

  # Get variant names and cell type names
  variant_names <- rownames(adjusted_pval)
  celltypes <- colnames(adjusted_pval)

  # Create a list to store significant variants per cell type
  sig_variants <- list()

  # Loop over each cell type and filter variants
  for (celltype_i in celltypes) {
    sig_idx <- which(adjusted_pval[, celltype_i] < alpha_threshold & coef_mat[, celltype_i] > 0)
    if (length(sig_idx) > 0) {
      sig_variants[[celltype_i]] <- variant_names[sig_idx]
    } else {
      sig_variants[[celltype_i]] <- character(0)
    }
  }

  # Write the significant variants to a text file
  message("Saving significant variants for each celltype..")
  file_conn <- file(output_file, "w")
  for (celltype_i in names(sig_variants)) {
    writeLines(paste("Cell Type:", celltype_i), file_conn)
    if (length(sig_variants[[celltype_i]]) > 0) {
      writeLines(sig_variants[[celltype_i]], file_conn)
    } else {
      writeLines("No significant variants found.", file_conn)
    }
    writeLines("----------", file_conn)
  }
  close(file_conn)

  # Return the list of significant variants for further use if desired
  return(sig_variants)
}

#' Power Analysis for Identifying Negative Mutations and Carrier Celltypes.
#' Simulate alternative allele count values Xi through Binomial distribution Xi ~ Binom(Ni,beta_0 + beta_1 * Wik)
#'
#' @param beta A numeric value. Equivalent as beta_1 in the simulation.
#' @param Nj A total coverage counts for a given variant j
#' @param X A MT variant x cell alternative allele count matrix. Optional if vaf matrix is provided.
#' @param vaf A MT variant x cell variant allele frequnecy (VAF) matrix. Equivalent as af.dm matrix.
#' @param Ws A cell x celltype weight matrix. Representing spatial decomposition results from tools like RCTD.
#' @param spatial_coords A cell x spatial coordinates matrix. Representing the actual spot/cell location in space.
#' @param test_type Type of Celltype test. Options from c("linear","weighted").
#' @param permute_num A numeric number representing permuation number for each variant.
#'
calc_power <- function(beta,Nj,Wj,vaf_j,alpha,n_sim,verbose=T,report_frac=3){
  rejections = 0
  p_vals = c()
  for(i in 1:n_sim){
    if(verbose){
      if(i%%(round(n_sim/report_frac))){
       print(paste0("sim:",i,"/",n_sim))
      }
    }
    X_sim = rbinom(n=length(Nj),size=Nj,prob=beta*Wj)
    data = data.frame(vaf_j=X_sim/(Nj+0.001),Nj,Wj=Wj)
    fit=lm(vaf_j ~ Wj,data=data)
    res=summary(fit)
    p_vals = append(p_vals,res$coefficients[2,4])
  }
  rejections  = sum(p_vals <= alpha,na.rm=T)
  power = rejections/(n_sim)
  return(power)
}

#' Power Analysis for Identifying Negative Mutations and Carrier Celltypes.
#' for all variants of interest & all celltypes
#' set permute_num=20 for now, shall be 1000
#' @export
power_analysis_all <- function(voi=NULL,celltypes=NULL,Ws=NULL,N=NULL,vaf=NULL,X=NULL,
                               sample_num=100,alpha=0.05,n_sim=10,beta_threshold=0.5,plot=T,save_path=NULL,
                               figure_height=8,figure_width=8){
  if(is.null(save_path)){
    save_path = "output/"
    message("Setting output folder to be 'outoput'")
    dir.create(save_path)

  }else{
    dir.create(save_path)
  }
  message("Running Power Analysis for Negative Celltypes and Mutations")
  intersect_bc = intersect(colnames(N),rownames(Ws))
  beta_list = c()

  if(plot){
    pdf( paste0(save_path,"/example_power_analysis.pdf"),height=figure_height,width=figure_width)
  }
  for(j in 1:length(voi)){
    var = voi[j];
    print(paste0("voi ",j,":",var))
    N_j = N[j,intersect_bc]
    vaf_j = vaf[j,intersect_bc]
    if(plot){
      # Dynamically set the plotting grid to fit the number of cell types
      par(mfrow = c(ceiling(sqrt(length(celltypes))), ceiling(sqrt(length(celltypes)))))
    }
    for(k in 1:length(celltypes)){
      sample_idx = sample(intersect_bc,sample_num)
      Nj = N_j[sample_idx] # sampled coverage
      Wj = Ws[sample_idx,k] # celltype ratio
      beta_grid <- seq(0.1,1,by=0.2)

      power_results = sapply(beta_grid,function(beta){
        #print(beta)
        calc_power(beta,Nj,Wj,vaf_j,alpha,n_sim,verbose=F)
      })
      names(power_results) = beta_grid
      beta_list[[paste0(voi[j],"_",celltypes[k])]] <- power_results

      if(plot){
        plot(beta_grid,power_results,type="l",col="blue",lwd=2,
             xlab="Beta",ylab="Power",main=paste0(var,", ",celltypes[k]),
             ylim=c(-0.1,1),xlim=c(0,1)
        )
        points(c(0.1,0.5,1),power_results[as.character(c(0.1,0.5,1))],pch=16)
        abline(h=0.8,col="orange",lty=2)
      }
    }
  }
  if(plot){
    dev.off()
  }
  return(beta_list)

}


# Grid based test - Scan statistics
#' Perform Celltype Test Within Spatial Grids
#'
#' This function partitions the spatial coordinates into a 20x20 grid and performs
#' the linear regression test within each grid to identify significant cell types.
#'
#' @param celltypes A string vector of unique cell types.
#' @param voi A string vector of variants of interest.
#' @param N A matrix of total coverage counts (MT variant x cell).
#' @param vaf A matrix of variant allele frequencies (MT variant x cell).
#' @param Ws A matrix of cell type weights (cell x celltype).
#' @param spatial_coords A matrix of spatial coordinates (cell x 2).
#' @param test_type A string specifying test type: "linear" or "weighted".
#' @param permute_num Number of permutations for weighted test.
#' @param grid_size Number of bins per axis (default = 20 for 20x20 grid).
#' @param method P-value correction methods. Options from c("FDR","FWER","Raw").
#' @param sample_idx Spots as control for each grid's regression. Usually low diseased celltype proportions and not low coverage spots.
#'
#' @return A nested list of p-values for each grid and variant.
#' @export
celltype_test_grid <- function(celltypes, voi, N_voi, vaf, Ws, spatial_coords,
                               test_type = c("linear", "weighted"),
                               permute_num = 1000, grid_size= 20,verbose=F,
                               vaf_cellprop=T,power=T,disease_celltype=NULL,
                               method="Raw",sample_idx = NULL, min_cell_grid=20) {
  # Find common cells across datasets
  intersect_bc <- intersect(intersect(colnames(N_voi), rownames(Ws)), rownames(spatial_coords))

  # Ensure spatial_coords is a matrix/dataframe
  spatial_coords <- as.data.frame(spatial_coords[intersect_bc,])
  colnames(spatial_coords) <- c("X", "Y") # Ensure correct column names

  # Create Grid Labels
  spatial_coords$X_bin <- cut(spatial_coords$X, breaks = grid_size, labels = FALSE)
  spatial_coords$Y_bin <- cut(spatial_coords$Y, breaks = grid_size, labels = FALSE)

  # Create Unique Grid IDs
  spatial_coords$Grid_ID <- paste0("Grid_", spatial_coords$X_bin, "_", spatial_coords$Y_bin)

  # Merge Spatial Data with Weights and VAF
  plot_df <- cbind(spatial_coords[intersect_bc, ],
                   Ws[intersect_bc, ],
                   vaf=as.data.frame(vaf[, intersect_bc]))
  colnames(plot_df)[length(colnames(plot_df))] = voi

  Ws = Ws[intersect_bc, ]
  vaf = vaf[, intersect_bc,drop=F]
  N_voi = N_voi[,intersect_bc,drop=F]

  # Add Coverage Data
  #cov <- N[, intersect_bc]
  cov <- N_voi[,intersect_bc,drop=F]
  rownames(cov) <- paste0(rownames(N_voi), "_cov")
  plot_df <- cbind(plot_df, t(cov))

  # Initialize Output
  results <- list()
  grid_vaf_cellprop_plot <- list()
  if(method !=" Raw"){
    all_pvals <- c()  # Store all p-values for correction later
  }

  # Loop through each grid
  for (grid in unique(spatial_coords$Grid_ID)) {
    if(verbose){
      message(paste("Processing:", grid))
    }

    # Subset data for the current grid
    grid_cells <- spatial_coords$Grid_ID == grid
    if (sum(grid_cells) < min_cell_grid) next  # Skip empty grids

    N_grid <- N_voi[, grid_cells, drop = FALSE]
    vaf_grid <- vaf[, grid_cells, drop = FALSE]
    Ws_grid <- Ws[grid_cells, , drop = FALSE]


    # Initialize Data Storage
    intercept_df <- Matrix(NA, nrow = length(voi), ncol = length(celltypes))
    coef_df <- Matrix(NA, nrow = length(voi), ncol = length(celltypes))
    pval_df <- Matrix(NA, nrow = length(voi), ncol = length(celltypes))

    rownames(intercept_df) <- voi
    colnames(intercept_df) <- celltypes
    rownames(coef_df) <- voi
    colnames(coef_df) <- celltypes
    rownames(pval_df) <- voi
    colnames(pval_df) <- celltypes

    # for adjusted p-value
    adjusted_pval_df <- Matrix(NA, nrow = length(voi), ncol = length(celltypes))
    rownames(adjusted_pval_df) <- voi
    colnames(adjusted_pval_df) <- celltypes


    # Loop through each variant
    for (j in 1:length(voi)) {
      var <- voi[j]
      if(verbose){
        message(paste("  Variant:", var))
      }

      N_j <- N_grid[j, ]
      vaf_j <- vaf_grid[j, ]

      # low diseased proportion but high coverage
      if(is.null(sample_idx)){
        sample_idx = intersect(rownames(Ws)[(Ws[,disease_celltype] < 0.3) & !is.na(Ws[,disease_celltype])],
                               colnames(N_voi)[N_voi[var,] > 1])
      }else{
        message("Using user input spots as control: ", length(sample_idx)," spots")
      }
      # if(verbose){
      #   print(paste0("low prop spots:",length(rownames(Ws)[(Ws[,disease_celltype] < 0.3) & !is.na(Ws[,disease_celltype])])))
      #   print(paste0("high cov spots:",length(colnames(N_voi)[N_voi[var,] > 1])))
      # }
      #intersect(rownames(Ws)[(Ws[,disease_celltype] < 0.2) & !is.na(Ws[,disease_celltype])],
      #          colnames(N)[N[var,] > 2])

        #colnames(vaf)[(N_voi[j,] <5) & (N_voi[j,] >0) ] # low coverage spots
        #sample(colnames(vaf)[(N_voi[j,] <5 )&(N_voi[j,] >0)],1000)


      N_j <- append(N_j,N_voi[j,sample_idx])
      vaf_j <- append(vaf_j,vaf[j,sample_idx])


      for (k in 1:length(celltypes)) { # Loop through each cell type
        Ws_jk = append(Ws_grid[, k], Ws[sample_idx,k])
        data <- data.frame(vaf_j = vaf_j, N_j = N_j, W_sk = Ws_jk )

        # Check if there is enough variation in the predictor
        if (length(unique(data$W_sk)) == 1 || length(unique(data$vaf_j)) == 1) {
          message(paste("Skipping grid", grid, "for", var, "-", celltypes[k], "due to insufficient variation"))
          next
        }

        # Perform Linear Regression Safely
        if (test_type == "linear") {
          fit <- tryCatch(lm(vaf_j ~ W_sk, data = data), error = function(e) NULL)
          permute=F
        } else if (test_type == "weighted") {
          fit <- tryCatch(lm(vaf_j ~ W_sk, data = data, weights = sqrt(N_j)), error = function(e) NULL)
          permute=T
        }

        # Skip if fit fails
        if (is.null(fit)) {
          message(paste("Skipping grid", grid, "for", var, "-", celltypes[k], "due to regression failure"))
          next
        }

        # Store Results
        res <- summary(fit)

        # Ensure coefficients exist
        if (nrow(res$coefficients) < 2) {
          message(paste("Skipping grid", grid, "for", var, "-", celltypes[k], "due to missing coefficients"))
          next
        }

        intercept_df[j, k] <- res$coefficients[1, 1]
        coef_df[j, k] <- res$coefficients[2, 1]
        pval_df[j, k] <- res$coefficients[2, 4]

        # Collect p-values for correction
        if(method!="Raw"){
          all_pvals <- c(all_pvals, pval_df[j, k])
        }

        #if(power){}
      }
      # plot vaf vs celltype diagnostic for grid
      if(vaf_cellprop){
        intercept=intercept_df[j,]
        coef=coef_df[j,]
        pval=pval_df[j,]
        print(pval)
        temp_plot_df = plot_df[grid_cells,];
        grid_vaf_cellprop_plot[[grid]] =
          plot_vaf_cellprop(j,vaf_grid,Ws_grid,temp_plot_df ,
                          intercept,coef,pval,permute,return_plot=T)
      }
    }

    # Store grid results
    results[[grid]] <- list(intercept = intercept_df, coef = coef_df, pval = pval_df,
                            adjusted_pval = adjusted_pval_df)
  }

  ### **Apply P-Value Adjustment After All Grids Have Been Processed**
  if (method == "FDR") {
    all_pvals_adj <- p.adjust(all_pvals, method = "fdr")  # FDR correction
    message("Applying FDR correction to p-values.")
  } else if (method == "FWER") {
    all_pvals_adj <- p.adjust(all_pvals, method = "bonferroni")  # Bonferroni for FWER
    message("Applying FWER correction to p-values.")
  } else if(method == "Raw"){
    all_pvals_adj <- all_pvals  # Keep raw p-values
    message("Using raw p-values.")
  }

  # Restore adjusted p-values into `adjusted_pval_df`
  index <- 1
  for (grid in names(results)) {
    results[[grid]]$adjusted_pval[1,1:length(celltypes)] <- as.numeric(all_pvals_adj[index:(index + length(celltypes) - 1)])
    index <- index + length(celltypes)
  }

  if(vaf_cellprop){
    return(list(results=results,vaf_cellprop_plots=grid_vaf_cellprop_plot,
                spatial_coords_grids=spatial_coords))
  }else{
    return(list(results=results,spatial_coords_grids=spatial_coords))
  }
}

#' Plot P-Values from Grid-Based Regression for Each Cell Type
#'
#' This function plots the -log10 p-values for each cell type and variant
#' within the spatial grid, using a continuous color scale.
#'
#' @param results_grid A list containing p-values for each spatial grid.
#' @param grid_size Number of bins per axis (default = 20 for 20x20 grid).
#' @param coef_plot_option A string from c("grey","negative"),
#' @param max_log10p_cap Maximum cap for -log10(p-values) (default = 5).
#' @param p_thresh Threshold for significance (default = 0.05).
#' indicating if plotting the negative coeffient p-value as 1 or reversed sign.
#' Default "grey", setting to p-val = 1 (-log10p-val = 0).
#'
#' @return A grid plot showing p-value distributions for each cell type.
#' @import ggplot2 dplyr tidyr cowplot
#' @export
#'
plot_pval_grid_per_celltype <- function(results_grid, grid_size = 20,
                                        coef_plot_option="grey",
                                        method="Raw",
                                        max_log10p_cap=5,
                                        p_thresh=0.05){

  # Extract Grid IDs
  grid_ids <- names(results_grid)
  all_celltypes <- colnames(results_grid[[grid_ids[1]]]$pval)  # Extract cell types

  # Initialize List to Store Plots
  plot_list <- list()

  # Loop Through Each Cell Type and Create a Plot
  for (celltype in all_celltypes) {
    if(method=="Raw"){
      grid_data <- data.frame(
        Grid_ID = rep(grid_ids, each = length(results_grid[[1]]$pval[, celltype])),
        Variant = rep(rownames(results_grid[[1]]$pval), times = length(grid_ids)),
        X_bin = as.numeric(sapply(strsplit(grid_ids, "_"), "[[", 2)),
        Y_bin = as.numeric(sapply(strsplit(grid_ids, "_"), "[[", 3)),
        coef = unlist(lapply(results_grid, function(res) res$coef[, celltype])),
        Pval = unlist(lapply(results_grid, function(res) res$pval[, celltype]))
      )
    }else if (method %in% c("FDR","FWER")){
      grid_data <- data.frame(
        Grid_ID = rep(grid_ids, each = length(results_grid[[1]]$pval[, celltype])),
        Variant = rep(rownames(results_grid[[1]]$pval), times = length(grid_ids)),
        X_bin = as.numeric(sapply(strsplit(grid_ids, "_"), "[[", 2)),
        Y_bin = as.numeric(sapply(strsplit(grid_ids, "_"), "[[", 3)),
        coef = unlist(lapply(results_grid, function(res) res$coef[, celltype])),
        Pval = unlist(lapply(results_grid, function(res) res$adjusted_pval[, celltype]))
      )
    }

    # Drop missing values (NA p-values)
    grid_data <- grid_data %>% filter(!is.na(Pval))

    # If no valid data, skip this celltype
    if (nrow(grid_data) == 0) {
      message(paste("Skipping cell type", celltype, "due to missing data"))
      next
    }

    # Transform p-values: -log10 and cap at 4
    grid_data$Log_Pval <- -log10(grid_data$Pval)
    grid_data$Log_Pval[grid_data$Log_Pval > max_log10p_cap] <- max_log10p_cap  # Cap log value
    if(coef_plot_option=="grey"){
      grid_data$Log_Pval[grid_data$coef < 0] = 0
    }else if(coef_plot_option=="negative"){
      grid_data$Log_Pval[grid_data$coef < 0] = -grid_data$Log_Pval[grid_data$coef < 0]
    }



    # Define Continuous Color Scale
    if(coef_plot_option=="grey"){
      pval_palette <- scale_fill_gradientn(
        colors = c("darkgrey","grey","red"),   # Grey for low p-values, red for high p-values
        values = scales::rescale(c(0, -log10(p_thresh), max_log10p_cap)),  # Map p-values to colors
        limits = c(0, max_log10p_cap),  # Ensures consistent scale
        na.value = "white"  # Assigns white for missing values
      )
    }else if(coef_plot_option=="negative"){
      pval_palette <- scale_fill_gradientn(
        colors = c("blue","grey", "red"),   # Grey for low p-values, red for high p-values
        values = scales::rescale(c(-max_log10p_cap, 0, max_log10p_cap)),  # Map p-values to colors
        limits = c(-max_log10p_cap, max_log10p_cap),  # Ensures consistent scale
        na.value = "white"  # Assigns white for missing values
      )
    }

    # Generate Plots for Each Variant
    for (variant in unique(grid_data$Variant)) {
      df_variant <- grid_data %>% filter(Variant == variant)

      # Skip empty data
      if (nrow(df_variant) == 0) {
        message(paste("Skipping", variant, "for", celltype, "due to missing data"))
        next
      }

      p <- ggplot(df_variant, aes(x = X_bin, y = Y_bin, fill = Log_Pval)) +
        geom_tile(color = "white") +  # Grid Cells
        pval_palette +
        labs(title = paste(variant, "-", celltype, "-",method, " P value"),
             x = "Grid X", y = "Grid Y", fill = "-log10(P-value)") +
        theme_minimal() +
        coord_fixed()  # Keep aspect ratio square

      # Store plot in list only if it has data
      plot_list[[paste(variant, celltype, sep = "_")]] <- p
    }
  }

  # Ensure at least one plot exists
  if (length(plot_list) == 0) {
    message("No valid plots available")
    return(NULL)
  }

  # Arrange Plots in a Grid
  combined_plot <- plot_grid(plotlist = plot_list, ncol = ceiling(sqrt(length(all_celltypes))))  # Adjust ncol if needed

  return(combined_plot)
}


#' Compute Gene Program Enrichment per Grid
#'
#' @param seu Seurat object containing gene expression data.
#' @param results_grid A list containing grid-based cell assignments.
#' @param gene_programs A named list where each element is a vector of genes for a specific program.
#'
#' @return A list containing enrichment scores for each grid.
#' @export
compute_gene_program_enrichment <- function(seu, results_grid, spatial_coords_grids, gene_programs) {

  # Initialize storage list
  enrichment_scores <- list()
  # Extract Grid IDs
  grid_ids <- names(results_grid)

  # Loop through each grid
  for (grid in grid_ids) {

    # Get barcodes in this grid
    barcodes <- rownames(spatial_coords_grids)[spatial_coords_grids$Grid_ID == grid]

    # Skip empty grids
    if (length(barcodes) == 0) {
      message(paste("Skipping", grid, "- No barcodes found"))
      next
    }

    # Subset Seurat object to these barcodes
    seu_subset <- subset(seu, cells = barcodes)

    # Calculate module scores for gene programs
    #genes=intersect(gene_programs, rownames(seu_subset))
    seu_subset <- AddModuleScore(seu_subset, features = gene_programs,
                                 name=names(gene_programs),
                                 #name = "GeneProgram",
                                 assay="RNA",nbin=10)

    # Store scores in list
    enrichment_scores[[grid]] <- FetchData(seu_subset, vars=paste0(names(gene_programs),seq_along(gene_programs) ))
                                           #vars =  paste0("GeneProgram", seq_along(gene_programs)))
  }

  # Convert enrichment results to a plottable data frame
  plot_data <- do.call(rbind, lapply(names(enrichment_scores), function(grid) {
    scores <- enrichment_scores[[grid]]
    if (is.null(scores)) return(NULL)  # Skip if empty

    # Compute mean enrichment per grid
    mean_scores <- colMeans(scores, na.rm = TRUE)

    # Extract grid positions
    x_bin <- as.numeric(strsplit(grid, "_")[[1]][2])
    y_bin <- as.numeric(strsplit(grid, "_")[[1]][3])

    # Convert to data frame
    data.frame(Grid_ID = grid, X_bin = x_bin, Y_bin = y_bin, t(mean_scores))
  }))

  # Reshape to long format for ggplot
  plot_data_long <- pivot_longer(plot_data, cols =paste0(names(gene_programs), seq_along(gene_programs))
                                 , names_to = "GeneProgram", values_to = "Score")

  # Generate heatmap of gene program enrichment
  enrichment_plot <- ggplot(plot_data_long, aes(x = X_bin, y = Y_bin, fill = Score)) +
    geom_tile(color = "white") +  # Grid cells
    scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, na.value = "white") +  # Color scale
    facet_wrap(~GeneProgram) +  # One plot per gene program
    labs(title = "Gene Program Enrichment Across Spatial Grid",
         x = "Grid X", y = "Grid Y", fill = "Enrichment Score") +
    theme_minimal() +
    coord_fixed()  # Keep aspect ratio square

  return(list(scores = enrichment_scores, plot = enrichment_plot))
}


# Grid based test
# Option 1 - move by a small fraction per time
celltype_test_continuous <- function(celltypes, voi, N_voi, vaf, Ws, spatial_coords,
                                     test_type = c("linear", "weighted"),
                                     permute_num = 1000, grid_size = 20, step_fraction = 1/50,
                                     verbose = FALSE, vaf_cellprop = TRUE, power = TRUE,
                                     disease_celltype = NULL, method = "Raw", sample_idx = NULL,
                                     min_cell_grid = 20) {
  # Find common cells across datasets
  intersect_bc <- intersect(intersect(colnames(N_voi), rownames(Ws)), rownames(spatial_coords))
  spatial_coords <- as.data.frame(spatial_coords[intersect_bc, ])
  colnames(spatial_coords) <- c("X", "Y")

  # Determine overall spatial range and window size
  x_min <- min(spatial_coords$X)
  x_max <- max(spatial_coords$X)
  y_min <- min(spatial_coords$Y)
  y_max <- max(spatial_coords$Y)

  #win_width_x <- (x_max - x_min) / grid_size
  #win_width_y <- (y_max - y_min) / grid_size

  # Define step sizes (moving window)
  #step_x <- win_width_x * step_fraction
  #step_y <- win_width_y * step_fraction
  step_x = (x_max - x_min) * step_fraction
  step_y = (y_max - y_min) * step_fraction

  # Initialize output list
  results <- list()
  grid_vaf_cellprop_plot <- list()  # (optional diagnostic plots)
  if(method != "Raw"){
    all_pvals <- c()  # to collect raw p-values (if correction is applied externally)
  }

  window_id <- 1
  # Loop over continuous window positions
  for (x_start in seq(x_min, x_max - win_width_x, by = step_x)) {
    for (y_start in seq(y_min, y_max - win_width_y, by = step_y)) {
      # Define window boundaries
      x_end <- x_start + win_width_x
      y_end <- y_start + win_width_y

      # Identify cells in the current window
      window_cells <- rownames(spatial_coords)[
        spatial_coords$X >= x_start & spatial_coords$X <= x_end &
          spatial_coords$Y >= y_start & spatial_coords$Y <= y_end
      ]
      # Skip if fewer than min_cell_grid cells
      if (length(window_cells) < min_cell_grid) next

      # Create a unique window ID (e.g., "Window_123.45_678.90")
      win_id <- paste0("Window_", format(x_start, digits = 3), "_", format(y_start, digits = 3))

      # Subset matrices for these cells
      N_win <- N_voi[, window_cells, drop = FALSE]
      vaf_win <- vaf[, window_cells, drop = FALSE]
      Ws_win <- Ws[window_cells, , drop = FALSE]

      # Initialize storage for regression results in this window
      intercept_df <- Matrix(NA, nrow = length(voi), ncol = length(celltypes))
      coef_df <- Matrix(NA, nrow = length(voi), ncol = length(celltypes))
      pval_df <- Matrix(NA, nrow = length(voi), ncol = length(celltypes))
      rownames(intercept_df) <- voi
      colnames(intercept_df) <- celltypes
      rownames(coef_df) <- voi
      colnames(coef_df) <- celltypes
      rownames(pval_df) <- voi
      colnames(pval_df) <- celltypes

      # Loop over each variant
      for (j in seq_along(voi)) {
        var <- voi[j]
        if(verbose) message(paste("Window", win_id, "Variant:", var))
        N_j <- N_win[j, ]
        vaf_j <- vaf_win[j, ]

        # If a set of control spots is desired, determine sample_idx_local
        if(is.null(sample_idx) && !is.null(disease_celltype)) {
          sample_idx_local <- intersect(
            rownames(Ws)[(Ws[, disease_celltype] < 0.3) & !is.na(Ws[, disease_celltype])],
            colnames(N_voi)[N_voi[var, ] > 1]
          )
        } else {
          sample_idx_local <- sample_idx
        }
        # Append control spots if available
        if(length(sample_idx_local) > 0){
          N_j <- append(N_j, N_voi[j, sample_idx_local])
          vaf_j <- append(vaf_j, vaf[j, sample_idx_local])
        }

        # Loop over each cell type
        for (k in seq_along(celltypes)) {
          # Append corresponding control weights
          Ws_jk <- append(Ws_win[, k], Ws[sample_idx_local, k])
          data <- data.frame(vaf_j = vaf_j, N_j = N_j, W_sk = Ws_jk)

          # Check for sufficient variation
          if (length(unique(data$W_sk)) == 1 || length(unique(data$vaf_j)) == 1) {
            if(verbose) message(paste("Skipping window", win_id, "for", var, "-", celltypes[k], "insufficient variation"))
            next
          }

          # Perform regression
          if (test_type == "linear") {
            fit <- tryCatch(lm(vaf_j ~ W_sk, data = data), error = function(e) NULL)
            permute <- FALSE
          } else if (test_type == "weighted") {
            fit <- tryCatch(lm(vaf_j ~ W_sk, data = data, weights = sqrt(N_j)), error = function(e) NULL)
            permute <- TRUE
          }
          if (is.null(fit)) {
            if(verbose) message(paste("Skipping window", win_id, "for", var, "-", celltypes[k], "regression failure"))
            next
          }
          res <- summary(fit)
          if(nrow(res$coefficients) < 2) {
            if(verbose) message(paste("Skipping window", win_id, "for", var, "-", celltypes[k], "missing coefficients"))
            next
          }
          intercept_df[j, k] <- res$coefficients[1, 1]
          coef_df[j, k] <- res$coefficients[2, 1]
          pval_df[j, k] <- res$coefficients[2, 4]
          if(test_type != "Raw"){
            all_pvals <- c(all_pvals, pval_df[j, k])
          }
        } # end cell type loop
        # (Optional: generate diagnostic plots for this variant in the window)
        if(vaf_cellprop) {
          # Example: call a function like plot_vaf_cellprop (not provided here)
          # grid_vaf_cellprop_plot[[win_id]] <- plot_vaf_cellprop(j, vaf_win, Ws_win, ... )
        }
      } # end variant loop

      # Store window results
      results[[win_id]] <- list(intercept = intercept_df,
                                coef = coef_df,
                                pval = pval_df)
      window_id <- window_id + 1
    } # end y_start loop
  } # end x_start loop

  ### Optionally, apply p-value adjustment across windows here if desired
  if(method != "Raw"){
    all_pvals_adj <- p.adjust(all_pvals, method = ifelse(method == "FDR", "fdr", "bonferroni"))
    if(verbose) message(paste("Applied", method, "correction"))
    # For each window, assign adjusted p-values
    index <- 1
    for(win_id in names(results)){
      n_rows <- nrow(results[[win_id]]$pval)
      n_cols <- ncol(results[[win_id]]$pval)
      total <- n_rows * n_cols
      # Fill in a new matrix in row-major order:
      adj_mat <- matrix(all_pvals_adj[index:(index + total - 1)], nrow = n_rows, ncol = n_cols, byrow = TRUE)
      results[[win_id]]$adjusted_pval <- adj_mat
      index <- index + total
    }
  }

  return(list(results = results,spatial_coords=spatial_coords,
              step_size=step_size,grid_size=grids_size))
}

#' Plot Spatial Distribution of Adjusted P-values on a Signed -log10 Scale
#'
#' This function takes a results object returned from a continuous grid-based
#' celltype test. Each grid (window) is assumed to be named in the format
#' "Grid_x_y" (where x and y are the grid indices). For each grid, the function
#' extracts the adjusted p-values and regression coefficients for each variant (rows)
#' and cell type (columns), computes the signed -log10(p-value) (using the sign of the coefficient),
#' and then aggregates the data into a continuous heatmap across the spatial grid.
#'
#' @param results A list returned from a grid-based celltype test. Each element (grid)
#'   should contain matrices \code{adjusted_pval} and \code{coef} with row names as variant IDs
#'   and column names as cell types.
#'
#' @return A \code{ggplot2} object showing the spatial distribution of adjusted p-values
#'   on a signed -log10 scale. Spots with negative coefficients are colored in blue,
#'   positive coefficients in red, and 0 in grey.
#'
#' @examples
#' \dontrun{
#'   # Assume 'results' is the object returned by your grid-based test function.
#'   p <- plot_spatial_logP(results)
#'   print(p)
#' }
#' @import distances
#' @export
plot_spatial_logP <- function(results,var,p_thresh=0.05, max_log10p_cap=6,coef_plot_option="negative") {
  # Get grid IDs; assume they are in the format "Grid_x_y"
  grid_ids <- names(results$results)
  if (is.null(grid_ids)) stop("Results object must have named elements corresponding to grid IDs.")

  # Determine overall spatial range and window size
  spatial_coords_plot=results$spatial_coords
  x_min <- min(spatial_coords_plot$X)
  x_max <- max(spatial_coords_plot$X)
  y_min <- min(spatial_coords_plot$Y)
  y_max <- max(spatial_coords_plot$Y)

  win_width_x <- (x_max - x_min) / grid_size
  win_width_y <- (y_max - y_min) / grid_size

  # Define step sizes (moving window)
  step_x <- win_width_x * step_fraction
  step_y <- win_width_y * step_fraction

  step_x = (x_max - x_min) * step_fraction
  step_y = (y_max - y_min) * step_fraction

  # Parse window centers from window IDs
  window_centers <- do.call(rbind, lapply(grid_ids, function(id) {
    parts <- unlist(strsplit(id, "_"))
    if(length(parts) < 3) stop("Window ID format not recognized. Expected format 'Window_x_y'.")
    x_center <- as.numeric(parts[2]) + step_x/2
    y_center <- as.numeric(parts[3]) + step_y/2
    data.frame(window_id = id, x_center = x_center, y_center = y_center)
  }))

  # For each spot, assign the nearest window (using Euclidean distance)
  assign_window <- function(spot_x, spot_y, centers_df) {
    dists <- sqrt((centers_df$x_center - spot_x)^2 + (centers_df$y_center - spot_y)^2)
    centers_df$window_id[which.min(dists)]
  }

  spatial_coords_plot$assigned_window <- apply(spatial_coords[, c("X", "Y")], 1, function(row) {
    assign_window(row[1], row[2], window_centers)
  })

  # For each spot, extract the adjusted p-value and coefficient for the specified variant and cell type
  get_spot_value <- function(win_id, results,celltype) {
    res <- results$results[[win_id]]
    p_val <- res$pval[, celltype]
    coef_val <- res$coef[, celltype]
    # Compute signed -log10(p)
    signed_logp <- sign(coef_val) * (-log10(p_val))
    return(signed_logp)
  }


  # Define Continuous Color Scale
  if(coef_plot_option=="grey"){
    pval_palette <- scale_color_gradientn(
      colors = c("darkgrey","grey","red"),   # Grey for low p-values, red for high p-values
      values = scales::rescale(c(0, -log10(p_thresh), max_log10p_cap)),  # Map p-values to colors
      limits = c(0, max_log10p_cap),  # Ensures consistent scale
      na.value = "white"  # Assigns white for missing values
    )
  }else if(coef_plot_option=="negative"){
    pval_palette <- scale_color_gradientn(
      colors = c("blue","grey", "red"),   # Grey for low p-values, red for high p-values
      values = scales::rescale(c(-max_log10p_cap, 0, max_log10p_cap)),  # Map p-values to colors
      limits = c(-max_log10p_cap, max_log10p_cap),  # Ensures consistent scale
      na.value = "white"  # Assigns white for missing values
    )
  }

  plot_list  = list()
  # Generate Plots for Each Variant
  for (celltype in celltypes) {
    # Apply the function to each spot
    #spatial_coords_plot$signed_logP
    spatial_coords_plot$signed_logP <- sapply(spatial_coords$assigned_window, function(win_id) {
      get_spot_value(win_id, results, celltype)
    })
    spatial_coords_plot$capped_signed_logP <- pmax(pmin(spatial_coords_plot$signed_logP,
                                                        max_log10p_cap), -max_log10p_cap)

    p <- ggplot(spatial_coords_plot, aes(x = X, y = Y,color = capped_signed_logP)) +
      geom_point(size=1) +
      pval_palette +
      labs(title = paste(var, "-", celltype, "-", " P value"),
           x = "X", y = "Y", fill = "-log10(P-value)") +
      theme_minimal() +
      coord_fixed()  # Keep aspect ratio square

    # Store plot in list only if it has data
    plot_list[[paste(variant, celltype, sep = "_")]] <- p
  }


  # Ensure at least one plot exists
  if (length(plot_list) == 0) {
    message("No valid plots available")
    return(NULL)
  }

  # Arrange Plots in a Grid
  combined_plot <- plot_grid(plotlist = plot_list, ncol = ceiling(sqrt(length(celltypes))))  # Adjust ncol if needed


  return(combined_plot)
}


#' Perform Celltype Test Using Nearest Neighbors for Each Spot
#'
#' This function performs a celltype test for each spot by selecting its nearest k neighbors
#' (including itself) and running a linear (or weighted) regression of variant allele frequency
#' versus cell type weights for each variant of interest. The regression is performed for each cell
#' type. Optionally, p-values across all tests are adjusted using the specified method.
#'
#' @param celltypes A string vector of unique cell types.
#' @param voi A string vector of variants of interest.
#' @param N A matrix of total coverage counts (MT variant x cell).
#' @param vaf A matrix of variant allele frequencies (MT variant x cell).
#' @param Ws A matrix of cell type weights (cell x celltype).
#' @param spatial_coords_temp A matrix or data frame of spatial coordinates (cell x 2) with row names as cell IDs.
#' @param test_type A string specifying the regression test: "linear" or "weighted".
#' @param permute_num Number of permutations for the weighted test (default 1000).
#' @param k_neighbors Number of nearest neighbors to use for each spot (default = 100).
#' @param method P-value correction method. Options are "FDR", "FWER", or "Raw" (default = "Raw").
#'
#' @return A list with one element per spot (named by the spot ID). Each element is a list containing:
#'   \item{intercept}{A matrix of intercept estimates (rows = variants, columns = cell types).}
#'   \item{coef}{A matrix of regression coefficients.}
#'   \item{pval}{A matrix of raw p-values.}
#'   \item{adjusted_pval}{A matrix of adjusted p-values (if method != "Raw"); otherwise, same as raw.}
#'
#' @examples
#' \dontrun{
#'   # results_knn <- celltype_test_knn(celltypes, voi, N, vaf, Ws, spatial_coords_temp,
#'   #                                  test_type = "linear", permute_num = 1000,
#'   #                                  k_neighbors = 100, method = "FDR")
#' }
#'
#' @export
#' @import cowplot
celltype_test_knn <- function(celltypes, vars, N_voi, vaf, Ws, spatial_coords_temp,
                              test_type = c("linear", "weighted"),permute=F,
                              permute_num = 1000, k_neighbors = 100,
                              method = "Raw", sample_idx=NULL,disease_celltype="BE",
                              ratio_threshold=0.03,
                              vaf_cellprop=F,exclude_plot_idx=NULL,
                              p_thresh=0.05,max_log10p_cap=6,coef_plot_option="negative") {
  # # Load necessary package for nearest neighbors
  # if (!requireNamespace("FNN", quietly = TRUE)) {
  #   stop("Please install the 'FNN' package to use this function.")
  # }

  #test_type <- match.arg(test_type)
  # Define Continuous Color Scale
  if(coef_plot_option=="grey"){
    pval_palette <- scale_color_gradientn(
      colors = c("darkgrey","grey","red"),   # Grey for low p-values, red for high p-values
      values = scales::rescale(c(0, -log10(p_thresh), max_log10p_cap)),  # Map p-values to colors
      limits = c(0, max_log10p_cap),  # Ensures consistent scale
      na.value = "white"  # Assigns white for missing values
    )
  }else if(coef_plot_option=="negative"){
    pval_palette <- scale_color_gradientn(
      colors = c("blue","grey", "red"),   # Grey for low p-values, red for high p-values
      values = scales::rescale(c(-max_log10p_cap, 0, max_log10p_cap)),  # Map p-values to colors
      limits = c(-max_log10p_cap, max_log10p_cap),  # Ensures consistent scale
      na.value = "white"  # Assigns white for missing values
    )
  }

  # Find common cells across datasets
  common_cells <- intersect(intersect(colnames(N_voi), rownames(Ws)), rownames(spatial_coords_temp))
  if(length(common_cells) == 0) stop("No common cells found across datasets.")

  N_voi <- N_voi[, common_cells, drop = FALSE]
  vaf <- vaf[, common_cells, drop = FALSE]
  Ws <- Ws[common_cells, , drop = FALSE]
  spatial_coords_temp <- as.data.frame(spatial_coords_temp[common_cells, ])
  colnames(spatial_coords_temp) = c("X","Y")

  # Use FNN::get.knn to obtain indices of k nearest neighbors for each cell
  knn_res <- FNN::get.knn(as.matrix(spatial_coords_temp[, c("X", "Y")]), k = k_neighbors)

  results_list <- list()  # will store results for each spot
  all_pvals <- c()        # to collect p-values for multiple testing correction

  # Loop over each spot (each row in spatial_coords_temp)
  for (i in seq_len(nrow(spatial_coords_temp))) {
    if(i%%1000==0 ){
      print(paste0(i,"/", nrow(spatial_coords_temp)," spots being processed."))
    }
    spot_id <- rownames(spatial_coords_temp)[i]
    neighbor_idx <- knn_res$nn.index[i, ]  # indices of nearest neighbors
    neighbor_ids <- rownames(spatial_coords_temp)[neighbor_idx]
    test_ids = c(neighbor_ids,spot_id ) # include the spot for test

    # Subset matrices for these neighbors
    N_neighbors <- N_voi[, test_ids, drop = FALSE]
    vaf_neighbors <- vaf[, test_ids, drop = FALSE]
    Ws_neighbors <- Ws[test_ids, , drop = FALSE]

    # Initialize matrices for this spot's regression results
    intercept_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(celltypes))
    coef_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(celltypes))
    pval_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(celltypes))
    rownames(intercept_mat) <- vars
    colnames(intercept_mat) <- celltypes
    rownames(coef_mat) <- vars
    colnames(coef_mat) <- celltypes
    rownames(pval_mat) <- vars
    colnames(pval_mat) <- celltypes

    # Loop over each variant
    for (j in seq_along(vars)) {
      var <- vars[j]
      N_j <- N_neighbors[j, ]
      vaf_j <- vaf_neighbors[j, ]
      Ws_j = Ws_neighbors

      if(sum(vaf_j > 0)/length(vaf_j) < ratio_threshold){ # excluding outliers
        intercept_mat[j, ] <- 0
        coef_mat[j, ] <- 0
        pval_mat[j, ] <- 1
        next
      }

      # If a set of control spots is desired, determine sample_idx_local
      if(is.null(sample_idx) && !is.null(disease_celltype)) {
        sample_idx_local <- intersect(
          rownames(Ws)[(Ws[, disease_celltype] < 0.3) & !is.na(Ws[, disease_celltype])],
          colnames(N_voi)[N_voi[var, ] > 1]
        )
      } else {
        sample_idx_local <- sample_idx
      }

      #print(paste0(length(sample_idx_local), " paired control spots."))
      # Append control spots if available
      if(length(sample_idx_local) > 0){
        N_j <- append(N_j, N_voi[j, sample_idx_local])
        vaf_j <- append(vaf_j, vaf[j, sample_idx_local])
        Ws_j = rbind(Ws_neighbors,Ws[sample_idx_local,])
      }

      for (k in seq_along(celltypes)) {
        data_df <- data.frame(vaf_j = vaf_j,
                              N_j = N_j,
                              W_sk = Ws_j[, celltypes[k]])

        # Check for variation in predictors
        if (length(unique(data_df$W_sk)) == 1 || length(unique(data_df$vaf_j)) == 1) {
          next
        }

        # Run regression (linear or weighted)
        if (test_type == "linear") {
          fit <- tryCatch(lm(vaf_j ~ W_sk, data = data_df), error = function(e) NULL)
        } else if (test_type == "weighted") {
          fit <- tryCatch(lm(vaf_j ~ W_sk, data = data_df, weights = sqrt(N_j)), error = function(e) NULL)
        }

        if (is.null(fit)) next

        res <- summary(fit)
        if (nrow(res$coefficients) < 2) next

        intercept_mat[j, k] <- res$coefficients[1, 1]
        coef_mat[j, k] <- res$coefficients[2, 1]
        pval_mat[j, k] <- res$coefficients[2, 4]

        all_pvals <- c(all_pvals, pval_mat[j, k])
      } # end celltype loop

      # plot vaf vs celltype diagnostic for grid
      if(vaf_cellprop){
        intercept=intercept_mat[j,]
        coef=coef_mat[j,]
        pval=pval_mat[j,]
        #print(pval)
        temp_plot_df = as.data.frame(cbind(spatial_coords_temp[c(test_ids,sample_idx_local),c("X","Y")],
                             vaf_j,N_j,Ws_j))
        colnames(temp_plot_df) = c("X","Y",var,paste0(var,"_cov"),celltypes)
        grid_vaf_cellprop_plot[[grid]] = plot_vaf_cellprop(j,
                             vaf[,c(test_ids,sample_idx_local),drop=F],Ws_j,temp_plot_df,
                            intercept,coef,pval,permute,return_plot=T)
      }
    } # end variant loop

    results_list[[spot_id]] <- list(intercept = intercept_mat,
                                    coef = coef_mat,
                                    pval = pval_mat)
  } # end spot loop

  plot_list  = list()
  # Generate Plots for Each celltype
  for (celltype in celltypes) {
    # For each spot (each element in results_list), extract the adjusted p-value and coefficient for the given variant and cell type.
    spot_values <- sapply(names(results_list), function(spot_id) {
      res <- results_list[[spot_id]]
      #if (is.null(res)) return(NA)
      #if (!(variant %in% rownames(res$adjusted_pval))) return(NA)
      #if (!(ct %in% colnames(res$adjusted_pval))) return(NA)

      p_val <- res$pval[, celltype]
      coef_val <- res$coef[, celltype]
      # Compute signed -log10(adjusted p)
      signed_logp <- sign(coef_val) * (-log10(p_val))
      return(signed_logp)

    })
    spatial_coords_temp$signed_logP = spot_values
    spatial_coords_temp$capped_signed_logP =  pmax(pmin(spot_values, max_log10p_cap),
                                                   -max_log10p_cap)
    # exclude spots with low diseased celltype proportions
    final_plot_temp = spatial_coords_temp
    final_plot_temp[exclude_plot_idx,"capped_signed_logP"] = 0

    p <- ggplot(final_plot_temp, aes(x = X, y = Y,color = capped_signed_logP)) +
      geom_point(size=0.1) +
      pval_palette +
      labs(title = paste(var, "-", celltype, "-", " P value"),
           x = "X", y = "Y", fill = "-log10(P-value)") +
      theme_minimal() +
      coord_fixed()  # Keep aspect ratio square

    # Store plot in list only if it has data
    plot_list[[paste(var, celltype, sep = "_")]] <- p
  }

  # Arrange Plots in a Grid
  combined_plot <- plot_grid(plotlist = plot_list, ncol = ceiling(sqrt(length(celltypes))))  # Adjust ncol if needed

  return(list(results=results_list,
              spatial_coords = spatial_coords_temp,
              combined_plot=combined_plot))
}



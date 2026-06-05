plot_empty <- function(title,msg=NULL) {
  if(!is.null(msg)){
    message(msg)
  }
  grid.newpage()
  grid.text(title, x = 0.5, y = 0.6, gp = gpar(fontsize = 14))
  grid.text("No genes selected", x = 0.5, y = 0.5, gp = gpar(fontsize = 11))
  
}

plot_DEG <- function(all_contrasts, # DEG results from SpotGLM
                     cell_types=NULL,
                     pdf_path=NULL,
                     effect_col="test_statistic",
                     top_k=20, # top genes to keep
                     q_thresh=0.1, # q-value for filtering DEGs
                     require_sig=TRUE, # TRUE = only keep genes with best q < q_thresh
                     use_q_col="qval",    
                     qcap_val = 10, # max capping value for plotting p-values
                     width=10,height=5){
  # Select top genes
  top_genes_by_ct <- all_contrasts %>%
    group_by(cell_type, gene) %>%
    summarise(
      best_q = suppressWarnings(min(qval, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(best_q = ifelse(is.infinite(best_q), NA_real_, best_q)) %>%  # convert Inf -> NA
    filter(!is.na(best_q)) %>%                                         # drop genes with no q
    { if (require_sig) filter(., best_q < q_thresh) else . } %>%
    arrange(cell_type, best_q) %>%
    group_by(cell_type) %>%
    slice_head(n = top_k) %>%
    summarise(genes = list(gene), .groups = "drop")
  
  
  # Plotting
  pdf(pdf_path, width = 10, height = 5)
  for (ct in cell_types) {
    print(paste0("Plotting ",ct, " ..."))
    genes_ct <- top_genes_by_ct %>% filter(cell_type == ct) %>% pull(genes)
    
    if(length(genes_ct) == 0){
      plot_empty(paste0("Cell type ", ct), msg = paste0(ct," - No genes pass q < ", q_thresh))
      next
    }else{
      genes_ct <- genes_ct[[1]]
    }
    
    # Effect matrix (e.g., estimate or test_statistic)
    df_eff <- all_contrasts %>%
      filter(cell_type == ct, gene %in% genes_ct) %>%
      select(gene, contrast, !!sym(effect_col)) %>%
      pivot_wider(names_from = contrast, values_from = !!sym(effect_col))
    
    mat_eff <- df_eff %>%
      column_to_rownames("gene") %>%
      as.matrix()
      
    # -log10(p) matrix
    df_p <- all_contrasts %>%
      filter(cell_type == ct, gene %in% genes_ct) %>%
      mutate(mlog10p = -log10(pmax(pval, .Machine$double.xmin))) %>%  # avoid Inf
      select(gene, contrast, mlog10p) %>%
      pivot_wider(names_from = contrast, values_from = mlog10p)
    
    mat_p <- df_p %>%
      column_to_rownames("gene") %>%
      as.matrix()
    
    # keep only contrasts that exist (drop all-NA columns)
    keep_cols <- intersect(
      colnames(mat_eff)[colSums(!is.na(mat_eff)) > 0],
      colnames(mat_p)[colSums(!is.na(mat_p)) > 0]
    )
    mat_eff <- mat_eff[, keep_cols, drop = FALSE]
    mat_p   <- mat_p[, keep_cols, drop = FALSE]
      
    # replace NAs so pheatmap doesn't crash
    mat_eff[is.na(mat_eff)] <- 0
    mat_p[is.na(mat_p)] <- 0
  
    can_cluster_rows <- (nrow(mat_eff) >= 2)
    can_cluster_cols <- (ncol(mat_eff) >= 2)
  
    # Make two pheatmaps and draw side-by-side
    ph1 <- pheatmap(
      mat_eff,
      #scale = "row",
      cluster_rows = can_cluster_rows,
      cluster_cols = can_cluster_cols,
      main = paste0(ct, " - ", effect_col),
      fontsize_row = 7,
      silent = TRUE
    )
    
    # P-value heatmap with same ordering as the effect size
    # cap -log10 p value by a max capping value.
    mat_p <- pmin(mat_p, qcap_val)
    ph2 <- pheatmap(
      mat_p,
      scale = "none",
      cluster_rows = if (can_cluster_rows) ph1$tree_row else FALSE,,
      cluster_cols = if (can_cluster_cols) ph1$tree_col else FALSE,,
      main = paste0(ct, " - -log10(p)"),
      fontsize_row = 7,
      silent = TRUE
    )
    
    #grid.newpage()
    grid.arrange(ph1$gtable, ph2$gtable, ncol = 2)
  }
  
  dev.off()
}

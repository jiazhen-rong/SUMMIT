#' Calculate Differentially Expressed Genes Between Two Groups of Cells
#'
#' This function performs differential expression analysis between two groups of cells,
#' defined by their cell barcodes. It creates a temporary metadata column in the Seurat object
#' to assign cells to either "Lineage1" or "Lineage2" based on the provided barcode lists,
#' subsets the object to only include those cells, and then runs Seurat's \code{FindMarkers()}
#' to identify differentially expressed genes.
#'
#' @param seu A Seurat object containing gene expression data.
#' @param lineage1_cells A character vector of cell barcodes belonging to the first group.
#' @param lineage2_cells A character vector of cell barcodes belonging to the second group.
#' @param assay A character string specifying which assay to use for DE analysis (default "RNA").
#' @param ... Additional arguments passed to \code{FindMarkers()} (e.g., test.use, logfc.threshold, etc.).
#'
#' @return A data frame containing the differentially expressed genes between the two groups.
#'
#' @examples
#' \dontrun{
#'   # Suppose you have two sets of cell barcodes from your AF.dm/N analysis:
#'   degs <- calculate_DEGs_from_cell_barcodes(seu, 
#'             lineage1_cells = c("Cell1", "Cell2", "Cell3"),
#'             lineage2_cells = c("Cell10", "Cell11", "Cell12"),
#'             assay = "RNA",
#'             test.use = "wilcox", logfc.threshold = 0.25)
#'   head(degs)
#' }
#'
#' @export
calculate_DEGs_from_cell_barcodes <- function(seu, lineage1_cells, lineage2_cells, assay = "RNA", ...) {
  # Create a temporary metadata column to hold lineage assignments
  all_cells <- colnames(seu)
  lineage_temp <- rep(NA, length(all_cells))
  names(lineage_temp) <- all_cells
  
  lineage_temp[all_cells %in% lineage1_cells] <- "Lineage1"
  lineage_temp[all_cells %in% lineage2_cells] <- "Lineage2"
  
  seu$lineage_temp <- lineage_temp
  
  # Subset Seurat object to only include cells that are assigned to a lineage
  # Identify the cells to keep (those with a non-NA assignment)
  cells_to_keep <- names(lineage_temp)[!is.na(lineage_temp)]
  # Subset the Seurat object by cell names
  seu_subset <- subset(seu, cells = cells_to_keep)
  
  # Set the identities to the temporary lineage assignment
  Idents(seu_subset) <- seu_subset$lineage_temp
  
  # Run differential expression analysis
  degs <- FindMarkers(seu_subset, ident.1 = "Lineage1", ident.2 = "Lineage2", assay = assay)

  return(degs)
}
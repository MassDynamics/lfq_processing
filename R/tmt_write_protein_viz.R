#' @import jsonlite
#' @export tmt_write_protein_viz


tmt_write_protein_viz <- function(prot, output_folder){
  protein_viz = list()
  comparisons <- get_comparisons(prot)
  i = 1
  for (comparison in comparisons){
    print(comparison)
    protein_viz[[i]] = list(
      "conditionComparison" = unbox(comparison),
      "fdrLimit" = unbox(0.05),
      "data" = prep_comparison_table(prot, comparison))
    i = 1 + i
  }

  write_json(protein_viz, file.path(output_folder, "protein_viz.json"))
}

prep_comparison_table <- function(prot, comparison){

  # disambiguate, assume no contaminants
  prot$Accession <- sapply(strsplit(prot$majority.protein.ids,";"), `[`, 1)
  prot$Gene <- get_genes_from_id(prot$Accession)
  prot$Accession_id <- str_extract(string = prot$Accession,
                                   pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")

  cols <- c("id", "Accession_id", "Gene", "fasta.headers", "q-value", "score",
            str_c(c("P.Value ", "adj.P.Val ", "logFC ", "CI.L ","CI.R "), comparison))

  comparison_dt = prot[,..cols]
  colnames(comparison_dt) <- c("ProteinGroupID", "ProteinId", "GeneName",
                               "FastaHeaders", "ProteinQValue", "ProteinScore",
                               "PValue", "AdjustedPValue", "FoldChange",
                               "ConfLow","ConfHigh")

  comparison_dt
}

get_comparisons <- function(prot){
  idx <- grep("logFC", colnames(prot))
  logFC_cols <- colnames(prot)[idx]
  comparisons <- gsub("logFC ","",logFC_cols)
  return(comparisons)
}

get_genes_from_id <- function(genes){
  genes <- gsub(".*\\|", "", genes)
  genes <- gsub("_.*", "", genes)
  genes
}

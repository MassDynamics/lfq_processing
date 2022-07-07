#' MaxQuant FQ protein-only QC plots names
#' @description QC report names of the plots created for a MaxQuant LFQ protein only analysis.

get_names_qc_lfq_protein_only <- function(){
  protein_only_split_qc <- c("PCA_proteins", "PCA_screeplot_proteins", "PCA_DE_proteins", 
                             "CV_proteins",
                             "samples_correlations_proteins",
                             "samples_correlations_DE_proteins", 
                             "samples_correlations_scatter_proteins",
                             "missing_by_proteins",
                             "missing_by_samples_proteins",
                             "identifications_proteins")
  return(protein_only_split_qc)
}

#' MaxQuant LFQ full QC plots names
#' @description QC report names of the plots created for a MaxQuant LFQ analysis including all files.

get_names_qc_lfq_all <- function(){
  all_lfq_split_qc <- c("PCA_proteins", "PCA_screeplot_proteins", "PCA_DE_proteins",
                             "PCA_peptides", "PCA_screeplot_peptides", "PCA_DE_peptides",
                             "PCA_modified_peptides", "PCA_screeplot_modified_peptides", "PCA_DE_modified_peptides",
                             "CV_proteins", "CV_peptides", "CV_modified_peptides",
                             "missing_by_proteins", "missing_by_peptides", "missing_by_modified_peptides",
                             "missing_by_samples_proteins", "missing_by_samples_peptides", "missing_by_samples_modified_peptides",
                             
                             "samples_correlations_proteins",
                             "samples_correlations_DE_proteins", 
                             "samples_correlations_scatter_proteins",
                             "samples_correlations_peptides",
                             "samples_correlations_DE_peptides", 
                             "samples_correlations_scatter_peptides",
                             "samples_correlations_modified_peptides",
                             "samples_correlations_DE_modified_peptides", 
                             "samples_correlations_scatter_modified_peptides",
                             
                             "identifications_proteins", "identifications_proteins_evidence",
                             "identifications_peptides", "identifications_peptides_evidence",
                             "identifications_modified_peptides", "identifications_modified_peptides_evidence",
                             
                             "missed_cleavages_evidence")
  return(all_lfq_split_qc)
}


#' MaxQuant TMT protein-only QC plots names
#' 
#' @param evidence bool. TRUE to return QC also available when the evidence.txt file is uploaded. 
#' @description QC report names of the plots created for a MaxQuant TMT protein-only analysis. 
#' If `evidence` file is provided, extra plots are produced. 

get_names_qc_tmt_protein_only <- function(evidence=FALSE){
  
  if(!evidence){
  tmt_split_qc_protein_only <- c("TMT_PCA_proteins", "TMT_PCA_screeplot_proteins", "TMT_PCA_DE_proteins", 
                        "TMT_CV_proteins",
                        "TMT_missing_by_proteins",
                        "TMT_missing_by_samples_proteins",
                        
                        "TMT_intensities_boxplots_proteins",
                        "TMT_intensities_normalised_boxplots_proteins",
                        
                        "TMT_imputed_proteins",
                        "TMT_identifications_proteins")
  } else {
    tmt_split_qc_protein_only <- c("TMT_parent_ion_frac_evidence",
                                   "TMT_missed_cleavages_evidence")
  }
                      
  return(tmt_split_qc_protein_only)
}
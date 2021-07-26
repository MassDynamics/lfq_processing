#' Run Pairwise Protein Quantification on each level of LFQ Maxquant output
#'
#' @param mq_folder A maxquant txt output folder.
#' @param output_folder An output folder to store produced files.
#' @param imputStDev The Standard Deviation parameter for MNAR Imputation
#' @param imputePosition The Position parameter for MNAR Imputation
#' @return A string describing the type of experiment
#' @examples
#' tmp =  lfq_transformer(mq_folder = upload_folder,
#'  output_folder = output_folder,
#'  imputeStDev=0.3,
#'  imputePosition=1.8)
#' # get each of the produced quantification files.
#' prot = tmp[[1]]
#' prot_int = tmp[[2]]
#' pept = tmp[[3]]
#' pept_int = tmp[[4]]
#' mod_pept = tmp[[5]]
#' mod_pept_int = tmp[[6]]
#' expdes = tmp[[7]]
#' evidence = tmp[[8]]
#' msms = tmp[[9]]
#'
#' @imports data.table
#' @export lfq_transformer
lfq_transformer <- function(mq_folder, output_folder, imputeStDev=0.3, imputePosition=1.8) {
  
  checkBy <- 'Experiment'
  
  #output_folder <- file.path(mq_folder, "transform")
  dir.create(output_folder, showWarnings = FALSE)
  
  
  # read data
  ma_tables <- lfq_read_data(mq_folder, experiment_type)
  msms <- ma_tables[[1]]
  prot <- ma_tables[[2]]
  pept <- ma_tables[[3]]
  mod_pept <- ma_tables[[4]]
  evidence <- ma_tables[[5]]
  expdes <- ma_tables[[6]]
  conditions_dict <- ma_tables[[7]]
  
  # Deal with data aggregated over fractions already
  if (length(unique(expdes$file_name))>length(unique(expdes$mqExperiment))){
    print("Detected Fractions. Dealing with post-frac-aggregate level data")
    expdes <- remove_fracs_files_from_des(expdes)
  }
  
  ###### ModPep ######
  mod_pep_list <- lfq_quant_analysis(mod_pept, expdes, id_var = "id", output_folder,
                                     "modificationSpecificPeptides_quant.txt",
                                     "modificationSpecificPeptides_quant_intensities.txt",
                                     conditions_dict=conditions_dict,
                                     imputeStDev = imputeStDev,
                                     imputePosition = imputePosition)
  mod_pept <- mod_pep_list[[1]]
  mod_pept_int <- mod_pep_list[[2]]
  
  ###### Peptides #######
  pept_list <- lfq_quant_analysis(pept, expdes, id_var = "id", output_folder,
                                  "peptides_quant.txt",
                                  "peptides_quant_intensities.txt",
                                  conditions_dict=conditions_dict,
                                  imputeStDev = imputeStDev,
                                  imputePosition = imputePosition)
  pept <- pept_list[[1]]
  pept_int <- pept_list[[2]]
  ###### Proteins ######
  prot_list <- lfq_quant_analysis(prot, expdes, id_var = "id", output_folder,
                                  "proteinGroups_quant.txt",
                                  "proteinGroups_quant_intensities.txt",
                                  conditions_dict=conditions_dict,
                                  imputeStDev = imputeStDev,
                                  imputePosition = imputePosition)
  prot <- prot_list[[1]]
  prot_int <- prot_list[[2]]
  conditionComparisonMapping <- prot_list[[3]]
  
  ###### phospho ######
  # phospho_enrichment_quant(upload_folder, expdes, output_folder)
  # id_table <- make_id_link_table(output_folder, pept)
  
  # CHange condition names in experiment design to match original names
  expdes <- condition_name_decode_intensity_data(dt = expdes, dict = conditions_dict, writerunid = FALSE)
  msms <- condition_name_decode_intensity_data(dt = msms, dict = conditions_dict, writerunid = FALSE)
  evidence <- condition_name_decode_intensity_data(dt = evidence, dict = conditions_dict, writerunid = FALSE)
  conditionComparisonMapping = decodeComparisonConditionMapping(conditionComparisonMapping, conditions_dict)
  
  return(list(prot, prot_int, pept, pept_int, mod_pept, mod_pept_int, expdes, evidence, msms, conditionComparisonMapping))
}


detect_fractions <- function(expdes){
  if (length(unique(expdes$file_name))>length(unique(expdes$mqExperiment))){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

remove_fracs_files_from_des <- function(expdes){
  expdes[,file_name := NULL]
  expdes <- unique(expdes[,list(mqExperiment, experiment, Replicate)])
}


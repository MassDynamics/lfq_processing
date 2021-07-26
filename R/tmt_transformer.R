#' Run Pairwise Protein Quantification on each level of TMT Maxquant output
#'
#' @param mq_folder A maxquant txt output folder.
#' @param output_folder An output folder to store produced files.
#' @param imputStDev The Standard Deviation parameter for MNAR Imputation
#' @param imputePosition The Position parameter for MNAR Imputation
#' @return A string describing the type of experiment
#' @examples
#' tmp =  tmt_transformer(mq_folder = upload_folder,
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
#' @export tmt_transformer
tmt_transformer <- function(mq_folder, output_folder, imputeStDev=0.3, imputePosition=1.8) {
  
  dir.create(output_folder, showWarnings = FALSE)
  
  protein_groups <- tmt_proteinGroup_txt_reader(mq_folder)
  
  des <- fread(file.path(mq_folder, "experimentDesign_original.txt"))
  des$reporter_channel <- as.character(des$reporter_channel)
  verify_tmt_des(des)
  
  protein_groups <- remove_not_annotated_channels(protein_groups, des)
  
  des_list <- condition_name_encoder_tmt(des = des)
  des <- des_list[[1]]
  conditions_dict <- des_list[[2]]
  
  tmp <- tmt_quant_analysis(protein_groups,
                            des,
                            "id",
                            imputeStDev=0.3,
                            imputePosition= 1.8)
  prot <- tmp[[1]]
  prot_int <- tmp[[2]]
  conditionComparisonMapping <- tmp[[3]]
  
  prot_int <- condition_name_decode_intensity_data_tmt(dt = prot_int, dict = conditions_dict)
  prot <- condition_name_decode_quantitative_data(dt = prot, dict = conditions_dict)
  des <- condition_name_decode_intensity_data_tmt(dt = des, dict = conditions_dict,  writerunid = FALSE)
  conditionComparisonMapping = decodeComparisonConditionMapping(conditionComparisonMapping, conditions_dict)
  
  write_output(prot, output_folder, "proteinGroups_quant.txt")
  write_output(prot_int, output_folder, "proteinGroups_quant_intensities.txt")
  
  return(list(prot, prot_int, des, conditionComparisonMapping))
}

# utilities
#' @export verify_tmt_des
verify_tmt_des <- function(des){
  stopifnot("reporter_channel" %in% tolower(colnames(des)))
  stopifnot("condition" %in% tolower(colnames(des)))
  stopifnot("replicate" %in% tolower(colnames(des)))
  stopifnot("experiment" %in% tolower(colnames(des)))
  
  if (detect_fractions(des)){
    stopifnot("mixture" %in% tolower(colnames(des)))
  }
}

#' @export remove_not_annotated_channels
remove_not_annotated_channels <- function(protein_groups, des){
  intensity_columns = colnames(protein_groups)[(grepl("reporter intensity corrected [0-9]* ", colnames(protein_groups)))]
  cols_to_keep <- str_c(str_c("reporter intensity corrected ", des$reporter_channel), tolower(des$experiment), sep = " ")
  cols_to_remove = setdiff(intensity_columns, cols_to_keep)
  
  for (col in cols_to_remove){
    protein_groups[,(col):=NULL]
  }
  
  
  return(protein_groups)
}

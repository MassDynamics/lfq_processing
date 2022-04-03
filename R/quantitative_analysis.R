#' Run Pairwise Protein Quantification on a level of LFQ Maxquant output
#'
#' @param dt A maxquant txt output text file (eg: proteinGroups.txt)
#' @param des An experiment design file.
#' @param id_var A column in dt identifying the unit (eg: majority.protein.id)
#' @param output_folder An output folder to store produced files
#' @param quant_fn A quantitative results file name
#' @param dt_int_fn An intensity results file name
#' @param conditions_dict A table with original and R-safe condition names
#' @param imputStDev The Standard Deviation parameter for MNAR Imputation
#' @param imputePosition The Position parameter for MNAR Imputation
#' @return a list containing the quantitative results first and intentisty results second.
#' @examples
#' \dontrun{
#' #producing results per modified peptide id
#' mod_pep_list <- lfq_quant_analysis(mod_pept, expdes, id_var = "id", output_folder,
#'  "modificationSpecificPeptides_quant.txt",
#'  "modificationSpecificPeptides_quant_intensities.txt",
#'  conditions_dict=conditions_dict,
#'  imputeStDev = imputeStDev,
#'  imputePosition = imputePosition)
#'
#' mod_pept <- mod_pep_list[[1]]
#' mod_pept_int <- mod_pep_list[[2]]
#' }
#' @export lfq_quant_analysis
lfq_quant_analysis <- function(dt, des, id_var, output_folder, quant_fn, dt_int_fn,
                               conditions_dict=conditions_dict,
                               imputeStDev=imputeStDev,
                               imputePosition=imputePosition) {
  # DEBUG
  # dt <- copy(prot)
  # des <- copy(expdes)
  # id_var = "id"
  # output_folder
  # quant_fn <- "proteinGroups_quant.txt"
  # dt_int_fn <- "proteinGroups_quant_intensities.txt"
  
  # quant_fn <- "modificationSpecificPeptides_quant.txt"
  # dt_int_fn <- "modificationSpecificPeptides_quant_intensities.txt"
  
  if (length(grep("^lfq intensity ", colnames(dt))) > 0) {
    dt_int <- melt.data.table(
      dt,
      id.vars = id_var,
      measure.vars = grep("lfq intensity ", colnames(dt), value = T),
      value.name = "intensity",
      variable.name = "mqExperiment"
    )
    dt_int[, mqExperiment := str_replace(mqExperiment, "lfq intensity ", "")]
    
  } else {
    dt_int <- melt.data.table(
      dt,
      id.vars = id_var,
      measure.vars = grep("intensity ", colnames(dt), value = T),
      value.name = "intensity",
      variable.name = "mqExperiment"
    )
    dt_int[, mqExperiment := str_replace(mqExperiment, "intensity ", "")]
  }
  
  dt_int <- merge(dt_int, des, by = "mqExperiment", all = T)
  dt_int[, Imputed := 0L]
  dt_int[intensity == 0, Imputed := 1L]
  dt_int[, log2NInt := 0.0]
  dt_int[intensity > 0 , log2NInt := log2(intensity)]
  
  # imputation
  dt_int <- impute_lfq(
    myQuantDT = dt_int,
    id_type = "id",
    int_type = "log2NInt",
    imputeStDev,
    imputePosition
  )
  dt_int[, run_id := str_c(experiment, Replicate, sep = ".")]
  
  
  results <- limma_stats_fun(
    ID_type = id_var,
    int_type = "log2NInt",
    condition_col_name = "experiment",
    run_id_col_name = "run_id",
    rep_col_name = "Replicate",
    funDT = dt_int
  )
  
  dt_quant <- results$stats 
  conditionComparisonMapping <- results$pairwise.comp
  
  dt_quant[, eval(id_var) := as.integer(get(id_var))]
  
  dt <- merge(dt, dt_quant, by = id_var, all.x = T)
  
  # recover original condition names
  dt <- condition_name_decode_quantitative_data(dt = dt, dict = conditions_dict)
  dt_int <- condition_name_decode_intensity_data(dt = dt_int, dict = conditions_dict)
  
  
  setnames(dt_int, "experiment", "condition")
  fwrite(dt, file.path(output_folder, quant_fn), row.names = FALSE, sep = "\t", quote = FALSE, showProgress = T, verbose = T, nThread = detectCores() - 1)
  fwrite(dt_int, file.path(output_folder, dt_int_fn), row.names = FALSE, sep = "\t", quote = FALSE, showProgress = T, verbose = T, nThread = detectCores() - 1)
  
  return(list(
    dt, dt_int, conditionComparisonMapping
  ))
}


#' @export tmt_quant_analysis
tmt_quant_analysis <- function(dt, des, id_var = "id",
                               imputeStDev=0.3,
                               imputePosition=1.8) {
  
  has_fractions <- detect_fractions(des)
  colnames(des) <- tolower(colnames(des))
  
  #colnames(dt) <- gsub(" ", ".", colnames(dt))
  measure_regex_channel = "reporter intensity corrected [0-9]*\\s\\s?"
  measure_vars = grep(measure_regex_channel, colnames(dt), value = T)
  
  if (length(measure_vars) == 0){
    print("Default backup: Try using single space between channel and experiment")
    measure_regex_channel = "reporter intensity corrected [0-9]* ?"
    measure_vars = grep(measure_regex_channel, colnames(dt), value = T)
  }
  
  # if measure_vars is empty, try changing the regex for not correct
  if (length(measure_vars) == 0){
    print("Default backup: Try using not corrected intensities")
    measure_regex_channel = "reporter intensity not corrected [0-9]*\\s\\s?"
    measure_vars = grep(measure_regex_channel, colnames(dt), value = T)
  }
  
  # if measure_vars is empty, try changing the regex for not correct
  if (length(measure_vars) == 0){
    print("Default backup: Try using not corrected intensities")
    measure_regex_channel = "reporter intensity not corrected [0-9]* ?"
    measure_vars = grep(measure_regex_channel, colnames(dt), value = T)
  }
  
  cat("==============\n")
  cat(paste("Number of identified measurement columns: ", length(measure_vars)),"\n")
  cat(paste("Number of rows in experimental design: "), nrow(des),"\n")
  
  stopifnot(length(measure_vars) == nrow(des)) # this has to be true for 
  # an experiment to work
  
  
  # melt data
  dt_int <- melt.data.table(
    dt,
    id.vars = id_var,
    measure.vars = measure_vars,
    value.name = "intensity",
    variable.name = "reporter_channel"
  )
  stopifnot("intensity" %in% colnames(dt_int))

  #print(grep(measure_regex_channel, colnames(dt), value = T))
  dt_int <- as.data.table(dt_int) # make sure it's a data table


  #get ready for imputation, log intensity column
  dt_int[, Imputed := 0L]
  dt_int[intensity == 0, Imputed := 1L]
  dt_int[, log2NInt := 0.0]
  dt_int[intensity > 0 , log2NInt := log2(intensity)]

  # extract info from columns
  dt_int[, experiment := mygsub("reporter intensity .?.?.? ?corrected [0-9]+? (.*)", reporter_channel)]
  dt_int[, experiment := tolower(experiment)]
  dt_int[, reporter_channel := mygsub("reporter intensity .?.?.? ?corrected ([0-9]+?) ?.*", reporter_channel)]

  # perform median subtraction for each reporter_channel group
  # reporter channel + condition level median subtraction
  dt_int[Imputed == F, ExpChannelMedian := median(log2NInt), by = list(reporter_channel,experiment)]
  dt_int[Imputed == F, log2NIntNorm := log2NInt - median(log2NInt), by = list(reporter_channel,experiment)]

  # imputation
  dt_int <- impute_lfq(
    myQuantDT = dt_int,
    id_type = id_var,
    int_type = "log2NIntNorm",
    f_imputeStDev = imputeStDev,
    f_imputePosition = imputePosition
  )
  
  
  # merge channel info with experiment design
  des$experiment <- tolower(des$experiment)
  dt_int <- merge(dt_int, des, by = c("experiment","reporter_channel"), all = T)
  
  # create run id's
  dt_int[, run_id := str_c(condition,experiment, replicate, sep = ".")]
  
  #call limma stats
  results <- limma_stats_fun(
    ID_type = id_var,
    int_type = "log2NIntNorm",
    condition_col_name = "condition",
    run_id_col_name = "run_id",
    rep_col_name = "replicate",
    funDT = dt_int
  )
  dt_quant <- results$stats 
  conditionComparisonMapping <- results$pairwise.comp
  
  dt$id <- as.character(dt$id)
  
  # check for duplicated column names
  dt <- merge(dt, dt_quant, by = id_var , all.x = T)
  
  
  return(list(dt, dt_int, conditionComparisonMapping))
}


get_channel_from_col_data <- function(data){
  tmp = gsub("reporter.intensity.corrected.","", data)
  channels = gsub(" .*","",tmp)
  return(channels)
}

get_experiment_from_col_data <- function(data){
  tmp = gsub("reporter.intensity.corrected.","", data)
  channels = gsub(".* ","",tmp)
  return(channels)
}

detect_fractions <- function(des){
  if (("fractions" %in% tolower(colnames(des)))){
    if (length(unique(des$fraction)) == 1){
      return(FALSE)
    } else {
      return(TRUE)
    }
  } else {
    return (FALSE)
  }
}


# return the regex capture if it is there, otherwise return an empty string
mygsub <- function( pattern, x){
  ans <- ifelse(grepl(pattern, x), 
                gsub(pattern, "\\1", x), 
                "")
  return(ans)
}


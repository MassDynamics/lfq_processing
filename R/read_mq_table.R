#' @import data.table
#' @param protein_only Boolean, TRUE means LFQ will only process protein level data
#' @export lfq_read_data
lfq_read_data <- function(upload_folder, experiment_type, protein_only) {
  
  # read experiment design
  expdes <- experiment_design_reader(upload_folder)
  expdes_list <- experiment_design_formatter(expdes)
  expdes <- expdes_list[[1]]
  conditions_dict <- expdes_list[[2]]
  
  cat('Reading Tables\n')
  
  
  # proteinGroups.txt
  prot <- generic_mq_table_reader(upload_folder, 'proteinGroups.txt')
  prot <- lfq_proteinGroup_txt_formatter(prot, expdes)
  
  if (!protein_only){ 
    # msms.txt
    msms <- msms_txt_reader(upload_folder, expdes)
    # peptides.txt
    pept <- lfq_peptides_txt_reader(upload_folder, expdes)
    # modificationSpecificPeptides.txt
    mod_pept <- lfq_modificationSpecificPeptides_txt_reader(upload_folder, expdes)
    # evidence.txt
    evidence <- evidence_txt_reader(upload_folder, expdes)
    
    return(list(
      msms = msms,
      prot = prot,
      pept = pept,
      mod_pept = mod_pept,
      evidence = evidence,
      expdes = expdes,
      conditions_dict = conditions_dict
    ))
  } else {
    return(list(
      prot = prot,
      expdes = expdes,
      conditions_dict = conditions_dict
    ))
  }
  


}

#' @export experiment_design_reader
experiment_design_reader <- function(folder) {
  expdes <- fread(file.path(folder, 'experimentDesign_original.txt'),
                  stringsAsFactors = F, header = T, verbose = F, 
                  keepLeadingZeros = TRUE)

  
  return(expdes)
}

#' @export experiment_design_formatter
experiment_design_formatter <- function(expdes) {
  setnames(expdes, "name", "file_name")
  expdes[, Replicate := bioRep]
  expdes[, `:=`(
    bioRep = NULL,
    techRep = NULL
  )]
  expdes[, file_name := str_replace(file_name, ".raw", "")]
  expdes[, mqExperiment := tolower(mqExperiment)]
  expdes[, experiment := as.character(experiment)]
  
  #remove rows with "TRUE" in remove column from experimental design.
  if ("remove" %in% colnames(expdes)){
    cat("remove TRUE remove")
    expdes$remove <- toupper(expdes$remove)
    nrows_removed = sum(as.character(expdes$remove) == "TRUE")
    expdes <- expdes[as.character(expdes$remove) != "TRUE"]
    cat(paste("Removed this many rows:" , nrows_removed))
  }
  
  expdes_list <- condition_name_encoder(des=expdes)
  expdes <- expdes_list[[1]]
  conditions_dict <- expdes_list[[2]]
  
  return(list(expdes, conditions_dict))
}


#' @export generic_mq_table_reader
generic_mq_table_reader <- function(folder, filename) {
  if (file.exists(
    file.path(folder, filename)
  )) {
    dt <- fread(file.path(folder, filename), stringsAsFactors = F, header = T, verbose = F, sep="\t")
    colnames(dt) <- tolower(colnames(dt))
    columns_to_filter <- c("reverse", "potential contaminant", "only identified by site", "contaminant")
    if (any(columns_to_filter %in% colnames(dt))) {
      columns_to_filter <- columns_to_filter[columns_to_filter %in% colnames(dt)]
      for (col in columns_to_filter) {
        dt <- dt[(get(col) != "+")|is.na(get(col))]
        dt <- dt[, (col) := NULL] #remove col
      }
    }
    return(dt)
  } else {
    return("File does not exist")
  }
  
}

#' @export get_lfq_intensity_columns
get_lfq_intensity_columns <- function(dt, des) {
  
  # if all lfq intensity columns present
  # and all intensity columns present
    # return those as well
  # if just intensity columns
  # return those
  # if just lfq intesnity columns
  # return those
  # else return na
  samples_measurement_cols_needed = unique(des$mqExperiment)
  
  
  # need to handle what happens when there's a column in protein groups we don't 
  # care about it being missing cause it got removed from teh experimental design. 
  
  if (all_columns_present(dt, samples_measurement_cols_needed, lfq = FALSE)){
    if (all_columns_present(dt, samples_measurement_cols_needed, lfq = TRUE)) {
      # got all the columns we need
      intensity_columns = str_c("intensity ",  des[, mqExperiment])
      lfq_intensity_columns = str_c("lfq intensity ",  des[, mqExperiment])
      intensity_columns = c(intensity_columns, lfq_intensity_columns)
    } else {
      # got just intensity columns
      intensity_columns = str_c("intensity ",  des[, mqExperiment])
    }
  } else { # try to save with lfq intensity columns
    if (all_columns_present(dt, samples_measurement_cols_needed, lfq = TRUE)) {
      # got all the columns we need
      intensity_columns = str_c("lfq intensity ",  des[, mqExperiment])
    }
  }
    
  if (exists("intensity_columns")) {
    return(unique(intensity_columns))
  } else {
    return(NA)
  }
}

all_columns_present <- function(dt, samples_measurement_cols_needed, lfq = F){
  if (lfq){
    return(
      all(samples_measurement_cols_needed %in% gsub("lfq intensity ", "", grep("^lfq intensity ", colnames(dt), value = T)))
    )
  } else {
    return(
      all(samples_measurement_cols_needed %in% gsub("intensity ", "", grep("^intensity ", colnames(dt), value = T)))
      )
  }
}

intensity_cols_to_double <- function(dt) {
  cols <- grep("intensity", colnames(dt), value = T)
  dt[, (cols) := lapply(.SD, as.double), .SDcols = cols]
  return(dt)
}

msms_txt_reader <- function(folder, des) {
  
  # DEBUG
  # folder = upload_folder
  # des = expdes
  
  msms <- generic_mq_table_reader(folder, 'msms.txt')
  msms$`raw file` <- as.character(msms$`raw file`)
  
  if (!is.null(des)){ #merge if you have one
    msms <- merge(msms, des, by.x = "raw file", by.y = "file_name", all = T)
  }
  msms <- intensity_cols_to_double(msms)
  return(msms)
}

#' @import dplyr
#' @export lfq_proteinGroup_txt_formatter
lfq_proteinGroup_txt_formatter <- function(dt, des) {

  intensity_columns = get_lfq_intensity_columns(dt, des)
  
  # filter rows by at least one non-zero value in intensity columns
  #cat("Removing rows with all zero values amongst used intensity columns")
  #cat(paste("There were ",sum(rowSums(dt[, ..intensity_columns])>0), "such rows"))
  #cat("----------------")
  #dt <- dt[rowSums(dt[, ..intensity_columns])>0,]
  if (!("score"  %in% colnames(dt))){
    dt$score <- NA
  }
  columns_whitelist <- c("protein ids", "majority protein ids", "q-value", "score", "peptide counts (all)", "peptide counts (razor+unique)", "peptide counts (unique)", "fasta headers", "number of proteins", "peptides",
                         "razor + unique peptides", "unique peptides", "sequence coverage [%]", "unique + razor sequence coverage [%]", "unique sequence coverage [%]", "mol. weight [kda]", "sequence length",
                         "id", "peptide ids", "evidence ids", "ms/ms ids", "best ms/ms", intensity_columns
  )
  
  dt <- dt[, columns_whitelist, with = FALSE]
  dt <- intensity_cols_to_double(dt)
  
  dt <- dt %>% mutate(
    `fasta headers` = case_when(
      `fasta headers`=="" ~ dt$`majority protein ids`,
      TRUE ~ dt$`fasta headers`)
  )
  
  return(dt)
}

lfq_peptides_txt_reader <- function(folder, des) {
  # DEBUG
  # folder = upload_folder
  # des = expdes
  
  dt <- generic_mq_table_reader(folder, 'peptides.txt')
  
  intensity_columns = get_lfq_intensity_columns(dt, des)
  
  # columns_whitelist <- c(
  #   "sequence", "amino acid before", "first amino acid", "second amino acid", "second last amino acid", "last amino acid", "amino acid after",
  #   "length", "missed cleavages", "mass", "proteins", "leading razor protein", "unique (groups)", "unique (proteins)", "charges", "pep", "score",
  #   "id", "protein group ids", "mod. peptide ids", "evidence ids", "ms/ms ids", "best ms/ms", "taxonomy ids",
  #   intensity_columns
  # )
  
  # filter columns
  columns_whitelist <- c(
    "sequence", "amino acid before", "amino acid after",
    "length", "missed cleavages", "mass", "proteins", "leading razor protein", "unique (groups)", "unique (proteins)", "charges", "pep", "score",
    "id", "protein group ids", "mod. peptide ids", "evidence ids", "ms/ms ids", "best ms/ms",
    intensity_columns
  )
  
  dt <- dt[, columns_whitelist, with = FALSE]
  dt <- intensity_cols_to_double(dt)
  return(dt)
}

lfq_modificationSpecificPeptides_txt_reader <- function(folder, des) {
  # DEBUG
  # folder = upload_folder
  # des = expdes
  
  dt <- generic_mq_table_reader(folder, 'modificationSpecificPeptides.txt')
  intensity_columns = get_lfq_intensity_columns(dt, des)
  
  columns_whitelist <- c(
    "sequence", "modifications", "mass", "protein groups", "proteins", "unique (groups)", "unique (proteins)",
    "missed cleavages", "retention time", "calibrated retention time", "charges", "pep", "ms/ms scan number", "score", "delta score",
    "id", "protein group ids", "peptide id", "evidence ids", "ms/ms ids", "best ms/ms",
    intensity_columns
  )
  
  dt <- dt[, columns_whitelist, with = FALSE]
  dt <- intensity_cols_to_double(dt)
  return(dt)
}

evidence_txt_reader <- function(folder, des) {
  # DEBUG
  # folder = upload_folder
  # des = expdes
  dt <- generic_mq_table_reader(folder, 'evidence.txt')
  columns_whitelist <- c(
    "sequence", "length", "modifications", "modified sequence","missed cleavages", "proteins", "leading proteins", "leading razor protein", "type", "raw file", "experiment",
    "charge", "m/z", "mass", "uncalibrated - calibrated m/z [ppm]", "uncalibrated - calibrated m/z [da]", "uncalibrated mass error [ppm]", "uncalibrated mass error [da]",
    "retention time", "retention length", "calibrated retention time", "calibrated retention time start", "calibrated retention time finish",
    "retention time calibration", "match time difference", "match m/z difference", "number of data points", "number of scans", "number of isotopic peaks",
    "pif", "fraction of total spectrum", "base peak fraction", "pep", "ms/ms count", "ms/ms scan number", "score", "delta score",
    "intensity", "id", "protein group ids", "peptide id", "mod. peptide id", "ms/ms ids", "best ms/ms"
  )
  
  dt <- dt[, columns_whitelist, with = FALSE]
  setnames(dt, "experiment", "mqExperiment")
  
  dt$`raw file` <- as.character(dt$`raw file`)
  dt <- merge(dt, des, by.x = "raw file", by.y = "file_name", all = T)
  
  dt <- intensity_cols_to_double(dt)
  return(dt)
}

phospho_sty_reader <- function(folder, des) {
  # DEBUG
  # folder = upload_foslder
  # des = expdes
  dt <- generic_mq_table_reader(folder, 'Phospho (STY)Sites.txt')
  columns_whitelist <- c(
    "sequence", "length", "modifications", "modified sequence","missed cleavages", "proteins", "leading proteins", "leading razor protein", "type", "raw file", "experiment",
    "charge", "m/z", "mass", "uncalibrated - calibrated m/z [ppm]", "uncalibrated - calibrated m/z [da]", "uncalibrated mass error [ppm]", "uncalibrated mass error [da]",
    "retention time", "retention length", "calibrated retention time", "calibrated retention time start", "calibrated retention time finish",
    "retention time calibration", "match time difference", "match m/z difference", "number of data points", "number of scans", "number of isotopic peaks",
    "pif", "fraction of total spectrum", "base peak fraction", "pep", "ms/ms count", "ms/ms scan number", "score", "delta score",
    "intensity", "id", "protein group ids", "peptide id", "mod. peptide id", "ms/ms ids", "best ms/ms"
  )
  
  dt <- dt[, columns_whitelist, with = FALSE]
  setnames(dt, "experiment", "mqExperiment")
  dt$`raw file` <- as.character(dt$`raw file`)
  dt <- merge(dt, des, by.x = "raw file", by.y = "file_name", all = T)
  
  dt <- intensity_cols_to_double(dt)
  return(dt)
}

# TMT Code

#' @export tmt_proteinGroup_txt_reader
tmt_proteinGroup_txt_reader <- function(upload_folder){
  
  proteinGroups <- generic_mq_table_reader(upload_folder, 'proteinGroups.txt')
  proteinGroups$id <- as.character(proteinGroups$id)
  return(proteinGroups)
}

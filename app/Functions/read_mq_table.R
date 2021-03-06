read_data <- function(upload_folder, expdes_filename) {
  
  # DEBUG
  # upload_folder = mq_folder
  # exp_des_filename = "experimentDesign.txt"
  
  # read experiment design
  expdes <- experiment_design_reader(upload_folder, expdes_filename)
  
  cat('Reading Tables\n')
  
  # msms.txt
  msms <- msms_txt_reader(upload_folder, expdes)
  # proteinGroups.txt
  prot <- proteinGroup_txt_reader(upload_folder, expdes)
  # peptides.txt
  pept <- peptides_txt_reader(upload_folder, expdes)
  # modificationSpecificPeptides.txt
  mod_pept <- modificationSpecificPeptides_txt_reader(upload_folder, expdes)
  # evidence.txt
  evidence <- evidence_txt_reader(upload_folder, expdes)
  
  return(list(
    msms,
    prot,
    pept,
    mod_pept,
    evidence,
    expdes
  ))
}


experiment_design_reader <- function(folder, exp_des_filename) {
  expdes <- fread(file.path(folder, exp_des_filename), stringsAsFactors = F, header = T, verbose = F)
  expdes[, Replicate := replicate]
  expdes[, `:=`(
    replicate = NULL
  )]
  expdes[, file_name := str_replace(file_name, ".raw", "")]
  expdes[, mqExperiment := tolower(mqExperiment)]
  return(expdes)
}

generic_mq_table_reader <- function(folder, filename) {
  if (file.exists(
    file.path(folder, filename)
  )) {
    dt <- fread(file.path(folder, filename), stringsAsFactors = F, header = T, verbose = F)
    colnames(dt) <- tolower(colnames(dt))
    columns_to_filter <- c("reverse", "potential contaminant", "only identified by site", "contaminant")
    if (any(columns_to_filter %in% colnames(dt))) {
      columns_to_filter <- columns_to_filter[columns_to_filter %in% colnames(dt)]
      for (col in columns_to_filter) {
        dt <- dt[get(col) != "+"]
        dt[, eval(col) := NULL]
      }
    }
    return(dt)
  }

}

get_intensity_columns <- function(dt, des) {
  if (any(grepl("^intensity ", colnames(dt)))) {
    intensity_columns = str_c("intensity ",  des[, mqExperiment])
  }
  
  if (any(grepl("^lfq intensity ", colnames(dt)))) {
    lfq_intensity_columns = str_c("lfq intensity ",  des[, mqExperiment])
    intensity_columns = c(intensity_columns, lfq_intensity_columns)
  } 
  if (exists("intensity_columns")) {
    return(intensity_columns)
  } else {
    return(NA)
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
  msms <- merge(msms, des, by.x = "raw file", by.y = "file_name", all = T)
  msms <- intensity_cols_to_double(msms)
  return(msms)
}

proteinGroup_txt_reader <- function(folder, des) {
  # DEBUG
  # folder = upload_folder
  # des = expdes
  
  dt <- generic_mq_table_reader(folder, 'proteinGroups.txt')
  intensity_columns = get_intensity_columns(dt, des)
  
  columns_whitelist <- c("protein ids", "majority protein ids", "q-value", "score", "peptide counts (all)", "peptide counts (razor+unique)", "peptide counts (unique)", "fasta headers", "number of proteins", "peptides",
                         "razor + unique peptides", "unique peptides", "sequence coverage [%]", "unique + razor sequence coverage [%]", "unique sequence coverage [%]", "mol. weight [kda]", "sequence length",
                         "id", "peptide ids", "evidence ids", "ms/ms ids", "best ms/ms", intensity_columns
                         )
  
  dt <- dt[, columns_whitelist, with = FALSE]
  dt <- intensity_cols_to_double(dt)

  return(dt)
}

peptides_txt_reader <- function(folder, des) {
  # DEBUG
  # folder = upload_folder
  # des = expdes
  
  dt <- generic_mq_table_reader(folder, 'peptides.txt')
  intensity_columns = get_intensity_columns(dt, des)

  # columns_whitelist <- c(
  #   "sequence", "amino acid before", "first amino acid", "second amino acid", "second last amino acid", "last amino acid", "amino acid after",
  #   "length", "missed cleavages", "mass", "proteins", "leading razor protein", "unique (groups)", "unique (proteins)", "charges", "pep", "score",
  #   "id", "protein group ids", "mod. peptide ids", "evidence ids", "ms/ms ids", "best ms/ms", "taxonomy ids", 
  #   intensity_columns
  # )
  
  columns_whitelist <- c(
    "sequence", "amino acid before", "amino acid after",
    "length", "missed cleavages", "mass", "proteins", "leading razor protein", "unique (groups)", "unique (proteins)", "charges", "pep", "score",
    "id", "protein group ids", "mod. peptide ids", "evidence ids", "ms/ms ids", "best ms/ms", "taxonomy ids", 
    intensity_columns
  )
  
  dt <- dt[, columns_whitelist, with = FALSE]
  dt <- intensity_cols_to_double(dt)
  return(dt)
}

modificationSpecificPeptides_txt_reader <- function(folder, des) {
  # DEBUG
  # folder = upload_folder
  # des = expdes
  
  dt <- generic_mq_table_reader(folder, 'modificationSpecificPeptides.txt')
  intensity_columns = get_intensity_columns(dt, des)
  
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
  dt <- merge(dt, des, by.x = "raw file", by.y = "file_name", all = T)
  
  dt <- intensity_cols_to_double(dt)
  return(dt)
}

table_quant_analysis <- function(dt, des, id_var, output_folder, quant_fn, dt_int_fn) {
  # DEBUG
  # dt <- copy(prot)
  # des <- copy(expdes)
  # id_var = "id"
  # output_folder
  # quant_fn <- "proteinGroups_quant.txt"
  # dt_int_fn <- "proteinGroups_quant_intensities.txt"
  
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
    int_type = "log2NInt"
  )
  dt_int[, run_id := str_c(experiment, Replicate, sep = ".")]
  dt_quant <- limma_stats_fun(
    ID_type = id_var,
    int_type = "log2NInt",
    condition_col_name = "experiment",
    run_id_col_name = "run_id",
    rep_col_name = "Replicate",
    funDT = dt_int
  )
  dt_quant[, eval(id_var) := as.integer(get(id_var))]
  
  dt <- merge(dt, dt_quant, by = id_var, all.x = T)
  
  setnames(dt_int, "experiment", "condition")
  fwrite(dt, file.path(output_folder, quant_fn), row.names = FALSE, sep = "\t", quote = FALSE, showProgress = T, verbose = T, nThread = detectCores() - 1)
  fwrite(dt_int, file.path(output_folder, dt_int_fn), row.names = FALSE, sep = "\t", quote = FALSE, showProgress = T, verbose = T, nThread = detectCores() - 1)
  
  return(list(
    dt, dt_int
  ))
}

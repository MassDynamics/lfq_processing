phospho_enrichment_quant <- function(upload_folder, des, output_folder) {
  
  
  ## === Phosphorylation enrichment experiment ======
  phos_enrich <- FALSE
  if (file.exists(
    file.path(upload_folder, 'Phospho (STY)Sites.txt')
  )) {
    phos_enrich <- TRUE
    phos <- fread(file.path(upload_folder, 'Phospho (STY)Sites.txt'), stringsAsFactors = F, header = T, verbose = F)
    colnames(phos) <- tolower(colnames(phos))
    
    # fix intensity column type
    cols <- grep("intensity", colnames(phos), value = T)
    phos[, (cols) := lapply(.SD, as.double), .SDcols = cols]
    
    phos <- phos[reverse != "+"]
    phos[, reverse := NULL]
    
    measure_vars <- grep("intensity .+___[1-3]", colnames(phos), value = T)
    #phos[, lapply(.SD, sum, na.rm=T), .SDcols = measure_vars]
    
    phos_int <- melt.data.table(
      phos,
      id.vars = "id",
      measure.vars = measure_vars,
      value.name = "intensity",
      variable.name = "mqExperiment"
    )
    
    phos_int[, mqExperiment := str_replace(mqExperiment, "intensity ", "")]
    phos_int[, id_phos := str_c(id, gsub(".*(___[1-3])", "\\1", mqExperiment))]
    phos_int[, mqExperiment := str_replace(mqExperiment, "___[1-3]", "")]
    phos_int <- merge(phos_int, sed, by = "mqExperiment", all.x = T)
    phos_int[, mqExperiment := NULL]
    
    phos_int[, Imputed := 0L]
    phos_int[intensity == 0, Imputed := 1L]
    
    phos_int[, log2NInt := 0.0]
    phos_int[intensity > 0 , log2NInt := log2(intensity)]
    
    phos_int_id_keep <- phos_int[, .(max_int = max(log2NInt)), by = .(id_phos)][max_int > 0, unique(id_phos)]
    phos_int <- phos_int[id_phos %in% phos_int_id_keep]
    
    
    phos_int <- impute_lfq(
      myQuantDT = phos_int,
      id_type = "id_phos",
      int_type = "log2NInt"
    )
    
    phos_int[, run_id := str_c(experiment, Replicate, sep = ".")]
    phos_int[, id_temp := .GRP, by = .(id_phos)]
    phos_int <- phos_int[order(id_phos, experiment, Replicate)]
    phos_quant <- limma_stats_fun(
      ID_type = "id_phos",
      int_type = "log2NInt",
      condition_col_name = "experiment",
      run_id_col_name = "run_id",
      rep_col_name = "Replicate",
      funDT = phos_int
    )
    
    phos_quant <- merge(phos_quant, unique(phos_int[, .(id, id_phos)]), by = "id_phos", all.x = T)
    setnames(phos_int, "experiment", "condition")
    
    fwrite(phos_quant, file.path(output_folder, 'Phospho (STY)Sites_quant.txt'), row.names = FALSE, sep = "\t", quote = FALSE, showProgress = T, verbose = T, nThread = detectCores() - 1)
    fwrite(phos_int, file.path(output_folder, 'Phospho (STY)Sites_quant_intensities.txt'), row.names = FALSE, sep = "\t", quote = FALSE, showProgress = T, verbose = T, nThread = detectCores() - 1)
    
  }
  
}
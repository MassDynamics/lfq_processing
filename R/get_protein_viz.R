#' Print a Protein Visualization JSON file.
#'
#' @param prot a processed proteinGroups table
#' @param prot_int a protein intensities table
#' @param output_folder a folder to write the json file into
#' @return Writes a json object to the output folder
#' @examples
#' get_protein_viz(prot, prot_int, "../datasets/iPRG2015/transform")
#' @import jsonlite
#' @export get_protein_viz
get_protein_viz <- function(prot, prot_int, output_folder, conditionComparisonMapping){
  
  # Fix TMT log2IntNorm Column name
  if ("log2NIntNorm" %in% colnames(prot_int)){
    prot_int$log2NInt <- prot_int$log2NIntNorm
  }
  
  ## === cleanUp info from fasta header ================
  cl <- makeCluster(detectCores() - 1, type = "FORK")
  
  # For dev on windows.
  #cl <- makeCluster(detectCores() - 1, type = "SOCK")
  
  registerDoParallel(cl)
  
  fasta_db <- prot[, .(`fasta headers`, id, `peptide counts (all)`, `protein ids`, `majority protein ids`)]
  multiple_prot_index <- fasta_db[, grepl(";", `fasta headers`)] | fasta_db[, `fasta headers` == ""]
  
  multiple_prot_dt <- fasta_db[multiple_prot_index]
  
  fasta_headers <- foreach(dt_line = seq(1, multiple_prot_dt[, .N])) %do% {
    dt <- multiple_prot_dt[dt_line]
    # Fasta headers
    dt_fasta_headers <- data.table(
      id = dt[, id],
      `fasta headers` = unlist(
        str_split(dt[, `fasta headers`], ";")
      ),
      `majority protein ids` = unlist(
        str_split(dt[, `majority protein ids`], ";")
      )
    )
    # Protein ids
    dt_protein_ids = data.table(
      `protein ids` =  unlist(
        str_split(dt[, `protein ids`], ";")
      ),
      `peptide counts (all)` =  as.integer(unlist(
        str_split(dt[, `peptide counts (all)`], ";")
      ))
    )
    dt_protein_ids <- dt_protein_ids[!grep("CON__", `protein ids`)]
    dt_protein_ids <- dt_protein_ids[order(-`peptide counts (all)`)][1]
    
    if (dt_fasta_headers[`majority protein ids` %in% dt_protein_ids[, `protein ids`] & `fasta headers` != '', .N] > 0) {
      dt_fasta_headers <- dt_fasta_headers[`majority protein ids` %in% dt_protein_ids[, `protein ids`] & `fasta headers` != '']
    } else {
      dt_fasta_headers[`fasta headers` == "", `fasta headers` := `majority protein ids`]
      dt_fasta_headers <- data.table(id=dt_fasta_headers[, id][1], `fasta headers`=dt_fasta_headers[, `fasta headers`][1])
      
    }
    
    dt_fasta_headers[, id := as.integer(id)]
    return(dt_fasta_headers[, .(id, `fasta headers`)])
    
  }
  fasta_headers <- rbindlist(fasta_headers)
  setnames(fasta_headers, "fasta headers", "fasta_headers")
  
  fasta_db <- fasta_db[!multiple_prot_index]
  fasta_db[, fasta_headers := `fasta headers`]
  
  fasta_db <- rbindlist(list(fasta_db[, .(id, fasta_headers)], fasta_headers), use.names = T)
  
  fasta_db[, ProteinId := ""]
  fasta_db[, gene_name := ""]
  fasta_db[, protein_description := ""]
  
  prot <- merge(prot, fasta_db[, .(id, ProteinId, gene_name, protein_description, fasta_headers)], by = "id", all.x = T)
  
  
  # Get name of columns and comparisons from protein table
  adj_pval_cols <- grep("adj.P.Val", colnames(prot), value = TRUE)
  adj_pval_cols <- adj_pval_cols[adj_pval_cols != "adj.P.Val"]
  pval_cols <- grep("P.Value", colnames(prot), value = TRUE)
  
  comparisons <- str_remove(pattern = "P.Value ", string = pval_cols)
  
  adj_test_pre <- "adj.P.Val " # unique(str_extract(pattern = "adj.P.Val ", string = pval_cols))
  test_pre <- "P.Value " # unique(str_extract(pattern = "adj.P.Val ", string = pval_cols))
  
  estimate_pre <- "logFC " # unique(str_extract(pattern = "logFC ", string = grep("logFC", colnames(prot), value = TRUE)))
  confLow_pre <- "CI.L " # unique(str_extract(pattern = "PH.conf.low - |t.test conf.low - ", string = grep("conf.low", colnames(prot), value = TRUE)))
  confHigh_pre <- "CI.R " # unique(str_extract(pattern = "PH.conf.high - |t.test conf.high - ", string = grep("conf.high", colnames(prot), value = TRUE)))
  
  FUN <- function(j, DT, comparison_pair) {
    #debug
    # j = 1
    # DT = copy(dt)
    # comparison_pair = comparison_pair
    
    ID_num <- DT[j, ProteinGroupId]
    return(
      list(
        ProteinGroupId = ID_num,
        ProteinId = DT[j, ProteinId],
        GeneName = DT[j, gene_name],
        ProteinDescription = DT[j, protein_description],
        FastaHeaders = DT[j, FastaHeaders],
        ProteinQValue = DT[j, ProteinQValue],
        ProteinScore = DT[j, ProteinScore],
        PValue = DT[j, PValue],
        AdjustedPValue = DT[j, FDR],
        FoldChange = DT[j, FoldChange],
        ConfLow = DT[j, ConfLow],
        ConfHigh = DT[j, ConfHigh]
      )
    )
  }
  
  to_our <- list()
  
  #i = 1
  for (i in 1:length(comparisons)) {
    
    pval_col = str_c(test_pre, comparisons[i])
    adj_pval_col = str_c(adj_test_pre, comparisons[i])
    estimate_col = str_c(estimate_pre, comparisons[i])
    confLow_col = str_c(confLow_pre, comparisons[i])
    confHigh_col = str_c(confHigh_pre, comparisons[i])
    
    dt <- copy(prot[, .(
      ProteinGroupId = id,
      ProteinId,
      FastaHeaders = fasta_headers,
      ProteinQValue = `q-value`,
      ProteinScore = score,
      PValue = get(pval_col),
      FDR = get(adj_pval_col),
      FoldChange = get(estimate_col),
      ConfLow = get(confLow_col),
      ConfHigh = get(confHigh_col),
      gene_name,
      protein_description
      
    )])
    
    
    pc_fdr_limit <- dt[FDR <= 0.05, max(PValue, na.rm = T)]
    
    condition1 = getUpCondition(conditionComparisonMapping, comparisons[i])
    condition2 = getDownCondition(conditionComparisonMapping, comparisons[i])
    pair_wise_data_list <- parLapply(cl, 1:dt[, .N], FUN, DT=dt, 
                                     comparison_pair=c(condition1, condition2))
    
    to_our[[i]] <- list(
      conditionComparison = comparisons[i],
      up.condition = condition1,
      down.condition = condition2,
      fdrLimit = pc_fdr_limit,
      data = pair_wise_data_list
    )
    
  }
  write_json(to_our, path = file.path(output_folder, "protein_viz.json"), pretty = FALSE, auto_unbox = TRUE, digits = NA)
  
  ## ==== Replicate count by protein and condition ============
  
  # Calculate %Replicate per sample in Protein Intensity
  if ("Replicate" %in% colnames(prot_int)){
    setnames(prot_int, "Replicate", "replicate")
  }
  
  num_eplicates_by_exp <- unique(prot_int[, .(condition, replicate)])[, .N, by = .(condition)]
  prot_int[, ProteinGroupId := id]
  replicate_count_pc <- prot_int[Imputed == 0, .(count = .N), by = .(ProteinGroupId, condition)]
  replicate_count_pc <- merge(replicate_count_pc, num_eplicates_by_exp, by = "condition", all = T)
  replicate_count_pc[, PrecentageOfReplicates := count/N]
  
  replicate_count_dt <- dcast.data.table(replicate_count_pc, ProteinGroupId ~ condition, value.var = "PrecentageOfReplicates", fill = 0)
  replicate_count_dt <- melt.data.table(replicate_count_dt, id.vars = "ProteinGroupId", variable.name = "condition", value.name = "PrecentageOfReplicates")
  replicate_count_dt <- merge(replicate_count_dt, replicate_count_pc[, .(condition, ProteinGroupId, count)], by = c("condition", "ProteinGroupId"))
  replicate_count_dt[is.na(count), count := 0]
  
  replicate_count_dt <- merge(replicate_count_dt, prot[, .(ProteinGroupId = id,  ProteinQValue = `q-value`,
                                                           ProteinScore = score)], by = "ProteinGroupId", all.x = T)
  replicate_count_dt[, condition := as.character(condition)]
  #replicate_count_dt[, ProteinId := as.integer(ProteinId)]
  replicate_count_dt <- merge(replicate_count_dt, fasta_db, by.x = "ProteinGroupId", by.y = "id", all.x = T)
  
  prot_int <- prot_int[order(ProteinGroupId, condition, replicate)]
  
  replicate_count_dt_out <- foreach(id_n = replicate_count_dt[, unique(ProteinGroupId)]) %dopar% {
    # id_n = 2
    # dtn <- copy(replicate_count_dt)
    # int_dt <- copy(prot_int)
    dtn <- replicate_count_dt[ProteinGroupId == id_n]
    
    
    get_intensity <- function(one_prot_condition_dt) {
      # one_prot_condition_dt <- prot_int[ProteinGroupId == id_n & condition == "UPS1"]
      int_list = foreach(line_i = 1:one_prot_condition_dt[, .N]) %do% {
        return(
          list(
            replicateNum = one_prot_condition_dt[line_i, replicate],
            centeredIntensity = one_prot_condition_dt[line_i, nRLE][[1]],
            z_norm = one_prot_condition_dt[line_i, z_norm][[1]],
            log2NInt_ProteinGroupId = one_prot_condition_dt[line_i, log2NInt],
            Imputed = one_prot_condition_dt[line_i, Imputed]
          )
        )
      }
      return(int_list)
    }
    
    condition_list = foreach (line_i = 1:dtn[, .N]) %do% {
      # print(line_i)
      # dtn[line_i]
      int_list <- get_intensity(
        prot_int[ProteinGroupId == dtn[line_i, ProteinGroupId] & condition == dtn[line_i, condition]]
      )
      return(
        list(
          name =  dtn[line_i, condition],
          precentageOfReplicates = dtn[line_i, PrecentageOfReplicates],
          numberOfReplicateCount = dtn[line_i, count],
          intensityValues = int_list
        )
      )
    }
    
    return(list(
      ProteinGroupId = id_n,
      ProteinId = dtn[, unique(ProteinId)],
      GeneName = dtn[, unique(gene_name)],
      ProteinDescription = dtn[, unique(protein_description)],
      FastaHeaders = dtn[, unique(fasta_headers)],
      ProteinQValue = dtn[, unique(ProteinQValue)],
      ProteinScore = dtn[, unique(ProteinScore)],
      conditions = condition_list
    ))
  }
  
  jsonlite::write_json(replicate_count_dt_out, path = file.path(output_folder, "protein_counts_and_intensity.json"),
                       pretty = FALSE, auto_unbox = TRUE, digits = NA)
  
}


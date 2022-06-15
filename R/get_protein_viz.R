#' Print a Protein Visualization JSON file.
#'
#' @param prot a processed proteinGroups table
#' @param prot_int a protein intensities table
#' @param output_folder a folder to write the json file into
#' @return Writes a json object to the output folder
#' @import jsonlite
#' @export get_protein_viz
get_protein_viz <- function(prot, prot_int, output_folder, conditionComparisonMapping){
  
  cat("Writing Protein Viz Output")
  start = Sys.time()
  write_protein_viz(prot, output_folder, conditionComparisonMapping)
  end = Sys.time()
  cat("\n")
  cat(end-start)
  cat("\n")
  
  cat("Writing Protein Counts and Intensities")
  start_time <- Sys.time()
  writeReplicateData(prot_int, prot, output_folder)
  end_time <- Sys.time()
  cat("\n")
  cat(end_time-start_time)
  cat("\n")
}

#' @param prot a processed proteinGroups table
#' @param prot_int a protein intensities table
#' @param output_folder a folder to write the json file into
write_protein_viz <- function(prot, output_folder, conditionComparisonMapping){
  protein_viz = list()
  comparisons = conditionComparisonMapping$comparison.string
  i=1
  for (comparison in comparisons){
    fdrLimit = max(prot[prot[[str_c("adj.P.Val ",comparison)]] <= 0.05,][[str_c("P.Value ",comparison)]], na.rm = T)
    
    protein_viz[[i]] = list(
      "conditionComparison" = unbox(comparison),
      "up.condition" = unbox(getUpCondition(conditionComparisonMapping, comparison)),
      "down.condition" = unbox(getDownCondition(conditionComparisonMapping, comparison)),
      "fdrLimit" = fdrLimit,
      "data" = prep_comparison_table(prot, comparison)
    )
    i = i+1
  }
  
  write_json(protein_viz, file.path(output_folder, "protein_viz.json"),
             pretty = FALSE, auto_unbox = TRUE, digits = NA, na = "string")
}

#' @export prep_comparison_table
prep_comparison_table <- function(prot, comparison){

  prot <- parse_id_columns(prot)
  
  cols <- c("id", "Accession_id", "Gene", "ProteinDescription","fasta headers", "q-value", "score",
            str_c(c("P.Value ", "adj.P.Val ", "logFC ", "CI.L ","CI.R "), comparison))
  
  comparison_dt = prot[,..cols]
  colnames(comparison_dt) <- c("ProteinGroupId", "ProteinId", "GeneName","ProteinDescription",
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


#' Write Protein counts and intensities json
#' @param prot a processed proteinGroups table
#' @param prot_int a protein intensities table
#' @param output_folder a folder to write the json file into
#' Columns required: `ProteinId`, `GeneName`, `Description`, `log2NInt`, `Condition`,
#'  `Replicate`, `Imputed`, `ProteinQValue`, `ProteinScore`,  `FastaHeaders`. 
#' @export oneProteinReplData
writeReplicateData <- function(prot_int, prot, outputFolder){
  
  if ("ProteinGroupId" %in% colnames(prot_int)){
    if (!("id" %in% colnames(prot_int))){
      setnames(prot_int, old = "ProteinGroupId", new = "id")
    }
  }
  
  if ("ProteinGroupId" %in% colnames(prot)){
    if (!("id" %in% colnames(prot))){
      setnames(prot, old = "ProteinGroupId", new = "id")
    }
  }
  
  stopifnot("id" %in% colnames(prot_int))
  stopifnot("id" %in% colnames(prot))
  
  prot_int = merge(prot_int, prot, by.x = "id", by.y = "id", all.x = T)

  prot_int = prot_int[order(prot_int$condition),]
  prot_int[, numberOfReplicateCount:= sum(Imputed==0), by = .(condition, id)]
  
  prot_int = prot_int[order(prot_int$id),]
  
  prot_int <- parse_id_columns(prot_int)

  prot_int$Condition = prot_int$condition
  prot_int$Replicate = prot_int$replicate
  
  # if lfq then replicate is the mass spec run (collapsed if fractions)
  if ("mqExperiment" %in% colnames(prot_int)){
    setnames(prot_int, old = "mqExperiment", new = "replicate")
  } else if ("reporter_channel" %in% colnames(prot_int)){ # then it's labelled and we should have an exp/channel combo
    prot_int$replicate = paste(prot_int$reporter_channel,prot_int$experiment)
  }

  setnames(prot_int, old = "id", new = "ProteinGroupId")
  setnames(prot_int, old = "condition", new = "Condition")
  setnames(prot_int, old = "Accession_id", new = "ProteinId")
  setnames(prot_int, old = "Gene", new = "GeneName")
  setnames(prot_int, old = "score", new = "ProteinScore")
  setnames(prot_int, old = "q-value", new = "ProteinQValue")
  setnames(prot_int, old = "fasta headers", new = "FastaHeaders")
  setnames(prot_int, old = "nRLE", new = "centeredIntensity")
  
  if ("log2NIntNorm" %in% colnames(prot_int)){
    setnames(prot_int, old = "log2NIntNorm", new = "log2NInt_ProteinGroupId")
  } else {
    setnames(prot_int, old = "log2NInt", new = "log2NInt_ProteinGroupId")
  }

  proteinSet <- unique(prot_int$ProteinGroupId)
  
  # work out how many cores to use
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    num_workers <- 2L
  } else {
    # use all cores in devtools::test()
    num_workers <- max(1,detectCores()-1)
  }

  protList <- mclapply(proteinSet, function(prot) oneProteinReplData(prot_int[ProteinGroupId == prot,]),
                       mc.cores = num_workers)
  protDF <- do.call(rbind, protList)
  
  dir.create(outputFolder, showWarnings = FALSE)
  outPath = file.path(outputFolder,"protein_counts_and_intensity.json")
  write_json(protDF, outPath, digits = NA, na = "null")
}


#' Transform data for one protein from a long format to the nested data structure needed for the Replicates tab.
#' @param oneProt data.table single protein information stored in long format. 
#' Columns required: `ProteinId`, `GeneName`, `Description`, `log2NInt`, `Condition`,
#'  `Replicate`, `Imputed`, `ProteinQValue`, `ProteinScore`,  `FastaHeaders`. 
#' @export oneProteinReplData
oneProteinReplData <- function(oneProt){
  infoProt <- unique(oneProt[,c("ProteinGroupId", "ProteinId","GeneName", "ProteinDescription", "FastaHeaders", "ProteinQValue","ProteinScore")])
  
  infoConds <- oneProt[, numberOfReplicateCount:= sum(Imputed==0), by = Condition ]
  infoConds <- infoConds[, precentageOfReplicates:= sum(Imputed==0)/length(replicate), by = Condition ]
  
  infoConds <- unique(infoConds[, c("Condition", "precentageOfReplicates","numberOfReplicateCount")])
  setnames(infoConds, old = "Condition", new = "name")
  
  conditions <- data.frame(matrix(NA, nrow = length(infoConds$name), ncol = 4))
  for(cond_idx in 1:length(infoConds$name)){
    cond <- infoConds$name[cond_idx]
    infoOneCond <- infoConds[name %in% cond, ]
    
    oneCondRepl <- data.table(oneProt)[Condition %in% cond, c("replicate","log2NInt_ProteinGroupId", "Imputed", "centeredIntensity", "z_norm")] 
    #oneCondRepl$replicateNum <- 1:nrow(oneCondRepl)
    oneCondRepl <- oneCondRepl[,c("replicate","centeredIntensity", "z_norm", "log2NInt_ProteinGroupId", "Imputed")]
    
    entryCond <- dplyr::tibble(infoOneCond, intensityValues=list(oneCondRepl))
    
    conditions[cond_idx, ] <- entryCond
  }
  
  colnames(conditions) <- c("name", "precentageOfReplicates", "numberOfReplicateCount", "intensityValues")
  # Combine with protein Infos
  oneProtNested <- dplyr::tibble(infoProt, conditions=list(conditions))
}
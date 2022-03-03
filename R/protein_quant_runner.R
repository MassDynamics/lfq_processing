#' Run Pairwise Protein Quantification and Quality Control on Maxquant output
#' @param folder A maxquant txt output folder.
#' @param output_folder An output folder to store produced files.
#' @param protein_only Boolean, TRUE means LFQ will only process protein level data
#' @return A string describing the type of experiment
#' @example protein_quant_runner("../data/iPRG2015/txt/", "../data/iPRG2015/txt/transform")
#' @import data.table
#' @export protein_quant_runner
protein_quant_runner <- function(upload_folder, output_folder, protein_only = FALSE) {
  
  dir.create(output_folder, showWarnings = FALSE)
  
  start_time <- Sys.time()
  
  print("Running MQ Transformer")
  experiment_type = detect_exp_type(upload_folder)
  print("Detected experiment type:")
  print(experiment_type)
  
  # hacky solution to avoid refactoring LFQ for now.
  if (experiment_type == "LFQ"){
    start_time <- Sys.time()
    
    ma_tables <- lfq_read_data(upload_folder, experiment_type, protein_only)
    
    tmp =  lfq_transformer(ma_tables = ma_tables,
                           output_folder = output_folder,
                           imputeStDev=0.3,
                           imputePosition=1.8,
                           protein_only = protein_only)
    
    end_time <- Sys.time()
    print(end_time - start_time)
    
    prot = tmp$prot
    prot_int = tmp$prot_int
    conditionComparisonMapping = tmp$conditionComparisonMapping
    expdes = tmp$expdes
    
    if (!protein_only){
      pept = tmp$pept
      pept_int = tmp$pept_int
      mod_pept = tmp$mod_pept
      mod_pept_int = tmp$mod_pept_int
      evidence = tmp$evidence
      msms = tmp$msms
    }
    
    rm(tmp)
    
    checkBy <- 'Experiment'
    
    # get qc report from package
    file.copy(from=system.file("rmd","QC_Report.Rmd", package = "LFQProcessing"),
              to=output_folder,
              overwrite = TRUE, recursive = TRUE,
              copy.mode = TRUE)
    
  } else if (experiment_type == "LABEL"){
    
    start_time <- Sys.time()
    
    protein_groups <- tmt_proteinGroup_txt_reader(upload_folder)
    des <- fread(file.path(upload_folder, "experimentDesign_original.txt"))
    
    if (all(is.na(des$experiment))){ #if there's no lcms-run-name  you're in trouble
      if (!any(duplicated(des$reporter_channel))){
        #if reporter channels are used once, you're ok
        # do a rescue
        cat("Single Run TMT experiment, adding an experiment name")
        
        print("Current column names:")
        print(colnames(protein_groups)[grep("reporter intensity corrected [0-9]*",colnames(protein_groups))])
        
        colnames(protein_groups)[grep("reporter intensity corrected [0-9]*",colnames(protein_groups))] = paste(
          colnames(protein_groups)[grep("reporter intensity corrected [0-9]*",colnames(protein_groups))], 
          "single_run"
        )
        
        # do this for not corrected as well
        colnames(protein_groups)[grep("reporter intensity not corrected [0-9]*",colnames(protein_groups))] = paste(
          colnames(protein_groups)[grep("reporter intensity not corrected [0-9]*",colnames(protein_groups))], 
          "single_run"
        )
        
        print("New column names:")
        print(colnames(protein_groups)[grep("reporter intensity corrected [0-9]*",colnames(protein_groups))])
        
        
        #update experiment
        des$experiment <- as.character(des$experiment)
        des$experiment = "single_run"
      } 
    }
    stopifnot(!(all(is.na(des$experiment))))
    
    des$reporter_channel <- as.character(des$reporter_channel)
    verify_tmt_des(des)
    
    tmp <- tmt_transformer(protein_groups,
                           des,
                           output_folder,
                           imputeStDev=0.3,
                           imputePosition=1.8)
    end_time <- Sys.time()
    print(end_time - start_time)
    
    prot = tmp[[1]]
    prot_int = tmp[[2]]
    des = tmp[[3]]
    conditionComparisonMapping = tmp[[4]]
    rm(tmp)
    
    # get qc report from package
    file.copy(from=system.file("rmd","tmt_qc_report.Rmd", package = "LFQProcessing"),
              to=file.path(output_folder,"QC_Report.Rmd"),
              overwrite = TRUE, recursive = TRUE,
              copy.mode = TRUE)
    
  } else {
    print(experiment_type)
  }
  
  # temporarily turn off qc report
  output_format = "html"
  rmarkdown::render(file.path(output_folder, "QC_Report.Rmd"))
  output_format = "pdf"
  rmarkdown::render(file.path(output_folder, "QC_Report.Rmd"), output_format="pdf_document")
  
  
  # clean up
  figs <- list.files(file.path(output_folder, "QC_Report_files/figure-html/"))
  dir.create(file.path(output_folder, "figure_html/"))
  
  for (fig in figs){
    # copy the figure files for some reason
    file.copy(from=file.path(output_folder, "QC_Report_files/figure-html/",fig),
              to=file.path(output_folder, "figure_html/"),
              overwrite = TRUE, recursive = TRUE,
              copy.mode = TRUE)
  }
  
 # unlink(file.path(output_folder, "QC_Report_files/"), recursive = TRUE)
  #unlink(file.path(output_folder, "qc_report_files"), recursive = TRUE)
  #unlink(file.path(output_folder, "QC_Report.Rmd"), recursive = TRUE)
  
  rm(pept)
  rm(pept_int)
  rm(mod_pept)
  rm(mod_pept_int)
  rm(expdes)
  rm(evidence)
  rm(msms)
  rm(cvdt)
  
  gc()
  
  save(prot, prot_int, conditionComparisonMapping, file = file.path(output_folder,"protein_data.rda"))
  
  #Write Protein Viz
  cat("Writing Protein Viz")
  start_time <- Sys.time()
  get_protein_viz(prot, prot_int, output_folder, conditionComparisonMapping)
  end_time <- Sys.time()
  cat(end_time - start_time)
  
  
}


#' Run Pairwise Protein Quantification and Quality Control on Maxquant output
#' @param folder A maxquant txt output folder.
#' @param output_folder An output folder to store produced files.
#' @return A string describing the type of experiment
#' @example
#' protein_quant_runner("../data/iPRG2015/txt/",
#'  "../data/iPRG2015/txt/transform")
#' @imports data.table
#' @export protein_quant_runner
protein_quant_runner <- function(upload_folder, output_folder) {
  
  dir.create(output_folder, showWarnings = FALSE)
  
  start_time <- Sys.time()
  
  print("Running MQ Transformer")
  experiment_type = detect_exp_type(upload_folder)
  print("Detected experiment type:")
  print(experiment_type)
  
  # hacky solution to avoid refactoring LFQ for now.
  if (experiment_type == "LFQ"){
    start_time <- Sys.time()
    
    ma_tables <- lfq_read_data(upload_folder, experiment_type)
    
    tmp =  lfq_transformer(ma_tables = ma_tables,
                           output_folder = output_folder,
                           imputeStDev=0.3,
                           imputePosition=1.8)
    
    end_time <- Sys.time()
    print(end_time - start_time)
    
    prot = tmp[[1]]
    prot_int = tmp[[2]]
    pept = tmp[[3]]
    pept_int = tmp[[4]]
    mod_pept = tmp[[5]]
    mod_pept_int = tmp[[6]]
    expdes = tmp[[7]]
    evidence = tmp[[8]]
    msms = tmp[[9]]
    conditionComparisonMapping = tmp[[10]]
    
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
  
  unlink(file.path(output_folder, "QC_Report_files/"), recursive = TRUE)
  unlink(file.path(output_folder, "qc_report_files"), recursive = TRUE)
  unlink(file.path(output_folder, "QC_Report.Rmd"), recursive = TRUE)
  
  rm(pept)
  rm(pept_int)
  rm(mod_pept)
  rm(mod_pept_int)
  rm(expdes)
  rm(evidence)
  rm(msms)
  rm(cvdt)
  
  gc()
  
  #Write Protein Viz
  print("Writing Protein Viz")
  start_time <- Sys.time()
  get_protein_viz(prot, prot_int, output_folder, conditionComparisonMapping)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
}


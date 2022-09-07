#' Run Pairwise Protein Quantification and Quality Control on Maxquant output
#' @param folder A maxquant txt output folder.
#' @param output_folder An output folder to store produced files.
#' @param protein_only Boolean, TRUE means LFQ will only process protein level data
#' @return A string describing the type of experiment
#' @import data.table
#' @export protein_quant_runner

protein_quant_runner <- function(upload_folder, output_folder, protein_only = FALSE, write_qc = TRUE) {
  
  dir.create(output_folder, showWarnings = FALSE)
  
  start_time <- Sys.time()
  
  print("Running MQ Transformer")
  print("getwd")
  print(upload_folder)
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
      
      # get qc report from package
      file.copy(from=system.file("rmd","QC_Report.Rmd", package = "LFQProcessing"),
                to=output_folder,
                overwrite = TRUE, recursive = TRUE,
                copy.mode = TRUE)
      
      print("\nCopy all separate QCs\n")
      all_qcs <- LFQProcessing:::get_names_qc_lfq_all()
      
      for(qc_name in all_qcs){
        qc_report_name <- paste0("QC_", qc_name, ".Rmd")
        cat(paste("Writing ", qc_report_name))
        file.copy(from=system.file("rmd", qc_report_name, package = "LFQProcessing"),
                  to=file.path(output_folder, qc_report_name),
                  overwrite = TRUE, recursive = TRUE,
                  copy.mode = TRUE)
        
        # I need to run it here as I need datasets loaded in the environment
        rmarkdown::render(file.path(output_folder, qc_report_name), 
                          params = list(output_figure = "figure_html_separate/"),
                          output_format = rmarkdown::html_document(
                            self_contained=FALSE, 
                            lib_dir="qc_report_files", 
                            theme = "united",
                            fig_caption = TRUE,
                            df_print = "paged"))
        print("done")
      }
      
    } else {
      # get qc report from package
      file.copy(from=system.file("rmd","QC_Report_protein_only.Rmd", package = "LFQProcessing"),
                to=file.path(output_folder,"QC_Report.Rmd"),
                overwrite = TRUE, copy.mode = TRUE)
      
      print("\nCopy protein only separate QCs\n")
      protein_only_split_qc <- LFQProcessing:::get_names_qc_lfq_protein_only()
      
      for(qc_name in protein_only_split_qc){
        qc_report_name <- paste0("QC_", qc_name, ".Rmd")
        print(paste0("RENDER:", qc_name))
        file.copy(from=system.file("rmd", qc_report_name, package = "LFQProcessing"),
                  to=file.path(output_folder,qc_report_name),
                  overwrite = TRUE, copy.mode = TRUE)
        
        # I need to run it here as I need datasets loaded in the environment
        rmarkdown::render(file.path(output_folder, qc_report_name), 
                          params = list(output_figure = "figure_html_separate/"),
                          output_format = rmarkdown::html_document(
                            self_contained=FALSE, 
                            lib_dir="qc_report_files")
        )
        
      }
      
      
    }
    
    # Cleanup the seraparet QCs
    #rmd_separate_qc_cleanup(output_folder)
    
    rm(tmp)
    
    checkBy <- 'Experiment'
    

    
  } else if (experiment_type == "LABEL"){
    
    start_time <- Sys.time()
    
    protein_groups <- tmt_proteinGroup_txt_reader(upload_folder)
    expdes <- fread(file.path(upload_folder, "experimentDesign_original.txt"))
    
    if (all(is.na(expdes$experiment))){ #if there's no lcms-run-name  you're in trouble
      if (!any(duplicated(expdes$reporter_channel))){
        #if reporter channels are used once, you're ok
        # do a rescue
        cat("Single Run TMT experiment. Experiment name is blank")
        des$experiment = ""
      } else {
        cat("Problem, no experiment or is.na experiment column without unique channels")
        stop()
      }
    }
    stopifnot(!(all(is.na(expdes$experiment))))
    
    expdes$reporter_channel <- as.character(expdes$reporter_channel)
    verify_tmt_des(expdes)
    
    tmp <- tmt_transformer(protein_groups,
                           expdes,
                           output_folder,
                           imputeStDev=0.3,
                           imputePosition=1.8)
    end_time <- Sys.time()
    print(end_time - start_time)
    
    prot = tmp[[1]]
    prot_int = tmp[[2]]
    expdes = tmp[[3]]
    conditionComparisonMapping = tmp[[4]]
    rm(tmp)
    
    # get qc report from package
    file.copy(from=system.file("rmd","tmt_qc_report.Rmd", package = "LFQProcessing"),
              to=file.path(output_folder,"QC_Report.Rmd"),
              overwrite = TRUE, recursive = TRUE,
              copy.mode = TRUE)
    
    print("\nCopy separate QCs\n")
    names_tmt_protein_only <- LFQProcessing:::get_names_qc_tmt_protein_only(evidence=FALSE)
    
    for(qc_name in names_tmt_protein_only){
      qc_report_name <- paste0("QC_", qc_name, ".Rmd")
      file.copy(from=system.file("rmd",qc_report_name, package = "LFQProcessing"),
                to=file.path(output_folder,qc_report_name),
                overwrite = TRUE, recursive = TRUE,
                copy.mode = TRUE)
      
      # I need to run it here as I need datasets loaded in the environment
      rmarkdown::render(file.path(output_folder, qc_report_name), 
                        params = list(output_figure = "figure_html_separate/"),
                        output_format = rmarkdown::html_document(
                          self_contained=FALSE, 
                          lib_dir=file.path(output_folder, "qc_report_files"))
      )
      
    }
    
    if (file.exists(file.path(upload_folder,"evidence.txt"))){
      evidence <- generic_mq_table_reader(upload_folder, 'evidence.txt')
      
      print("Evidence file was uploaded")

      names_tmt_evidence <- LFQProcessing:::get_names_qc_tmt_protein_only(evidence=TRUE)
      
      for(qc_name in names_tmt_evidence){
        qc_report_name <- paste0("QC_", qc_name, ".Rmd")
        file.copy(from=system.file("rmd",qc_report_name, package = "LFQProcessing"),
                  to=file.path(output_folder,qc_report_name),
                  overwrite = TRUE, copy.mode = TRUE)
        
        # I need to run it here as I need datasets loaded in the environment
        rmarkdown::render(file.path(output_folder, qc_report_name), 
                          params = list(output_figure = "figure_html_separate/"),
                          output_format = rmarkdown::html_document(
                            self_contained=FALSE,  
                            lib_dir="qc_report_files")
        )
      }
      
    } else {
      cat("Please upload your evidence.txt file to include Missed Cleavage/Parent Ion Fraction Statistics in this QC report.")
    }
    
  } else {
    print(experiment_type)
  }
  
  if (write_qc){
    output_format = "html"
    rmarkdown::render(file.path(output_folder, "QC_Report.Rmd"),
                      params = list(output_figure = "figure_html/"),
                      output_format = rmarkdown::html_document(
                        self_contained=TRUE,
                        code_folding= "hide",
                        theme="united",
                        toc = TRUE,
                        toc_float = TRUE,
                        fig_caption= TRUE,
                        df_print="paged",
                        lib_dir="qc_report_files"))

     output_format = "pdf"
     rmarkdown::render(file.path(output_folder, "QC_Report.Rmd"),
                       output_format=rmarkdown::pdf_document(
                         toc = TRUE,
                         fig_caption= TRUE),
                       params = list(output_figure = "figure_pdf/"))

  }

  # Remove folders
  list_lib_files <- list.files(output_folder, pattern = "qc_report_files", full.names = TRUE)
  unlink(list_lib_files, recursive = TRUE)
  list_rmd <- list.files(output_folder, pattern = "*Rmd$")
  unlink(file.path(output_folder, list_rmd), recursive = TRUE)
  
  fwrite(expdes, file.path(output_folder,"experimentDesign_output.txt"), sep ="\t")
  
  suppressWarnings(rm(pept))
  suppressWarnings(rm(pept_int))
  suppressWarnings(rm(mod_pept))
  suppressWarnings(rm(mod_pept_int))
  suppressWarnings(rm(expdes))
  suppressWarnings(rm(evidence))
  suppressWarnings(rm(msms))
  suppressWarnings(rm(cvdt))

  gc()
  
  cat("Saving protein data rda")
  save(prot, prot_int, conditionComparisonMapping, file = file.path(output_folder,"protein_data.rda"))
  
  #Write Protein Viz
  cat("Writing Protein Viz")
  start_time <- Sys.time()
  get_protein_viz(prot, prot_int, output_folder, conditionComparisonMapping)
  end_time <- Sys.time()
  cat(end_time - start_time)
  
  
}


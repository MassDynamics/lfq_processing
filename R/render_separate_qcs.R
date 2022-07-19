#' Render separate QC reports

render_separate_qcs <- function(experiment_type, 
                                protein_only, 
                                evidence_exists, 
                                output_folder){
  if (experiment_type == "LFQ"){
    
    if(!protein_only){
      
      all_qcs <- LFQProcessing:::get_names_qc_lfq_all()
      for(qc_name in all_qcs){
        qc_report_name <- paste0("QC_", qc_name, ".Rmd")
        rmd_separate_qc_call(output_folder, qc_report_name)
      }
      
    } else {
      
      protein_only_split_qc <- LFQProcessing:::get_names_qc_lfq_protein_only()
      
      for(qc_name in protein_only_split_qc){
        print(qc_name)
        qc_report_name <- paste0("QC_", qc_name, ".Rmd")
        rmd_separate_qc_call(output_folder, qc_report_name)
      }
    }
    
  } else if (experiment_type == "LABEL"){
    
    names_tmt_protein_only <- LFQProcessing:::get_names_qc_tmt_protein_only(evidence=FALSE)
    
    for(qc_name in names_tmt_protein_only){
      qc_report_name <- paste0("QC_", qc_name, ".Rmd")
      rmd_separate_qc_call(output_folder, qc_report_name)
    }
    
    if (evidence_exists){
      names_tmt_evidence <- LFQProcessing:::get_names_qc_tmt_protein_only(evidence=TRUE)
      
      for(qc_name in names_tmt_evidence){
        qc_report_name <- paste0("QC_", qc_name, ".Rmd")
        rmd_separate_qc_call(output_folder, qc_report_name)
      }
    }
  } else {
    print(paste0("QC reports not available for:", experiment_type))
  }
  
}



rmd_separate_qc_call <- function(output_folder, qc_report_name){
  rmarkdown::render(file.path(output_folder, qc_report_name), 
                    params = list(output_figure = file.path(output_folder, "figure_html_separate/")),
                    output_format = rmarkdown::html_document(
                      self_contained=FALSE,
                      lib_dir = file.path(output_folder,"qc_report_files"))
  )
}


rmd_separate_qc_cleanup <- function(output_folder){
  # clean up
  figs <- list.files(file.path(output_folder, "figure_html_separate/"))
  dir.create(file.path(output_folder, "figure_html_separate/"))
  
  for (fig in figs){
    # copy the figure files for some reason
    file.copy(from=file.path(output_folder, "figure_html_separate/",fig),
              to=file.path(output_folder, "figure_html_separate/"),
              overwrite = TRUE, recursive = TRUE,
              copy.mode = TRUE)
  }
}

#' @export clean_up_lfq_qc
clean_up_lfq_qc <- function(output_folder){

  # delete artifacts we don't need
  rm(pept)
  rm(pept_int)
  rm(mod_pep_list)
  rm(mod_pept_int)
  rm(cvdt)
  rm(pca_dt)
  rm(res.pca)
  rm(samples.pca)
  rm(screeplot)
  rm(eig.val)
  rm(DT_corMatrix)
  rm(dt)
  rm(g)
  rm(p)
  rm(eig.val)
  rm(int_corr_dt)
  rm(run_per_condition)
  rm(samples.coord)
  rm(scree_plot)

  # Organise output directory
  file.copy(
    file.path("inst","rmd", "QC_Report.html"),
    output_folder
  )

  dir.create(
    file.path(output_folder, "figure_html"),
    showWarnings = FALSE
  )

  dir.create(
    file.path(output_folder, "figure-html"),
    showWarnings = FALSE
  )

  file.copy(
    file.path(app_folder, "QC_Report_files", "figure-html"),
    file.path(output_folder),
    recursive = TRUE
  )

  file.copy(
    from = list.files(
      path = file.path(output_folder, "figure-html"),
      full.names = TRUE,
      recursive = TRUE
    ),
    to = file.path(output_folder, "figure_html"),
    overwrite = TRUE,
    recursive = FALSE,
    copy.mode = TRUE
  )

  unlink(
    file.path(output_folder, "figure-html"),
    recursive=TRUE
  )

}

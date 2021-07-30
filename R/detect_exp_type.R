#' Detect Experiment Type
#'
#' @param folder A maxquant txt output folder
#' @return A string describing the type of experiment
#' @export detect_exp_type
detect_exp_type <- function(folder){
  con <- file(file.path(folder, "proteinGroups.txt"),"r")
  first_line <- readLines(con,n=1)
  close(con)

  # for dev
  #print(first_line)

  if (grepl("LFQ intensity",first_line, fixed = TRUE)){
    return("LFQ")
  } else if (grepl("Reporter",first_line, fixed = TRUE)){
    return("LABEL")
  } else {
    return("Unknown - check experiment type manually.")
  }
}

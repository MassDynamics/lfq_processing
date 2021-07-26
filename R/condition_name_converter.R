#' @import data.table
#' @export
condition_name_encoder <- function(des=expdes) {
  # DEBUG
  # des = copy(expdes)


  conditions_dict <- data.table(original=des[, unique(experiment)])

  set.seed(255)
  conditions_dict[, safe := stri_rand_strings(.N, 5, pattern = "[A-Za-z]")]

  des[, original := experiment]
  des <- merge(des, conditions_dict, by="original", all = T)
  des[, experiment := safe]
  des[, `:=`(original = NULL, safe = NULL)]

  return(
    list(des, conditions_dict)
  )
}

#' @export
condition_name_encoder_tmt <- function(des=des) {
  # DEBUG
  # des = copy(expdes)


  conditions_dict <- data.table(original=des[, unique(condition)])

  conditions_dict[, safe := stri_rand_strings(.N, 5, pattern = "[A-Za-z]")]

  des[, original := condition]
  des <- merge(des, conditions_dict, by="original", all = T)
  des[, condition := safe]
  des[, `:=`(original = NULL, safe = NULL)]

  return(
    list(des, conditions_dict)
  )
}


#' @export
condition_name_decode_intensity_data <- function(dt, dict, writerunid=TRUE){
  # DEBUG
  #dt = copy(dt_int)
  #dict = copy(conditions_dict)

  dt <- merge(dt, dict, by.x = "experiment", by.y = "safe", all = T)
  dt[, experiment := original]
  dt[, original := NULL]
  if (writerunid) {
    dt[, run_id := str_c(experiment, Replicate, sep = ".")]
  }
  return(dt)
}

#' @export
condition_name_decode_intensity_data_tmt <- function(dt, dict, writerunid=TRUE){
  # DEBUG
  #dt = copy(dt_int)
  #dict = copy(conditions_dict)

  dt <- merge(dt, dict, by.x = "condition", by.y = "safe", all = T)
  dt[, condition := original]
  dt[, original := NULL]
  if (writerunid) {
    dt[, run_id := str_c(condition, replicate, sep = ".")]
  }
  return(dt)
}


#' @export
condition_name_decode_quantitative_data <- function(dt, dict){
  # DEBUG
  # dt = copy(pept)
  #dict = copy(conditions_dict)

  dt_columns <- colnames(dt)
  for (i in 1:dict[, .N]) {
    #print(i)
    #print(dict[i, ])
    original <- dict[i, original]
    safe <- dict[i, safe]

    dt_columns <- str_replace(dt_columns, pattern = safe, replacement = original)
  }
  colnames(dt) <- dt_columns
  return(dt)
}


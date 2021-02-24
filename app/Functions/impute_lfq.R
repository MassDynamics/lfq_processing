impute_lfq <- function(myQuantDT,
                       id_type,
                       int_type,
                       f_imputeStDev = imputeStDev,
                       f_imputePosition = imputePosition) {
  # DEBUG
  # myQuantDT <- copy(mod_pept_int)
  # f_imputeStDev <- imputeStDev
  # f_imputePosition <- imputePosition
  # id_type = "id"
  # int_type = "log2NInt"
  
  myQuantDT[, int_impute := get(int_type)]
  myQuantDT <- myQuantDT[, colnames(myQuantDT) != int_type, with = F]
  impute.StdDev <-
    f_imputeStDev * myQuantDT[Imputed == 0, sd(int_impute)]
  impute.position <-
    myQuantDT[Imputed == 0, mean(int_impute)] - f_imputePosition * myQuantDT[Imputed == 0, sd(int_impute)]
  set.seed(255)
  myQuantDT[Imputed == 1, int_impute := rnorm(.N, (impute.position), (impute.StdDev))]
  
  # Add columns for centered and centered + scaled intensities
  myQuantDT[, nRLE := scale(int_impute, center = TRUE, scale = FALSE), by = .(get(id_type))]
  myQuantDT[, z_norm := scale(int_impute, center = TRUE, scale = TRUE), by = .(get(id_type))]
  setnames(myQuantDT, "int_impute", "log2NInt")
  
  return(myQuantDT)
}
bread.lmerMod <- function (object){  
  ## get all elements by getME and exclude multiple random effect models.
  parts <- getME(object, "ALL")
  if (length(parts$l_i) > 1) stop("Multiple cluster variables detected. Robust SEs are unavailable.")
  
  vcov.full.lmerMod(object) * parts$l_i
}

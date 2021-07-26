# functions for dealing with condition names and comparison mapping coming out of limma
# hold over until we deal with this better in MassExpression

#' This is called in Limma to create a mapping for comparison strings to conditions
#' used later when writing output.
#' @export assembleComparisonConditionMapping
assembleComparisonConditionMapping <- function(conditionComparisonMapping, seperator = " - "){

  colnames(conditionComparisonMapping) = c("up.condition", "down.condition")
  conditionComparisonMapping$comparison.string = str_c(conditionComparisonMapping$up.condition,
                                          seperator,
                                          conditionComparisonMapping$down.condition)

  return(conditionComparisonMapping)
}


#' This function is a decoder utility that uses the conditions dict
#' @export getOriginalFromSafe
getOriginalFromSafe <- function(safe_word, conditions_dict){

  safe_word= as.character(safe_word)
  conditions_dict = as.data.frame(conditions_dict)

  match <- conditions_dict$original[conditions_dict$safe == safe_word]
  stopifnot(length(match) == 1)
  as.character(match)[1]
}

#' This function orchestrates the conversion of encoded strings
#' in the comparison condition mapping to the original conditions
#' @export decodeComparisonConditionMapping
decodeComparisonConditionMapping <- function(conditionComparisonMapping, conditions_dict){

  for (i in 1:(dim(conditionComparisonMapping)[1])){
    safe_up = conditionComparisonMapping[i, "up.condition"]
    safe_down = conditionComparisonMapping[i, "down.condition"]

    original_up = getOriginalFromSafe(safe_up, conditions_dict)
    original_down = getOriginalFromSafe(safe_down, conditions_dict)

    conditionComparisonMapping[i, "up.condition"] =  original_up
    conditionComparisonMapping[i, "down.condition"] = original_down

    safe_comparison = conditionComparisonMapping[i, "comparison.string"]
    # sthese gsubs are safe because of the no special char and non-overlap
    # properties of the safe words
    hybrid_comparison = gsub(safe_up, original_up, safe_comparison)
    original_comparison = gsub(safe_down, original_down, hybrid_comparison)
    conditionComparisonMapping[i, "comparison.string"] = original_comparison

  }

  conditionComparisonMapping
}

#' This function uses the condition comparison mapping to get the up condition
#' @export getUpCondition
getUpCondition <- function(conditionComparisonMapping, comparison.string){
  comparison.string <- as.character(comparison.string)
  relevant_comparison_index = conditionComparisonMapping$comparison.string == comparison.string
  stopifnot(sum(relevant_comparison_index)==1)
  conditionComparisonMapping$up.condition[relevant_comparison_index]
}

#' This function uses the condition comparison mapping to get the down condition
#' @export getDownCondition
getDownCondition <- function(conditionComparisonMapping, comparison.string){
  comparison.string <- as.character(comparison.string)
  relevant_comparison_index = conditionComparisonMapping$comparison.string == comparison.string
  stopifnot(sum(relevant_comparison_index)==1)
  conditionComparisonMapping$down.condition[relevant_comparison_index]
}


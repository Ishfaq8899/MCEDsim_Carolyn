#' Extract rate matrices and cancer sites based on OMST and LMST specifications.
#' 
#' @param the_omsts vector of overall mean sojourn times for each of the cancer specified cancer sites
#' @param the_lmsts vector of late mean sojourn times for each of the cancer specified cancer sites
#' @param all_meta_data dataframe that indexes the rate matrices corresponding to specified OMST and LMST entries
#' @param all_rates array of rate matrices corresponding to specified OMST and LMST entries
#' @return A list with two items: 
#'           1) A list of rate matrices that correspond to the selected OMST and LMST specifications for each the selected cancer sites. 
#'           2) vector of the cancer sites  
#' @export
#' @import purrr 
get_filtered_rates <- function(the_omsts, the_lmsts, all_meta_data, all_rates, the_cancer_sites) {
  the_indices <- all_meta_data %>%
    filter(OMST %in% the_omsts, LMST %in% the_lmsts, cancer_site %in% the_cancer_sites) %>%
    select("index")

  ###########################
  #rates_list <- lapply(all_fits[unlist(the_indices)], "[[", "rate.matrix")
  #Extract the rate matrices corresponding to the specified OMSTs and LMSTs 

  if(length(the_indices$index)>1){
  rates_list <- purrr::array_branch(all_rates[,,unlist(the_indices)],3)
  }else{
    rates_list <- list(all_rates[,,unlist(the_indices)])
    
  }

  ###########################
  
  cancer_sites <- all_meta_data %>%
    filter(OMST %in% the_omsts, LMST %in% the_lmsts, cancer_site %in% the_cancer_sites) %>%
    select("cancer_site")
  
  return(list(rates_list = rates_list, cancer_sites = cancer_sites$cancer_site))
}

#' Get the initial natural history state based on the rate matrix and starting age, conditional on no clinical diagnoses before starting age.
#' 
#' @param rate.matrix the rate matrix used to simulate cancer natural history
#' @param a1 starting age at first screen 
#' @return the initial state
get_init <- function(rate.matrix, a1) {
  k <- dim(rate.matrix)[1]
  init <- numeric(k)
  init[k] <- 0
  init[k - 1] <- 0
  prob_mat_a1 <- MatrixExp(t = a1, mat = rate.matrix)
  denom <- 1 - prob_mat_a1[1, k - 1] - prob_mat_a1[1, k]
  
  for (j in 1:(k - 2)) {
    init[j] <- prob_mat_a1[1, j] / denom
  }
  
  init_state <- sample(1:k, size = 1, prob = init)
  return(init_state)
}

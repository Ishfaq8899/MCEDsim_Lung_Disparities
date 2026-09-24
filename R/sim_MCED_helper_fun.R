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
#' @import dplyr
get_filtered_rates <- function(the_omsts, the_lmsts, all_meta_data, all_rates, the_cancer_sites) {

  target_combos <- tibble(
    OMST = the_omsts,
    LMST = the_lmsts,
    cancer_site = the_cancer_sites
  )

  #  browser()

  the_indices <- all_meta_data %>%
    semi_join(target_combos, by = c("OMST", "LMST", "cancer_site")) %>%
    select(index)

  #   browser()

  if(length(the_indices$index) > 1){
    rates_list <- purrr::array_branch(all_rates[,,unlist(the_indices)], 3)
  } else if(length(the_indices$index) == 1){
    rates_list <- list(all_rates[,,unlist(the_indices)])
  } else {
    rates_list <- list()
  }

  cancer_sites <- all_meta_data %>%
    semi_join(target_combos, by = c("OMST", "LMST", "cancer_site")) %>%
    distinct(index, cancer_site) %>%
    arrange(index) %>%
    select(cancer_site)

  #   browser()

  return(list(rates_list = rates_list, cancer_sites = cancer_sites$cancer_site))
}

#########################################
#' Get the initial natural history state based on the rate matrix and starting age, conditional on no clinical diagnoses before starting age.
#'
#' @param rate.matrix the rate matrix used to simulate cancer natural history
#' @param a1 starting age at first screen
#' @return the initial state
#' @import msm
#' @export
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

#' @export
match_individual<- function(i,  combined_additional_results, no_primary_cancer){

  out<-no_primary_cancer %>% filter(sex==combined_additional_results[i,"sex"]&age_OC_death_cat==combined_additional_results[i,"age_OC_death_cat"])

  if(nrow(out)>0){
    return(out[1,])
  }else{
    print("No match")
    return("no_match")
  }
}






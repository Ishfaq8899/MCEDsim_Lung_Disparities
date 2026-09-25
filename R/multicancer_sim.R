#############################################################################
# Function to simulate multiple cancer sites for an individual
#############################################################################
#' Simulate multiple cancer sites for an individual
#'
#' This function simulates the progression of multiple cancer sites for a single individual over a specified time period.
#'
#' @param ID                        Numeric identifier for the individual.
#' @param cancer_sites              Character vector specifying the names of cancer sites.
#' @param rate_matrices             List of rate matrices corresponding to each cancer site.
#' @param early_sensitivities       Numeric vector of early sensitivities for each cancer site.
#' @param late_sensitivities        Numeric vector of late sensitivities for each cancer site.
#' @param obs.times                 Numeric vector specifying observation times.
#' @param end.time                  Numeric value specifying the end time of the simulation.
#' @param start.time                Numeric value specifying the start time of the simulation.
#' @param start.states              Numeric value specifying the initial states for each of the cancer types.
#'
#' @return A data frame containing simulated observations for each cancer site.
#' @export
#'
#' @examples
#' sim_multiple_cancer_indiv(ID = 1, cancer_sites = c("Lung", "Liver"),
#'                            rate_matrices = list(lungrate, liverate),
#'                            early_sensitivities = c(0.3, 0.2),
#'                            late_sensitivities = c(0.9, 0.8),
#'                            obs.times = c(0, seq(50, 75, by = 2)), end.time = 100,
#'                            start.time = 0,
#'                            start.state = 1)
#'
sim_multiple_cancer_indiv <- function(ID,
                                      rate_matrices,
                                      early_sensitivities,
                                      late_sensitivities,
                                      specificities,
                                      start.time,
                                      end.time,
                                      obs.times,
                                      start.states,
                                      cancer_sites)  {

  emission_matrices <- mapply(FUN = "create_emission_matrix", rate_matrices,
                              early_sensitivities, late_sensitivities, specificities, SIMPLIFY = FALSE)

  get.obs.data.individual.temp <- function(rate.matrix, emission.matrix, start.state, cancer_site, ID,
                                           obs.times = obs.times, end.time = end.time, start.time = start.time) {
    get.obs.data.individual(ID = ID,
                            cancer_site = cancer_site,
                            rate.matrix = rate.matrix,
                            emission.matrix = emission.matrix,
                            obs.times = obs.times,
                            end.time = end.time,
                            start.time = start.time,
                            start.state = start.state)
  }

  out <- mapply(FUN = get.obs.data.individual.temp,
                rate_matrices,
                emission_matrices,
                start.states,
                cancer_sites,
                MoreArgs = list(ID = ID,
                                obs.times = obs.times,
                                end.time = end.time,
                                start.time = start.time),
                SIMPLIFY = FALSE)


  results <- do.call(rbind, out)
  results$cancer_site <- cancer_sites

  return(results)
}


############################################################################
# Function to simulate multiple cancer sites for the number of individuals
#############################################################################
#' Simulate multiple cancer sites for multiple individuals
#'
#' This function simulates the progression of multiple cancer sites for multiple individuals over a specified time period.
#'
#' @param num_individuals          Number of individuals to simulate.
#' @param cancer_sites             Character vector specifying the names of cancer sites.
#' @param rate_matrices            List of rate matrices corresponding to each cancer site.
#' @param early_sensitivities      Numeric vector of early sensitivities for each cancer site.
#' @param late_sensitivities       Numeric vector of late sensitivities for each cancer site.
#' @param specificities            Numeric vector of specificities for each cancer site.
#' @param obs.times                Numeric vector specifying observation times.
#' @param end.time                 Numeric value specifying the end time of the simulation.
#' @param start.time               Numeric value specifying the start time of the simulation.
#' @param start.state              Numeric value specifying the initial state of the simulation.
#'
#' @return A data frame containing simulated observations for each cancer site across multiple individuals.
#' @export
#'
#' @examples
#' sim_multiple_cancer_multiple_individuals(num_individuals = 50, cancer_sites = c("Lung", "Liver"),
#'                                           rate_matrices = list(lungrate, liverate),
#'                                           early_sensitivities = c(0.3, 0.2),
#'                                           late_sensitivities = c(0.9, 0.8),
#'                                            specificities = c(0.9, 0.85),
#'                                           obs.times = c(0, seq(50, 75, by = 2)), end.time = 100,
#'                                           start.time = 0,
#'                                           start.state = 1)
#'
sim_multiple_cancer_multiple_individuals <- function(num_individuals,
                                                     cancer_sites,
                                                     rate_matrices,
                                                     early_sensitivities,
                                                     late_sensitivities,
                                                     specificities,
                                                     obs.times,
                                                     end.time,
                                                     start.time,
                                                     start.state) {

  # obs.times = seq(1, 30, 2), end.time = 30,
  # start.time = 0, start.state = 1) {

  # browser()
  outlist = lapply(seq(1:num_individuals), FUN = function(i) sim_multiple_cancer_indiv(ID = i,
                                                                                       cancer_sites = cancer_sites,
                                                                                       rate_matrices = rate_matrices,
                                                                                       early_sensitivities = early_sensitivities,
                                                                                       late_sensitivities = late_sensitivities,
                                                                                       specificities = specificities,
                                                                                       obs.times = obs.times,
                                                                                       end.time = end.time,
                                                                                       start.time = start.time,
                                                                                       start.state = start.state))

  out = do.call(rbind, outlist)
  return(out)
}


############## Function to create emission matrix for a given rate matrix and sensitivities #####################
#' Create emission matrix for a given rate matrix and sensitivities
#'
#' This function creates an emission matrix for a given rate matrix and sensitivities.
#'
#' @param rate_matrix           Numeric matrix. Rate matrix for the cancer site.
#' @param early_sensitivity     Numeric. Early sensitivity for the cancer site.
#' @param late_sensitivity      Numeric. Late sensitivity for the cancer site.
#' @param specificities         Numeric. Specificities for each cancer site.
#'
#' @return A numeric matrix representing the emission matrix.
#'
#' @export
create_emission_matrix <- function(rate_matrix, early_sensitivity, late_sensitivity, specificity) {
  num_states <- nrow(rate_matrix)
  emission_matrix <- matrix(0, nrow = num_states, ncol = 6)

  # Define specific emission probabilities based on early and late sensitivities and specificity
  emission_matrix[1:(num_states - 4), 1] <- specificity
  emission_matrix[1:(num_states - 4), 6] <- 1 - specificity
  emission_matrix[(num_states - 3), 1] <- 1 - early_sensitivity
  emission_matrix[(num_states - 3), 2] <- early_sensitivity
  emission_matrix[(num_states - 2), 3] <- late_sensitivity
  emission_matrix[(num_states - 2), 1] <- 1 - late_sensitivity

  emission_matrix[(num_states - 1), 4] <- 1
  emission_matrix[num_states, 5] <- 1
  return(emission_matrix)
}


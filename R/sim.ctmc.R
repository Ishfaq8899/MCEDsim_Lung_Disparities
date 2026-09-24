########################################################################################
#Author: Jane Lange
#This function simulates from a time-homogeneous CTMC
#INPUTS: rate.matrix=rate matrix, start.state=starting state for CTMC, end.time=time
#         to stop data simulations; start.time= time to start data simulations
#OUTPUTS: a list with two objects: "times"=transition times and "states"=transition states
##########################################################################################
#' Simulate from a time-homogeneous continuous-time Markov chain (CTMC)
#'
#' This function simulates from a time-homogeneous CTMC characterized by
#'
#' @param start.state      Starting state for the CTMC.
#' @param rate.matrix      Rate matrix for the CTMC.
#' @param start.time       Time to start simulation (default is 0).
#' @param end.time         Time to stop simulations.
#' @param absorbing.states  Absorbing state .
#'
#' @return A list of two objects: "times" containing transition times and "states" containing transition states.
#'
#' @examples
#' # simulate from 2-state CTMC
#' rate.matrix <- matrix(c(-0.1, 0.1, 0.2, -0.2), nrow = 2)
#' sim_results <- sim.ctmc(start.state = 1, rate.matrix = rate_matrix, end.time = 10)
#'
#' @export
sim.ctmc <- function(start.state, rate.matrix, end.time, start.time = 0, absorbing.states = 0){

  state.space <- seq(1:dim(rate.matrix)[1])
  size <- dim(rate.matrix)[1]
  cur.state <- start.state
  times <- vector()
  states <- vector()

  times[1] <- start.time
  states[1] <- cur.state
  cur.time <- start.time

  k <- 2

  while ((cur.time < end.time) & !(cur.state %in% c(absorbing.states))) {
    exp.rate <- (-1) * rate.matrix[cur.state, cur.state]
    if (exp.rate == 0) {
      cur.time <- end.time
    } else {
      cur.time <- cur.time + rexp(n = 1, rate = exp.rate)
      if (cur.time < end.time) {
        times[k] <- cur.time

        if (size == 2){
          cur.state = as.numeric(state.space[-cur.state])
        } else {
          cur.state = sample(state.space[-cur.state], size = 1, prob = rate.matrix[cur.state, -cur.state])
        }

        states[k] <- cur.state
        k <- k + 1

      }
    }
  }

  return.list <- list(times, states)
  names(return.list) <- c("times", "states")
  return(return.list)
}

###################################################################################################
# This function gets the state of a CTMC at different discrete observation times
#INPUTS: ctmc.times=the transition times for the CTMC
#        ctmc.states=the states at each of the transition times
#        obs.times = the discrete observation times
#OUTPUTS: a dataframe with two columns: obs.times= observation times, states=value of state at obs.times
#WARNING: if the observation times are outside of max and min transition time, then
#         the state is assumed to be unchanged from the closest recorded transition time
#################################################################################################
#' Discretize a continuous-time Markov chain (CTMC) trajectory
#'
#' This function discretizes a CTMC trajectory given the transition times and states.
#'
#' @param ctmc.times    Transition times of the CTMC
#' @param ctmc.states   States of the CTMC at corresponding transition times.
#' @param obs.times.    Times at which observations are desired.
#'
#' @return A data frame containing the oberved times and corresponding states.
#'
#' @examples
#' # Discretize a CTMC trajectory
#' times <- c(0, 1, 2, 3, 4)
#' states <- c(1, 2, 1, 2, 1)
#' obs_times <- c(0.5, 1.5, 2.5, 3.5)
#' discretized_data <- discrete.ctmc(ctmc.times = times, ctmc.states = states, obs.times = obs_times)
#'
#' @export

discrete.ctmc<-function(ctmc.times, ctmc.states, obs.times){

  out<- data.frame(approx(x=ctmc.times, y=ctmc.states, xout=obs.times, rule=2, f=0, method="constant"))
  colnames(out)<-c("obs.times","states")
  return(out)
}


###################################################################################################
# This function gets an observed data point in a HMM based on an underlying state an emission matrix
#INPUTS: underlying.state = unobserved underlying state in HMM,
#        emmision.matrix=a matrix with the emission probablities.
#        the ith row corresponds to the hidden value X(t)=i, and the kth column to O(t)=k|X(t)=i
#        thus the rows sum to 1, and k columns correspond to the k possible observed states
#OUTPUTS: the observed data point
#################################################################################################
#' Get an observed data point from an underlying state using an emission matrix
#'
#' This function retrieves an observed data point from an underlying state using an emission matrix.
#'
#' @param underlying.state The underlying state from which observed data point is obtained.
#' @param emission.matrix The emission matrix representing the probablities of observing each state given an underlying state.
#'
#' @return An observed data point.
#'
#' @examples
#' #Get an observed data point
#' emission_matrix <- matrix(c(0.2, 0.5, 0.1, 0.4, 0.5), nrow = 2)
#' observed_data_point <- get.observed.datapoint(underlying.state = 1, emission.matrix = emission_matrix)
#'
#' @export

get.observed.datapoint<-function(underlying.state, emission.matrix){

  states<-seq(1:dim(emission.matrix)[2])
  probs<-emission.matrix[underlying.state,]
  sample(x=states, size=1, prob=probs)
}

###################################################################################################
# This function the observed states in an HMM for multiple observation times
# INPUTS: obs.times=a vector of observation times; underlying states=corresponding underlying states
#         emission.matrix=emission matrix for observed data
# OUTPUTS: dataframe with two columns: obs.times and obs.data (data at corresponding times)
#################################################################################################
#' Generate observed data from a Hidden Markov Model (HMM)
#'
#' This function generates observed data from a Hidden Markov Model (HMM) given the observation times, underlying states, and emission matrix.
#'
#' @param obs.times The observation times.
#' @param underlying.states The underlying states of the HMM.
#' @param emission.matrix The emission matrix representing the probabilities of observing each state given an underlying state.
#'
#' @return A data frame containing the observation times and corresponding observed data points.
#'
#' @examples
#' # Generate observed data from an HMM
#' obs_times <- c(0.5, 1.5, 2.5, 3.5)
#' underlying_states <- c(1, 2, 1, 2)
#' emission_matrix <- matrix(c(0.2, 0.3, 0.5, 0.1, 0.4, 0.5), nrow = 2)
#' observed_data <- observed.data.hmm(obs.times = obs_times, underlying.states = underlying_states, emission.matrix = emission_matrix)
#'
#' @export

observed.data.hmm <- function(obs.times, underlying.states, emission.matrix) {
  obs.data <- sapply(underlying.states, FUN = "get.observed.datapoint", emission.matrix)
  out <- data.frame(obs.times, obs.data)
  colnames(out) <- c("obs.times", "obs.data")
  return(out)
}

###################
# NEW: We added 4 cancer subtypes (Lung_Adeno, Lung_LargeCell,
#         - "Lung_SmallCell", "Lung_Squamous")
##################
# function for make_ctmc_seed (added 4 cancer subtypes)
make_ctmc_seed <- function(ID, cancer_site) {
  site_index <- match(cancer_site,c("Anus", "Bladder", "Esophagus", "Gastric", "Headandneck",
                                    "Liver", "Lung", "Lymphoma", "Ovary", "Pancreas", "Renal", "Uterine",
                                    "Lung_Adeno", "Lung_LargeCell", "Lung_SmallCell", "Lung_Squamous"))

  if (is.na(site_index)) {
    stop("Unknown cancer_site: ", cancer_site)
  }

  as.integer(ID * 100 + site_index)
}

#############################################################################################################
# This function obtains simulated data from an HMM at discrete observation times for a single individual
# INPUTS: rate.matrix=transition intensity matrix for underlying states, emission.matrix=emissin matrix for observed states
#         obs.times=discrete observation times; start.time=the start time for the CTMC simulation, start.state, the start state (single or vector)
#         for the CTMC, num.individuals=number of individuals to simulate data for
# OUTPUTS: dataframe with 5 columns: ID, screen_diagnosis_time, screen_diagnosis_stage, clinical_diagnosis_time, clinical_diagnosis_stage
##############################################################################################################
#' Get observed data for an individual from a CTMC trajectory and Hidden Markov Model (HMM)
#'
#' This function generates observed data for an individual from a continuous-time Markov chain (CTMC) trajectory and Hidden Markov Model (HMM).
#'
#' @param ID The identifier for the individual.
#' @param rate.matrix The rate matrix of the CTMC.
#' @param emission.matrix The emission matrix of the HMM.
#' @param obs.times The screening times (default is seq(1, 30, 2)).
#' @param end.time The end time for the CTMC trajectory (default is 30).
#' @param start.time The start time for the CTMC trajectory (default is 0).
#' @param start.state The start state for the CTMC trajectory (default is 1).
#'
#' @return A data frame containing the identifier, screen diagnosis time, screen diagnosis stage, clinical diagnosis time, and clinical diagnosis stage.
#'
#' @examples
#' # Get observed data for an individual
#' rate_matrix <- matrix(c(-0.1, 0.1, 0.2, -0.2), nrow = 2)
#' emission_matrix <- matrix(c(0.2, 0.3, 0.5, 0.1, 0.4, 0.5), nrow = 2)
#' obs_data_individual <- gets.obs.data.individual(ID = 1, rate.matrix = rate_matrix, emission.matrix = emission_matrix)
#'
#' @export
get.obs.data.individual <- function(ID, rate.matrix, emission.matrix,
                                    obs.times = seq(1, 30, 2), end.time = 30,
                                    start.time = 0, start.state = 1, cancer_site=NA) {

  n_states <- dim(rate.matrix)[1]
  clin_dx_late_state <- n_states
  clin_dx_early_state <- n_states - 1
  screen_early_state <- 2
  screen_late_state <- 3


  pre_clin_early_state <- n_states - 3
  pre_clin_late_state <- n_states -2

  #  set.seed(ID)
  if(!is.na(cancer_site)){
    # set random seed based on ID and cancer site
    the_seed <- make_ctmc_seed(ID, cancer_site)

    #  cat("SLURM_ARRAY_TASK_ID =", Sys.getenv("SLURM_ARRAY_TASK_ID"),
    #      "| ID =", ID,
    #      "| cancer_site =", cancer_site,
    #      "| seed =", the_seed, "\n")

    set.seed(the_seed)
  }

  #%%%%%%%%%%%%%

  # Simulate CTMC trajectory for individual
  trajectory <- sim.ctmc(rate.matrix = rate.matrix, start.state = start.state,
                         end.time = end.time, start.time = start.time,absorbing.states=c(n_states-1,n_states))

  onset_time=NA
  late_onset_time=NA

  if(pre_clin_early_state%in%trajectory$states){
    onset_time=min(trajectory$times[trajectory$states==pre_clin_early_state])
    # browser()
  }

  if(pre_clin_late_state%in%trajectory$states){
    late_onset_time=min(trajectory$times[trajectory$states==pre_clin_late_state])
  }

  # Discretize states at observation times
  discrete.states <- discrete.ctmc(ctmc.times = trajectory$times,
                                   ctmc.states = trajectory$states,
                                   obs.times = obs.times)

  # Generate observed data using HMM
  observed.data <- observed.data.hmm(obs.times = discrete.states$obs.times,
                                     underlying.states = discrete.states$states,
                                     emission.matrix = emission.matrix)

  # Clinical diagnosis time and stage
  clinical_diagnosis_time <- NA
  clinical_diagnosis_stage <- NA

  # Screen diagnosis time and stage
  screen_diagnosis_time <- NA
  screen_diagnosis_stage <- NA

  #Get clinical diagnosis time
  if (clin_dx_early_state %in% c(trajectory$states) || clin_dx_late_state %in% c(trajectory$states)) {

    clinical_diagnosis_index <- which(trajectory$states %in% c(clin_dx_early_state, clin_dx_late_state) == T)
    clinical_diagnosis_time <- trajectory$times[clinical_diagnosis_index]
    clinical_diagnosis_stage <- trajectory$states[clinical_diagnosis_index]

    clinical_diagnosis_stage <- ifelse(clinical_diagnosis_stage == clin_dx_early_state, "Early", "Late")
  }

  #Get screen diagnosis time
  if (screen_early_state %in% c(observed.data$obs.data) || screen_late_state %in% c(observed.data$obs.data)) {

    screen_diagnosis_index <- min(which(observed.data$obs.data %in% c(screen_early_state, screen_late_state) == T))
    screen_diagnosis_time <- observed.data$obs.times[screen_diagnosis_index]
    screen_diagnosis_stage <- observed.data$obs.data[screen_diagnosis_index]

    screen_diagnosis_stage <- ifelse(screen_diagnosis_stage == screen_early_state, "Early", "Late")
  }

  #Get cumulative number screens without cancer present
  total_no_canc_screens=sum(discrete.states$states<pre_clin_early_state)





  return(data.frame(ID, screen_diagnosis_time, screen_diagnosis_stage,
                    clinical_diagnosis_time, clinical_diagnosis_stage,onset_time,late_onset_time,total_no_canc_screens))
}


##############################################################################################################
#' Get observed data for an individual from of an individual in the control arm based on the CTMC natural history trajectory
#'
#' This function generates observed data for an individual from a continuous-time Markov chain (CTMC) trajectory and Hidden Markov Model (HMM).
#'
#' @param ID The identifier for the individual.
#' @param rate.matrix The rate matrix of the CTMC.
#' @param end.time The end time for the CTMC trajectory (default is 30).
#' @param start.time The start time for the CTMC trajectory (default is 0).
#' @param start.state The start state for the CTMC trajectory (default is 1).
#'
#' @return A data frame containing the identifier, screen diagnosis time, screen diagnosis stage, clinical diagnosis time, and clinical diagnosis stage.
#'
#' @examples
#' # Get observed data for an individual
#' rate_matrix <- matrix(c(-0.1, 0.1, 0.2, -0.2), nrow = 2)
#' obs_data_individual <- gets.obs.data.individual.control(ID = 1, rate.matrix = rate_matrix, end.time=30,start.time=0,start.state=1)
#'
#' @export
get.obs.data.individual.control <- function(ID, rate.matrix, end.time = 30,
                                            start.time = 0, start.state = 1) {

  n_states <- dim(rate.matrix)[1]
  clin_dx_late_state <- n_states
  clin_dx_early_state <- n_states - 1
  screen_early_state <- 2
  screen_late_state <- 3

  # Simulate CTMC trajectory for individual
  trajectory <- sim.ctmc(rate.matrix = rate.matrix, start.state = start.state,
                         end.time = end.time, start.time = start.time)


  # Clinical diagnosis time and stage
  clinical_diagnosis_time <- NA
  clinical_diagnosis_stage <- NA

  # Screen diagnosis time and stage
  screen_diagnosis_time <- NA
  screen_diagnosis_stage <- NA

  if (clin_dx_early_state %in% c(trajectory$states) || clin_dx_late_state %in% c(trajectory$states)) {

    clinical_diagnosis_index <- which(trajectory$states %in% c(clin_dx_early_state, clin_dx_late_state) == T)
    clinical_diagnosis_time <- trajectory$times[clinical_diagnosis_index]
    clinical_diagnosis_stage <- trajectory$states[clinical_diagnosis_index]

    clinical_diagnosis_stage <- ifelse(clinical_diagnosis_stage == clin_dx_early_state, "early", "late")
  }

  return(data.frame(ID, screen_diagnosis_time, screen_diagnosis_stage,
                    clinical_diagnosis_time, clinical_diagnosis_stage))
}

##########################################################################################################
# This function obtains simulated data from an HMM at discrete observation times for multiple subjects
# INPUTS: rate.matrix=transition intensity matrix for underlying states, emission.matrix=emissin matrix for observed states
#         obs.times=discrete observation times; start.time=the start time for the CTMC simulation, start.state, the start state (single or vector)
#         for the CTMC, num.individuals=number of individuals to simulate data for
# OUTPUTS: dataframe with 5 columns: ID, screen_diagnosis_time, screen_diagnosis_stage, clinical_diagnosis_time, clinical_diagnosis_stage
#########################################################################################################

#' Generate observed data for multiple individuals from a CTMC trajectory and Hidden Markov Model (HMM)
#'
#' This function generates observed data for multiple individuals from a continuous-time Markov chain (CTMC) trajectory and Hidden Markov Model (HMM).
#'
#' @param num.individuals The number of individuals for which to generate observed data.
#' @param rate.matrix The rate matrix of the CTMC.
#' @param emission.matrix The emission matrix of the HMM.
#' @param obs.times The observation times (default is seq(1, 30, 2)).
#' @param end.time The end time for the CTMC trajectory (default is 30).
#' @param start.time The start time for the CTMC trajectory (default is 0).
#' @param start.state The start state for the CTMC trajectory (default is 1).
#'
#' @return A data frame containing the observed data for multiple individuals.
#'
#' @examples
#' # Generate observed data for multiple individuals
#' rate_matrix <- matrix(c(-0.1, 0.1, 0.2, -0.2), nrow = 2)
#' emission_matrix <- matrix(c(0.2, 0.3, 0.5, 0.1, 0.4, 0.5), nrow = 2)
#' observed_data_many <- get.obs.data.many(num.individuals = 10, rate.matrix = rate_matrix, emission.matrix = emission_matrix)
#'
#' @export
get.obs.data.many<-function(num.individuals,rate.matrix, emission.matrix, obs.times = seq(1, 30, 2), end.time = 30,
                            start.time = 0, start.state = 1){
  outlist=lapply(seq(1:num.individuals),FUN="gets.obs.data.individual", rate.matrix=rate.matrix,
                 emission.matrix=emission.matrix, obs.times = obs.times, end.time = end.time, start.time = start.time,
                 start.state = start.state)

  out=do.call("rbind",outlist)
}






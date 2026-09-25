#####################################################
#' Simulate an individual with and without Multicancer Early Detection (MCED) screening
#'
#' This function simulates cancer outcomes for a single individual in two parallel scenarios:
#' with scheduled MCED screening and without screening. It generates cancer onset, screen detection,
#' clinical diagnosis, cancer death, and other-cause death, retaining the first cancer by earliest onset.
#'
#' Natural history models for each cancer site are specified based on a list of transition rate matrices for each site.
#' The function tracks only the first cancer diagnosis based on pre-clinical onset.
#' Cancer-specific mortality in the screen arm assumes a stage-shift benefit of screening. That is, individuals
#' who are diagnosed in early stage under screening but who would have been diagnosed in late stage clinically, are assume to remain in early stage post lead-time.
#' For both screened and unscreened individuals, cancer mortality is projected from the point of clinical diagnosis (i.e., post lead-time) to prevent
#' lead-time bias.
#'
#' @param ID Numeric ID for the individual.
#' @param cancer_sites Vector of cancer sites.
#' @param rates_list List of transition rate matrices for each cancer site.
#' @param test_performance A list with early_sens, late_sens, and specificities for the test.
#' @param MCED_specificity Overall specificity for the MCED test (not specific to cancer site)
#' @param other_cause_death_dist A table representing other-cause mortality.
#' @param starting_age Starting age for simulation.
#' @param num_screens Number of screening rounds.
#' @param screen_interval Interval between screening rounds.
#' @param end_time Ending time (age) of simulation.
#' @param sex "Male" or "Female".
#' @param race Single race value for this simulation run.
#' @param treatment_lookup Sex x Race x Histology x Stage table with
#'   novel_treatment_prob and novel_treatment_hr, pre-filtered to one race.
#' @param baseline_surv_table Unadjusted (HR=1) survival curve table (both
#'   Early and Late stage), used to compute novel-treatment death times live.
#' @param surv_param_table Survival parameters table (for simulating cancer death).
#' @param optimistic_surv_param_table Survival optimistic parameters table (for simulating cancer death).
#'
#' @export
#'
#' @return A data frame with the individual's first-cancer outcomes.
sim_individual_MCED<-function( ID,
                               cancer_sites,
                               rates_list,
                               test_performance,
                               MCED_specificity,
                               other_cause_death_dist,
                               starting_age,
                               num_screens,
                               screen_interval,
                               end_time,
                               sex,
                               race,
                               treatment_lookup,
                               baseline_surv_table,          #
                               adherence_lookup_table,       # NEW
                               surv_param_table,
                               optimistic_surv_param_table=NULL){

  set.seed(ID)
  # simulate time of other cause death
  other_cause_death = sim_othercause_death(other_cause_death_dist,ID=ID)

  # Get the starting states for each cancer
  start_states= unlist(lapply(rates_list, FUN="get_init",a1=starting_age))

  ### get the screening times
  screen_times = seq((starting_age),(num_screens + starting_age-1), by=screen_interval)

  # ==============================================================================
  # NEW: Screening adherence.
  # Everyone gets the first screen by default. Whether they get screens 2+
  # depends on a race-specific adherence probability (adherence_lookup).
  # ==============================================================================
  the_adherence_row <- adherence_lookup_table %>% filter(Race==race)

  if(nrow(the_adherence_row) != 1){
    stop("Expected exactly 1 matching row in adherence_lookup for Race=", race,
         " -- found ", nrow(the_adherence_row))
  }

  the_adherence_prob <- the_adherence_row$adherence_prob

  #  set.seed(ID * 1000 + 3)   # different from the above
  the_adherent <- rbinom(n = 1, size = 1, prob = the_adherence_prob)

  if(the_adherent == 0){
    # Not adherent: only the first scheduled screen actually happens.
    screen_times <- screen_times[1]
  }

  # =================================
  # END adherence block. screen_times
  # =================================

  #number of cancer sites
  num_sites=length(cancer_sites)

  #set cancer site specificity=1 because we are not tracking FP on a per cancer basis
  the_specificities=rep(1,times=num_sites)

  result <- sim_multiple_cancer_indiv(ID = ID,
                                      cancer_sites = cancer_sites,
                                      rate_matrices = rates_list,
                                      early_sensitivities = test_performance$early_sens,
                                      late_sensitivities =  test_performance$late_sens,
                                      specificities = the_specificities,
                                      obs.times = screen_times,
                                      start.time = starting_age,
                                      end.time = end_time,
                                      start.states =  start_states)

  # Add other-cause death info
  result$other_cause_death_status <- other_cause_death$status
  result$other_cause_death_time <- other_cause_death$time
  result$sex=sex
  result$adherent <- the_adherent          # NEW: record whether this person was adherent

  stored_result=result
  #-----------------------
  # Identify first cancer by onset time
  #-----------------------
  result$is_first_cancer=FALSE
  result$cancer_death_time_no_screen=NA
  result$cancer_death_time_screen=NA
  result$novel_tx_received=NA
  result$novel_tx_hr_used=NA

  if (sum(!is.na(result$onset_time)>=1)) {

    # selects the single row with the earliest (smallest) onset_time.
    first_cancer_row <- result %>%filter(!is.na(onset_time)) %>% slice_min(order_by = onset_time, n = 1, with_ties = FALSE)%>%mutate(is_first_cancer=TRUE)

    if(!is.null(optimistic_surv_param_table)){
      screen_surv_param_table=optimistic_surv_param_table
    }else{screen_surv_param_table=surv_param_table}


    if(!is.na(first_cancer_row$clinical_diagnosis_stage)){

      # ==============================
      # Treatment-assignment step
      # ============================
      the_tx_row <- treatment_lookup %>% filter(Sex == sex, Histology==first_cancer_row$cancer_site, Stage == first_cancer_row$clinical_diagnosis_stage)

      if(nrow(the_tx_row) != 1){stop("Expected exactly 1 matching row in treatment_lookup for Sex=", sex,
                                     ",Histology=", first_cancer_row$cancer_site,
                                     ", Stage=", first_cancer_row$clinical_diagnosis_stage,
                                     " -- found ", nrow(the_tx_row))
      }

      the_novel_tx_prob <- the_tx_row$novel_treatment_prob
      the_novel_tx_hr   <- the_tx_row$novel_treatment_hr

      set.seed(ID)
      the_novel_tx_received <- rbinom(n = 1, size = 1, prob = the_novel_tx_prob)

      first_cancer_row <- first_cancer_row %>% mutate(novel_tx_received = the_novel_tx_received, novel_tx_hr_used  = the_novel_tx_hr)

      # ==========================================================================
      # NEW: convert novel_tx_received into the_hr to pass into
      # sim_cancer_death_param(). If treated, pass their HR; otherwise NULL
      # (standard path). Replaces the old novel_tx_received/treatment_surv_table
      # pass-through.
      # ==========================================================================
      the_hr_to_use <- if(!is.na(first_cancer_row$novel_tx_received) && first_cancer_row$novel_tx_received == 1){
        first_cancer_row$novel_tx_hr_used
      } else {
        1
      }

      # NEW: overwrite novel_tx_hr_used with the value actually applied to
      # this person's death-time calculation (1 = no treatment, real HR =
      # treated), so the output data reflects reality instead of the raw
      # lookup-table value (which was NA for Early stage, and always 0.7
      # for anyone in a Late-stage group regardless of whether they were
      # actually selected by the Binomial draw).
      first_cancer_row <- first_cancer_row %>% mutate(novel_tx_hr_used = the_hr_to_use)



      first_cancer_row <-first_cancer_row  %>%mutate(cancer_death_time_no_screen=clinical_diagnosis_time+sim_cancer_death_param(the_stage=clinical_diagnosis_stage,
                                                                                                                                the_cancer_site=cancer_site,
                                                                                                                                the_sex=sex,
                                                                                                                                the_race=race,
                                                                                                                                the_model_type="Loglogistic",
                                                                                                                                param_table=surv_param_table,
                                                                                                                                the_hr=the_hr_to_use,
                                                                                                                                baseline_surv_table=baseline_surv_table,
                                                                                                                                ID=ID))


      #for non-optimistic scenario, do what we have previously done
      if(is.null(optimistic_surv_param_table)){
        first_cancer_row <-first_cancer_row  %>%
          mutate(cancer_death_time_screen=ifelse((screen_diagnosis_stage!=clinical_diagnosis_stage&!is.na(screen_diagnosis_stage))&!
                                                   is.na(clinical_diagnosis_stage),
                                                 clinical_diagnosis_time+
                                                   sim_cancer_death_param(the_stage="Early",
                                                                          the_cancer_site=cancer_site,
                                                                          the_sex=sex,
                                                                          the_race=race,
                                                                          the_model_type="Loglogistic",
                                                                          param_table=screen_surv_param_table,
                                                                          the_hr=the_hr_to_use,
                                                                          baseline_surv_table=baseline_surv_table,
                                                                          ID=ID),
                                                 cancer_death_time_no_screen))
        # This call now works safely even when the_hr_to_use is not NULL,
        # since baseline_surv_table has Early-stage rows too.

      }else{ #optimistic scenario

        first_cancer_row <-first_cancer_row  %>%
          mutate(cancer_death_time_screen=ifelse(!is.na(screen_diagnosis_stage)&!
                                                   is.na(clinical_diagnosis_stage),
                                                 clinical_diagnosis_time+
                                                   sim_cancer_death_param(the_stage=screen_diagnosis_stage,
                                                                          the_cancer_site=cancer_site,
                                                                          the_sex=sex,
                                                                          the_race=race,
                                                                          the_model_type="Loglogistic",
                                                                          param_table=screen_surv_param_table,
                                                                          the_hr=the_hr_to_use,
                                                                          baseline_surv_table=baseline_surv_table,
                                                                          ID=ID),
                                                 cancer_death_time_no_screen))

      }
    }
    result<-first_cancer_row

  }else{
    result<-slice_head(result,n=1)%>%mutate(cancer_site=NA)
  }

  stored_result <- stored_result %>% filter(!cancer_site==result$cancer_site)

  # set.seed(ID)
  result<-result %>% mutate(FP_tot=rbinom(n(),size=total_no_canc_screens,prob=1-MCED_specificity))

  return(list(first_result=result,stored_result=stored_result))
}


###########################################################
#' Simulate a cohort of individuals with and without Multicancer Early Detection (MCED) screening
#'
#' @description
#' This function simulates cancer outcomes in a population with a designated starting age using a "parallel universe" approach.
#' That is, cancer outcomes in the population are simulated with and without MCED screening.  The natural history (i.e., the times of cancer onset and clinical diagnosis) for
#' each individual is the same in both screening and no-screening scenarios.
#'
#' The user specifies cancer sites in the MCED screening test.
#' The user also provides sensitivity of the tests for early and late-stage disease, where early refers
#' to AJCC 7 stages I-II and late, III-IV, except for pancreas cancer, where early is stage I and late, II-IV.
#'
#' Natural history models are based on built-in fitted models that are calibrated to SEER 2015-2021 data by age and sex and can be specified based on
#' user-provided inputs about the overall mean sojourn time (OMST) and the late mean sojourn time (LMST) for each cancer site.
#' The function tracks only the first cancer diagnosis based on pre-clinical onset.
#'
#' Other-cause mortality is based on all cause mortality tables from the Human Mortality Database that have been adjusted to remove the mortality due to
#' cancers included in the MCED tests.  Cancer-specific mortality in the screen arm assumes a stage-shift benefit of screening. That is, individuals
#' who are diagnosed in early stage under screening but who would have been diagnosed in late stage clinically, are assume to remain in early stage post lead-time.
#' For both screened and unscreened individuals, cancer mortality is projected from the point of clinical diagnosis (i.e., post lead-time) to prevent
#' lead-time bias.
#'
#' @param cancer_sites Vector of cancer sites (allowable values include
#'   "Anus",  "Bladder", "Esophagus", "Gastric", "Headandneck",
#'  "Liver", "Lung", "Lymphoma", "Ovary", "Pancreas", "Renal", "Uterine")
#'
#' @param LMST_vec Numeric vector of late mean sojourn times (years) for each cancer site.
#' @param OMST_vec Numeric vector of overall mean sojourn times (years) for each cancer site.
#' @param test_performance_dataframe Data frame with test sensitivity/specificity info.
#' @param MCED_specificity Overall specificity for the MCED test (not specific to cancer site)
#' @param starting_age Numeric starting age for simulation
#' @param ending_age Numeric ending age for simulation
#' @param num_screens Number of screening rounds.
#' @param screen_interval Interval between screening rounds.
#' @param num_males Number of male individuals to simulate.
#' @param num_females Number of female individuals to simulate.
#' @param all_rates_male List of transition matrices for males.
#' @param all_rates_female List of transition matrices for females.
#' @param all_meta_data_female Metadata for female cancer sites.
#' @param all_meta_data_male Metadata for male cancer sites.
#' @param cdc_data CDC mortality data.
#' @param hmd_data Human Mortality Database data.
#' @param MCED_cdc CDC data for MCED.
#' @param surv_param_table Data frame of cancer-specific survival parameters.
#' @param simulation_seed Set seed for each simulation
#' @export
#'
#' @return A data frame with combined simulated results for all individuals.
#' The function returns each individual's first cancer site,  age and stage of clinical diagnosis in
#' absence of screening, age and stage at screen diagnosis, the time of other-cause mortality, and the time of cancer-specific mortality in absence
#' and presence of screening.   Cancer diagnosis and death times are presented both with and without competing other-cause mortality.
#'
#' @examples
#'library(LungDisparitiesMCEDsim)
#'
#'# Load the other-cause mortality tables
#'data("cdc_hmd_data")
#'# Load the prefitted natural history models
#'data("combined_fits")
#'# Load the prefitted cause-specific survival models
#'data("parametric_surv_fits")
#'#load("/home/groups/CEDAR/MCED_sim/parametric_surv_fits.rda")
#'
#'theseed      <- 1
#'scenario_no  <- 3
#'
#'cancer_sites_vec <- c(
#'  "Anus",  "Bladder", "Esophagus", "Gastric", "Headandneck",
#'  "Liver", "Lung", "Lymphoma", "Ovary", "Pancreas", "Renal", "Uterine")
#'
#'OMST_vec <- rep(2, 12)
#'LMST_vec <- rep(0.5, 12)
#'
#'early_sens <- c(
#'  0.5,  0.18, 0.48, 0.33, 0.72, 0.81,
#'  0.40, 0.61, 0.60, 0.61, 0.07, 0.18)
#'
#'late_sens <- c(
#'  1.00, 0.83, 0.97, 0.94, 0.93, 1.00,
#'  0.93, 0.94, 0.90, 0.94, 0.45, 0.81)
#'
#'test_performance_dataframe <- data.frame(early_sens  = early_sens,
#'                                         late_sens   = late_sens,
#'                                         cancer_site = cancer_sites_vec)
#'
#'set.seed(123)
#'results <- sim_multiple_individuals_MCED_parallel_universe(cancer_sites          = cancer_sites_vec,
#'                                                 LMST_vec              = LMST_vec,
#'                                                 OMST_vec              = OMST_vec,
#'                                                 test_performance_dataframe = test_performance_dataframe,
#'                                                 starting_age          = 45,
#'                                                 ending_age            = 500,
#'                                                 num_screens           = 30,
#'                                                 screen_interval       = 1,
#'                                                 num_males             = 50,
#'                                                 num_females           = 50,
#'                                                 all_rates_male        = all_rates_male,
#'                                                 all_rates_female      = all_rates_female,
#'                                                 all_meta_data_female  = all_meta_data_female,
#'                                                 all_meta_data_male    = all_meta_data_male,
#'                                                 cdc_data              = all_cause_cdc,
#'                                                 hmd_data              = hmd_data,
#'                                                 MCED_cdc              = MCED_cdc,
#'                                                 surv_param_table      = param_table,
#'                                                 MCED_specificity      = 0.995,
#'                                                 simulation_seed       = theseed)
sim_multiple_individuals_MCED_parallel_universe <- function(cancer_sites,
                                                              LMST_vec,
                                                              OMST_vec,
                                                              test_performance_dataframe,
                                                              MCED_specificity,
                                                              starting_age,
                                                              ending_age,
                                                              num_screens,
                                                              screen_interval,
                                                              num_males,
                                                              num_females,
                                                              all_rates_male,
                                                              all_rates_female,
                                                              all_meta_data_female,
                                                              all_meta_data_male,
                                                              cdc_data,
                                                              hmd_data,
                                                              MCED_cdc,
                                                              surv_param_table,
                                                              race,                          # single race for this run
                                                              treatment_lookup,              # Sex x Race x Histology x Stage novel-treatment table
                                                              baseline_surv_table,
                                                              adherence_lookup_table,        # NEW
                                                              optimistic_surv_param_table=NULL,
                                                              simulation_seed){



    total_individuals=num_males+num_females


    start_male=(simulation_seed-1)*(total_individuals)+1

    end_male=start_male+num_males-1
    # Create a vector of IDs
    IDs_male <- start_male:end_male

    # Female IDs: continue sequentially after males
    IDs_female <-(end_male+1):(num_females+end_male)


    # ---- Extract Sex-Specific Rate Matrices ----
    # Extract rate matrices matrices based on OMST and LMST specs (Male)
    rates_list_male = get_filtered_rates(the_omsts = OMST_vec, the_lmsts = LMST_vec,
                                         all_meta_data = all_meta_data_male,
                                         all_rates = all_rates_male, the_cancer_sites = cancer_sites)

    sites_male = rates_list_male$cancer_sites
    rates_list_male = rates_list_male$rates_list


    # Extract rate matrices matrices based on OMST and LMST specs (Female)
    rates_list_female = get_filtered_rates(the_omsts = OMST_vec, the_lmsts = LMST_vec,
                                           all_meta_data = all_meta_data_female,
                                           all_rates = all_rates_female, the_cancer_sites = cancer_sites)
    sites_female = rates_list_female$cancer_sites
    rates_list_female = rates_list_female$rates_list

    # ---- Extract Test Performance Parameters ----
    # Extract sensitivities and specificity based on selected cancer sites
    test_performance_male = test_performance_dataframe %>% filter(cancer_site %in% as.vector(sites_male))
    test_performance_female = test_performance_dataframe %>% filter(cancer_site %in% as.vector(sites_female))

    #Get the other-cause death tables for men and women
    other_cause_death_male=make_othercause_death_table(cdc_data=cdc_data,
                                                       MCED_cdc=MCED_cdc,
                                                       hmd_data=hmd_data,
                                                       the_starting_age = starting_age,
                                                       the_sex="Male",
                                                       selected_cancers=sites_male,
                                                       the_year=2018)

    other_cause_death_female=make_othercause_death_table(cdc_data=cdc_data,
                                                         MCED_cdc=MCED_cdc,
                                                         hmd_data=hmd_data,
                                                         the_starting_age = starting_age,
                                                         the_sex="Female",
                                                         selected_cancers=sites_female,
                                                         the_year=2018)

    # ---- Simulate Individual Outcomes ----
    # Use mapply to apply the sim_individual_MCED function to each ID (males)
    results_list_male <- mapply(sim_individual_MCED,
                                ID = IDs_male,
                                MoreArgs = list(rates_list=rates_list_male,
                                                cancer_sites=sites_male,
                                                test_performance=test_performance_male,
                                                other_cause_death_dist=other_cause_death_male,
                                                starting_age=starting_age,
                                                num_screens=num_screens,
                                                screen_interval=screen_interval,
                                                end_time=ending_age,
                                                surv_param_table=surv_param_table,
                                                optimistic_surv_param_table=optimistic_surv_param_table,
                                                race=race,                        # NEW
                                                treatment_lookup=treatment_lookup, # NEW
                                                baseline_surv_table=baseline_surv_table,   # RENAMED
                                                adherence_lookup_table=adherence_lookup_table,   # NEW
                                                sex="Male",MCED_specificity=MCED_specificity),
                                SIMPLIFY = FALSE)

    # Use mapply to apply the sim_individual_MCED function to each ID (females)
    results_list_female <- mapply(sim_individual_MCED,
                                  ID = IDs_female,
                                  MoreArgs = list(rates_list=rates_list_female,
                                                  cancer_sites=sites_female,
                                                  test_performance=test_performance_female,
                                                  other_cause_death_dist=other_cause_death_female,
                                                  starting_age=starting_age,
                                                  num_screens=num_screens,
                                                  screen_interval=screen_interval,
                                                  end_time=ending_age,
                                                  surv_param_table=surv_param_table,
                                                  optimistic_surv_param_table=optimistic_surv_param_table,
                                                  race=race,                        # NEW
                                                  treatment_lookup=treatment_lookup, # NEW
                                                  baseline_surv_table=baseline_surv_table,   # RENAMED
                                                  adherence_lookup_table=adherence_lookup_table,   # NEW
                                                  sex="Female",
                                                  MCED_specificity=MCED_specificity),
                                  SIMPLIFY = FALSE)


    #Get the first cancer and additional cancers for all individuals (female)
    first_site_female=lapply(results_list_female,"[[","first_result")
    additional_sites_female=lapply(results_list_female,"[[","stored_result")

    #Get the first cancer and additional cancers for all individuals (male)
    first_site_male=lapply(results_list_male,"[[","first_result")
    additional_sites_male=lapply(results_list_male,"[[","stored_result")


    # Combine all individual results (first cancers)
    combined_first_results_males <- do.call(rbind, first_site_male)%>%mutate(sex="Male")
    combined_first_results_females <- do.call(rbind, first_site_female)%>%mutate(sex="Female")
    combined_first_results=bind_rows(combined_first_results_males,combined_first_results_females)%>%
      mutate(start_age=starting_age,end_time=ending_age)


    # ==============================================================================
    # NEW: ported from the original MCEDsim package's
    # sim_multiple_individuals_MCED_parallel_universe(). Converts the raw,
    # uncensored simulated times above into actual usable outcomes, accounting
    # for the competing risk of other-cause death and study end. This is the
    # step that produces the real, reportable analysis variables --
    # life_years_diff and overdiagnosis in particular -- that everything we
    # computed by hand earlier in this project (e.g. cancer_death_time_no_screen
    # - clinical_diagnosis_time) was only a rough stand-in for.
    #
    # Column requirements confirmed compatible with our current
    # combined_first_results (other_cause_death_time, clinical_diagnosis_time,
    # end_time, clinical_diagnosis_stage, screen_diagnosis_time,
    # screen_diagnosis_stage, cancer_death_time_no_screen,
    # cancer_death_time_screen -- all already present). No changes needed to
    # the logic itself; ported as-is from the original.
    # ==============================================================================
    combined_first_results = combined_first_results %>% mutate(
      clin_dx_age = pmin(other_cause_death_time, clinical_diagnosis_time, end_time, na.rm = T),
      clin_dx_event = case_when(
        clin_dx_age == other_cause_death_time ~ "other_cause_death",
        clin_dx_age == end_time ~ "censor",
        clin_dx_age == clinical_diagnosis_time ~ "clin_cancer_diagnosis",
        .default = NA
      ),
      clin_dx_event_stage = case_when(
        clin_dx_event == "clin_cancer_diagnosis" & clinical_diagnosis_stage == "Early" ~ 1,
        clin_dx_event == "clin_cancer_diagnosis" & clinical_diagnosis_stage == "Late" ~ 2,
        .default = 3
      ),
      screen_dx_age = pmin(other_cause_death_time, screen_diagnosis_time, end_time, na.rm = T),
      screen_dx_event = case_when(
        screen_dx_age == other_cause_death_time ~ "other_cause_death",
        screen_dx_age == end_time ~ "censor",
        screen_dx_age == screen_diagnosis_time ~ "screen_cancer_diagnosis",
        .default = NA
      ),
      screen_dx_event_stage = case_when(
        screen_dx_event == "screen_cancer_diagnosis" & screen_diagnosis_stage == "Early" ~ 1,
        screen_dx_event == "screen_cancer_diagnosis" & screen_diagnosis_stage == "Late" ~ 2,
        .default = 3
      ),
      death_age_no_screen = pmin(other_cause_death_time, cancer_death_time_no_screen, end_time, na.rm = T),
      death_age_screen = pmin(other_cause_death_time, cancer_death_time_screen, end_time, na.rm = T),
      death_event_no_screen = case_when(
        death_age_no_screen == other_cause_death_time ~ "other_cause_death",
        death_age_no_screen == end_time ~ "censor",
        death_age_no_screen == cancer_death_time_no_screen ~ "cancer_death",
        .default = NA
      ),
      death_event_screen = case_when(
        death_age_screen == other_cause_death_time ~ "other_cause_death",
        death_age_screen == end_time ~ "censor",
        death_age_screen == cancer_death_time_screen ~ "cancer_death",
        .default = NA
      ),
      diagnosis_age_screen_scenario = pmin(clin_dx_age, screen_dx_age, na.rm = T),
      diagnosis_event_screen_scenario = ifelse(screen_dx_age <= clin_dx_age, screen_dx_event, clin_dx_event),
      diagnosis_event_stage_screen_scenario = case_when(
        screen_dx_event == "screen_cancer_diagnosis" & screen_diagnosis_stage == "Early" ~ 1,
        screen_dx_event == "screen_cancer_diagnosis" & screen_diagnosis_stage == "Late" ~ 2,
        (screen_dx_event != "screen_cancer_diagnosis" & clin_dx_event == "clin_cancer_diagnosis") & clinical_diagnosis_stage == "Early" ~ 1,
        (screen_dx_event != "screen_cancer_diagnosis" & clin_dx_event == "clin_cancer_diagnosis") & clinical_diagnosis_stage == "Late" ~ 2,
        .default = 3
      ),
      life_years_diff = death_age_screen - death_age_no_screen,
      overdiagnosis = ifelse(screen_dx_event == "screen_cancer_diagnosis" & clin_dx_event == "other_cause_death", 1, 0)
    )
    # ==============================================================================
    # END ported block.
    # ==============================================================================


    # Combine all individual results (additional cancers)
    combined_additional_results_males <- do.call(rbind, additional_sites_male)%>%mutate(sex="Male")
    combined_additional_results_females <- do.call(rbind, additional_sites_female)%>%mutate(sex="Female")
    combined_additional_results=bind_rows(combined_additional_results_males,combined_additional_results_females)%>%
      mutate(start_age=starting_age,end_time=ending_age)

    return(list(
      combined_additional_results=combined_additional_results,
      combined_first_results=combined_first_results
    ))
  }












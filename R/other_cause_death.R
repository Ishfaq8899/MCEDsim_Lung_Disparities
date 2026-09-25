# Load necessary libraries
#library(dplyr)
#library(tidyr)


#' Filter MCED cancer data for selected cancer types
#'
#' This function filters the input MCED cancer data to include only the cancer types
#' specified in the `cancer_sites` vector. It uses the `filter` function from `dplyr`
#' to select rows where the cancer type matches any of the types listed in `cancer_sites`.
#'
#' @param data A data frame containing MCED cancer mortality data. This data frame should
#'             include a column named `ICD-10 113 Cause List` which specifies the cancer types.
#' @param cancer_sites A vector of selected cancer type names to filter. These names
#'                     should match the values in the `ICD-10 113 Cause List` column of the data frame.
#'
#' @return A data frame containing only rows where cancer type is in `cancer_sites`.
#' @export
#' @examples
#' selected_cancers <- c(
#' "Malignant neoplasm of prostate (C61)",
#' "Malignant neoplasm of bladder (C67)"
#' )
#' filtered_data <- filter_multiple_cancers(MCED_cdc, selected_cancers)
# filter_multiple_cancers <- function(data, cancer_sites) {
#   filtered_data <-  data %>% mutate(cancer_site=case_when(
#     `ICD-10 113 Cause List`==   "Malignant neoplasms of anus(C21)" ~ "Anus",
#     `ICD-10 113 Cause List`== "Malignant neoplasm of bladder (C67)" ~ "Bladder",
#     `ICD-10 113 Cause List`==  "Malignant neoplasm of breast (C50)" ~ "Breast",
#     `ICD-10 113 Cause List`==   "Malignant neoplasms of colon, and rectum(C18, C19)" ~ "Colorectal",
#     `ICD-10 113 Cause List`==    "Malignant neoplasm of esophagus (C15)" ~ "Esophagus",
#     `ICD-10 113 Cause List`==   "Malignant neoplasm of stomach (C16)" ~ "Gastric",
#     `ICD-10 113 Cause List`=="Malignant neoplasms of lip, oral cavity and pharynx (C00-C14)" ~"Headandneck",
#     `ICD-10 113 Cause List`==  "Malignant neoplasms of liver and intrahepatic bile ducts (C22)" ~ "Liver",
#     `ICD-10 113 Cause List`==   "Malignant neoplasms of trachea, bronchus and lung (C33-C34)" ~ "Lung",
#     `ICD-10 113 Cause List`== "Malignant neoplasm of ovary (C56)" ~ "Ovary",
#     `ICD-10 113 Cause List`==   "Malignant neoplasm of pancreas (C25)" ~ "Pancreas",
#     `ICD-10 113 Cause List`== "Malignant neoplasm of prostate (C61)" ~ "Prostate",
#     `ICD-10 113 Cause List`== "Malignant neoplasms of kidney and renal pelvis (C64-C65)" ~ "Renal",
#     `ICD-10 113 Cause List`==  "Malignant neoplasms of corpus uteri and uterus, part unspecified (C54-C55)" ~ "Uterine",
#     .default=NA)) %>%
#     filter(cancer_site %in% cancer_sites)
#   return(filtered_data)
# }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NEW: CDC/HMD mortality data has no subtype distinction
# everything is just "Lung" at the ICD-10 level. So before filtering, collapse any
# Lung_* subtype names (Lung_Adeno, Lung_LargeCell, Lung_SmallCell,  Lung_Squamous)
#  back down to "Lung" for matching purposes only. This does NOT change what gets returned
# in the data -- cancer_site is still  assigned as "Lung" by the case_when below, exactly as before.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filter_multiple_cancers <- function(data, cancer_sites) {
  cancer_sites_collapsed <- ifelse(startsWith(cancer_sites, "Lung_"), "Lung", cancer_sites)
  filtered_data <-  data %>% mutate(cancer_site=case_when(
    `ICD-10 113 Cause List`==   "Malignant neoplasms of anus(C21)" ~ "Anus",
    `ICD-10 113 Cause List`== "Malignant neoplasm of bladder (C67)" ~ "Bladder",
    `ICD-10 113 Cause List`==  "Malignant neoplasm of breast (C50)" ~ "Breast",
    `ICD-10 113 Cause List`==   "Malignant neoplasms of colon, and rectum(C18, C19)" ~ "Colorectal",
    `ICD-10 113 Cause List`==    "Malignant neoplasm of esophagus (C15)" ~ "Esophagus",
    `ICD-10 113 Cause List`==   "Malignant neoplasm of stomach (C16)" ~ "Gastric",
    `ICD-10 113 Cause List`=="Malignant neoplasms of lip, oral cavity and pharynx (C00-C14)" ~"Headandneck",
    `ICD-10 113 Cause List`==  "Malignant neoplasms of liver and intrahepatic bile ducts (C22)" ~ "Liver",
    `ICD-10 113 Cause List`==   "Malignant neoplasms of trachea, bronchus and lung (C33-C34)" ~ "Lung",
    `ICD-10 113 Cause List`== "Malignant neoplasm of ovary (C56)" ~ "Ovary",
    `ICD-10 113 Cause List`==   "Malignant neoplasm of pancreas (C25)" ~ "Pancreas",
    `ICD-10 113 Cause List`== "Malignant neoplasm of prostate (C61)" ~ "Prostate",
    `ICD-10 113 Cause List`== "Malignant neoplasms of kidney and renal pelvis (C64-C65)" ~ "Renal",
    `ICD-10 113 Cause List`==  "Malignant neoplasms of corpus uteri and uterus, part unspecified (C54-C55)" ~ "Uterine",
    .default=NA)) %>%
    filter(cancer_site %in% cancer_sites_collapsed)
  return(filtered_data)
}

#' Adjust all-cause mortality to estimate other-cause mortality rates
#'
#' This function calculates mortality rates for causes other than the selected (MCED) cancers.
#' All cause and cause-specific (MCED) death rates for 2018-2022 were obtained from CDC (Centers
#' for Disease Control and Prevention 2025) and HMD (Human Mortality Database 2025) by age and gender.
#'
#' The function performs several steps:
#'
#' 1. Filters the MCED cancer data for the selected cancers.
#' 2. Converts HMD data to a long format for easier processing.
#' 3. Process all cause CDC and HMD data to generate all cause death rates
#' 4. Convert all-cause death rates to probabilities using formula of Rosenberg (2006).
#' 5. Process MCED CDC and HMD cancer data to generate MCED cancer death rates
#' 6. Convert MCED cancer death rates to probabilities using formula of Rosenberg (2006).
#' 7. Compute probability of death from causes other than MCED cancers
#' 8. Convert other cause death probabilities to rates using inverse formula of Rosenberg (2006).
#'
#'
#' @param cdc_data A data frame with all-cause mortality data from the CDC. This data frame should
#'                      include columns: "Single-Year Ages Code", "Gender", "Year", "Deaths", "Population", "Crude Rate".
#' @param MCED_cdc A data frame with cancer-specific mortality data from the CDC. This data frame should
#'                 include columns: "Single-Year Ages Code", "Sex", "Year", "ICD-10 113 Cause List",
#'                 "ICD-10 113 Cause List Code", "Deaths", "Population", "Crude Rate".
#' @param hmd_data A data frame with population data from HMD. This data frame should include columns:
#'                 "Year", "Age", "Female", "Male", "Total".
#' @param selected_cancers A character vector specifying the cancer types to analyze.
#'  These names include: "Bladder","Breast","Colorectal", "Esophagus","Gastric", "Headandneck",
#'  "Liver", "Lung","Ovary", "Pancreas", "Prostate", "Renal", "Uterine"
#'
#' @return A data frame with age-specific and sex-specific mortality rates due to causes
#'         other than the selected cancers. Includes columns: age, sex, year, other_cause_rate.
#' @export
#' @import tidyr
#' @examples
#' cdc_data <- read.csv("all_cause_cdc.csv")
#' MCED_cdc <- read.csv("MCED_cdc.csv")
#' hmd_data <- read.csv("hmd_data.csv")
#' selected_cancers <- c(
#' "Prostate",
#' "Bladder"
#' )
#' get_other_cause_mortality <- adjusted_all_cause_mortality(cdc_data, MCED_cdc, hmd_data, selected_cancers)
get_other_cause_mortality <- function(cdc_data, MCED_cdc, hmd_data, selected_cancers) {

  # Filter the MCED data for the selected cancers
  filtered_cancer_data <- filter_multiple_cancers(data = MCED_cdc, cancer_sites = selected_cancers)

  # Convert HMD data to long format for easier processing
  hmd_data_long <- hmd_data %>%
    pivot_longer(cols = c(Female, Male), names_to = "Gender", values_to = "Population") %>%
    filter(Year %in% 2018:2022) %>%
    select(Year, Age, Gender, Population) %>%
    mutate(Age = as.numeric(Age))

  # Process all cause CDC and HMD data to generate all cause death rates
  all_cause_death_rate <- cdc_data %>%
    rename(Age = `Single-Year Ages Code`) %>%
    rename(all_crude_Rate = `Crude Rate`) %>%
    select(Notes, Age, Gender, Year, Deaths, Population, all_crude_Rate)
  # %>% filter(!is.na(all_crude_Rate), !is.na(Gender), !is.na(Year))


  all_cause_death_rate<- all_cause_death_rate %>%
    left_join(hmd_data_long %>% filter(Age > 84), by = c("Year", "Age", "Gender"), suffix = c("", "_hmd")) %>%
    mutate(Population = ifelse(is.na(Population), Population_hmd, Population)) %>%
    select(-Population_hmd) %>%
    mutate(all_crude_Rate = ifelse(Age > 84, (Deaths / Population) * 100000, all_crude_Rate)) %>%
    rename(sex = Gender, age = Age, year = Year) %>%
    mutate(Type = "All-cause death rates") %>%
    filter(!is.na(all_crude_Rate), !is.na(sex), !is.na(year))

  # Convert all-cause death rates to probabilities
  all_cause_death_probability <- all_cause_death_rate %>%
    mutate(all_death_probability = (all_crude_Rate / 100000) / (1 + 0.5 * all_crude_Rate / 100000)) %>%
    mutate(Type = "All-cause death probability")

  # Process MCED CDC and HMD cancer data to generate MCED cancer death rates
  MCED_cancer_death_rate <- filtered_cancer_data %>%
    rename(Age = `Single-Year Ages Code`, Gender = Sex) %>%
    select(Year, Age, Gender, cancer_site, Deaths, Population) %>%
    left_join(hmd_data_long %>% filter(Age > 84), by = c("Year", "Age", "Gender"), suffix = c("", "_hmd")) %>%
    mutate(Population = ifelse(is.na(Population), Population_hmd, Population)) %>%
    select(-Population_hmd) %>%
    mutate(crude_Rate =  (Deaths / Population) * 100000) %>%
    mutate(crude_Rate = as.numeric(crude_Rate)) %>%
    rename(sex = Gender, age = Age, year = Year) %>%
    group_by(year, age, sex) %>%
    summarize(
      MCED_Deaths = sum(Deaths, na.rm = TRUE),
      Population = first(Population),
      MCED_crude_Rate = sum(crude_Rate, na.rm = TRUE),
      .groups = 'drop'
    )

  # Convert MCED cancer death rates to probabilities
  MCED_cause_death_probability <- all_cause_death_rate %>%
    left_join(MCED_cancer_death_rate %>% select(year, age, sex, MCED_crude_Rate), by = c("year", "age", "sex")) %>%
    mutate(MCED_death_probability = (MCED_crude_Rate / 100000) / (1 + 0.5 * all_crude_Rate / 100000)) %>%
    select(age, sex, year, all_crude_Rate, MCED_crude_Rate, MCED_death_probability)


  # Compute probability of death from causes other than MCED cancers
  other_cause_death_probability <- all_cause_death_probability %>%
    left_join(MCED_cause_death_probability %>% select(year, age, sex, MCED_death_probability), by = c("year", "age", "sex")) %>%
    mutate(other_cause_probability = if_else(is.na(MCED_death_probability),
                                             all_death_probability,
                                             all_death_probability - MCED_death_probability)) %>%
    select(age, sex, year, all_death_probability, MCED_death_probability, other_cause_probability) %>%
    mutate(Type = "Other-cause death probability")


  # Convert other cause death probabilities to rates
  other_cause_death_rate <- other_cause_death_probability %>%
    mutate(other_cause_rate = (other_cause_probability / (1 - 0.5 * other_cause_probability)) * 100000) %>%
    select(age, sex, year, other_cause_rate) %>%
    mutate(Type = "Other-cause death rate")

  return(other_cause_death_rate)
}

#################################################
#' Create a table of other-cause death survival probabilities
#'
#' This function generates a table of survival probabilities for causes other than selected (MCED) cancers.
#' Its uses the output from the `get_other_cause_mortality` function and filters it based on the specified
#' sex, year, and starting age.
#'
#' @param cdc_data A data frame with all-cause mortality data from the CDC. This data frame should
#'                      include columns: "Single-Year Ages Code", "Gender", "Year", "Deaths", "Population",
#'                      "Crude Rate".
#' @param MCED_cdc A data frame with cancer-specific mortality data from the CDC. This data frame should
#'                    include columns: "Single-Year Ages Code", "Sex", "Year", "ICD-10 113 Cause List",
#'                    "ICD-10 113 Cause List Code", "Deaths", "Population", "Crude Rate".
#' @param hmd_data A data frame with population data from HMD. This data frame should include columns:
#'                  "Year", "Age", "Female", "Male", "Total".
#' @param selected_cancers A character vector specifying the cancer types to analyze.
#'                         These names include: "Anus", "Bladder", "Breast", "Colorectal", "Esophagus", "Gastric",
#'                         "Headandneck", "Liver", "Lung", "Ovary", "Pancreas", "Prostate", "Renal", "Uterine".
#'
#'
#' @param the_sex  A character string specifying the sex to filter by ("Male" or "Female").
#' @param the_year An integer specoifying the year to filter by.
#' @param the_starting_age An integer specifying the starting age for the survival probabilities.
#'
#' @return A data frame with survival probabilities for causes other than the selected cancers.
#'         Includes columns: age, sex, year, other_cause_rate, surv.
#' @export
#' @examples
#' cdc_data <- read_excel("/path/to/cdc_all_cause_2018_2022.xlsx")
#' MCED_cdc <- read_excel("/path/to/modified_MCED_data.xlsx")
#' hmd_data <- read_excel("/path/to/hmd_population_1933_2023.xlsx")
#' selected_cancers <- c("Anus", "Bladder", "Breast", "Colorectal", "Esophagus", "Gastric", "Headandneck",
#'                       "Liver", "Lung", "Ovary", "Pancreas", "Prostate", "Renal", "Uterine")
#' othercause_surv <- make_othercause_death_table(cdc_data, MCED_cdc, hmd_data, selected_cancers, "Female", 2019, 60)
make_othercause_death_table <- function(cdc_data, MCED_cdc, hmd_data, selected_cancers,
                                        the_sex, the_year, the_starting_age) {

  othercause_surv=get_other_cause_mortality(cdc_data = cdc_data, MCED_cdc = MCED_cdc,
                                            hmd_data = hmd_data, selected_cancers = selected_cancers) %>%
    filter(sex == the_sex, year == the_year)%>%
    mutate(surv=cumprod(1-other_cause_rate/100000))%>%
    mutate(surv=ifelse(age<=the_starting_age, 1, surv/surv[age==the_starting_age]))

  return(othercause_surv)

}

#' Simulate other-cause death time
#'
#' This function simulates the time to death from causes other than the selected (MCED) cancers
#' using the survival probabilities from the `make_othercause_death_table` function.
#'
#' @param othercause_death_table A data frame with survival probabilities for causes other than the selected cancers.
#'                               Includes columns: age, sex, year, other_cause_rate, surv.
#' @param ID optional ID to set random seed.
#'
#' @return A numeric value representing the simulated time to death from other causes.
#' @export
#' @examples
#' othercause_death_table <- make_othercause_death_table(cdc_data, MCED_cdc, hmd_data, selected_cancers, "Female", 2019, 60)
#' sim_time <- sim_othercause_death(othercause_death_table)
sim_othercause_death <- function(othercause_death_table,ID=NA) {
  if(!is.na(ID)){
    set.seed(ID)
  }
  the_time=gettime(time = othercause_death_table$age, surv = othercause_death_table$surv)

  return(the_time)

}


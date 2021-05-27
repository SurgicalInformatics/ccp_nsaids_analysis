#set up analysis
set.seed(100) #reproduce bootstrapping / MI
# explan_severity = c('news', 'news_cat', 'qsofa', 'qsofa_cat', 'hr_vsorres', 'rr_vsorres', 'sysbp_vsorres', 'admission_diabp_vsorres', 'oxy_vsorres')

source('02_functions.R') # functions for tables etc.

#Load in libraries
library(tidyverse)
library(mice)
library(MatchThem)
library(knitr)
library(kableExtra)
library(finalfit)
library(survival)
library(ggplot2)
library(lubridate)
library(survminer)
library(ggsci)
library(patchwork)
library(janitor)

#First remove any sort of cream / gels from the analysis
nsaid_cream_gel = medication_out %>% 
  ungroup() %>% 
  filter(match > 0.84) %>% 
  select(subjid, entered_drug_name, medication_label) %>% 
  mutate(is_gel = ifelse(grepl('gel|cream|%', entered_drug_name, ignore.case = T) & (medication_label == 'ibuprofen' |
                                                                                       medication_label == 'diclofenac') , 'Is gel', 'Not gel/cream')) %>% 
  select(subjid, is_gel) %>% 
  group_by(subjid) %>% 
  mutate(is_gel_cream = ifelse(any(is_gel == 'Is gel'), 'gel/cream', NA)) %>% select(subjid, is_gel_cream) %>%  distinct(subjid, .keep_all = T)


# Select those with at least 14 days follow-up and remove implausible / age errors
topline = topline %>% 
  #select(-hosttim) %>% 
  filter(subjid %in% keep_14) %>%
  filter(!is.na(death)) %>%
  ff_relabel_df(topline) %>% 
  mutate_if(is.factor, fct_recode, NULL = "Unknown") %>% 
  mutate_if(is.factor, fct_recode, NULL = "Not specified") %>%
  mutate_if(is.factor, fct_recode, NULL = "N/K") %>%
  mutate_if(is.factor, fct_recode, NULL = "N/A") %>% 
  mutate_if(is.factor, fct_explicit_na) %>% 
  mutate_if(is.factor, fct_relevel, 'NO') %>% 
  mutate(age_estimateyears = ifelse(age_estimateyears < 120, age_estimateyears, NA))

# Do the same for the survival data so this has the same patients as the main topline dataset
surv_data = surv_data %>% 
  select(-hosttim) %>% 
  filter(subjid %in% keep_14) %>%
  filter(!is.na(status)) %>%
  ff_relabel_df(surv_data) %>%
  mutate_if(is.factor, fct_recode, NULL = "Unknown") %>% 
  mutate_if(is.factor, fct_recode, NULL = "Not specified") %>%
  mutate_if(is.factor, fct_recode, NULL = "N/K") %>%
  mutate_if(is.factor, fct_recode, NULL = "N/A") %>% 
  mutate_if(is.factor, fct_explicit_na) %>% 
  mutate_if(is.factor, fct_relevel, 'NO') %>% 
  mutate(age_estimateyears = ifelse(age_estimateyears < 120, age_estimateyears, NA))

# Do the same for the main, very large data so this has the same patients as the main topline dataset
ccp_data = ccp_data %>% 
  select(-hosttim) %>% 
  mutate(age_estimateyears = ifelse(age_estimateyears < 120, age_estimateyears, NA))


#Put reasonable limits on data to correct any input errors

topline = topline %>% 
  mutate(
    age = case_when(
      age > 110 ~ NA_real_,
      TRUE ~ age
    ), 
    temp_vsorres = case_when(
      temp_vsorres > 3000 ~ temp_vsorres / 100,
      temp_vsorres > 300 ~ temp_vsorres / 10,
      temp_vsorres  > 200 ~ temp_vsorres - 200,
      temp_vsorres > 100 ~ temp_vsorres - 100,
      temp_vsorres > 45 ~ NA_real_,
      temp_vsorres < 0 ~ abs(temp_vsorres),
      temp_vsorres < 3 ~ NA_real_,
      temp_vsorres < 4 ~ temp_vsorres * 10,
      temp_vsorres < 25 ~ NA_real_,
      TRUE ~ temp_vsorres
    ), 
    hr_vsorres = case_when(
      hr_vsorres > 300 ~ NA_real_,
      hr_vsorres < 30 ~ NA_real_,
      TRUE ~ hr_vsorres
    ),
    rr_vsorres = case_when(
      rr_vsorres > 150 ~ NA_real_,
      TRUE ~ rr_vsorres
    ),
    daily_gcs_vsorres = case_when(
      daily_gcs_vsorres < 3 ~ 3,
      daily_gcs_vsorres > 15 ~ 15,
      TRUE ~ daily_gcs_vsorres
    ),
    oxy_vsorres  = case_when(
      oxy_vsorres < 0 ~ abs(oxy_vsorres), 
      oxy_vsorres < 50 ~ NA_real_,
      oxy_vsorres > 100 ~ NA_real_,
      TRUE ~ oxy_vsorres),
    daily_plt_lborres = ifelse(daily_plt_lborres > 2000, NA_real_,  daily_plt_lborres),
    admission_diabp_vsorres = case_when(
      admission_diabp_vsorres < 20 ~ NA_real_,
      admission_diabp_vsorres > 200 ~ NA_real_,
      TRUE ~ admission_diabp_vsorres),
    sysbp_vsorres = case_when(
      sysbp_vsorres < 10 ~ NA_real_,
      sysbp_vsorres > 300 ~ NA_real_,
      TRUE ~ sysbp_vsorres),
    daily_fio2_lborres = case_when(
      daily_fio2_lborres > 1 ~ daily_fio2_lborres / 100,
      TRUE ~ daily_fio2_lborres),
    daily_pao2_lborres = case_when(
      daily_pao2_lborres > 100 ~ NA_real_,
      TRUE ~ daily_pao2_lborres),
    daily_ph_lborres  = case_when(
      daily_ph_lborres < 6 ~ NA_real_,
      TRUE ~ daily_ph_lborres),
    daily_urine_lborres = case_when(
      daily_urine_lborres < 0 ~ abs(daily_urine_lborres),
      TRUE ~ daily_urine_lborres),
    daily_bun_lborres = case_when(
      daily_bun_lborres > 100 ~ 100,
      TRUE ~ daily_bun_lborres),
    daily_hb_lborres = case_when(
      daily_hb_lborres > 300 ~ NA_real_,
      daily_hb_lborres < 20 ~ NA_real_,
      TRUE ~ daily_hb_lborres),
    daily_neutro_lborres = case_when(
      daily_neutro_lborres >= 100 ~ NA_real_,
      TRUE ~ daily_neutro_lborres),
    daily_lymp_lborres = case_when(
      daily_lymp_lborres >= 100 ~ NA_real_,
      TRUE ~ daily_lymp_lborres),
    daily_pt_lborres_add_inr = case_when(
      daily_pt_lborres_add_inr > 150 ~ NA_real_,
      daily_pt_lborres_add_inr < 9 ~ NA_real_,
      TRUE ~ daily_pt_lborres_add_inr),
    daily_aptt_lborres = case_when(
      daily_aptt_lborres > 150 ~ NA_real_,
      daily_aptt_lborres < 4 ~ daily_aptt_lborres * 22,
      TRUE ~ daily_aptt_lborres),
    daily_sodium_lborres = case_when(
      daily_sodium_lborres > 180 ~ NA_real_,
      daily_sodium_lborres < -100 ~ abs(daily_sodium_lborres),
      daily_sodium_lborres < 100 ~ NA_real_,
      TRUE ~ daily_sodium_lborres),
    daily_potassium_lborres = case_when(
      daily_potassium_lborres < 0.55 ~ daily_potassium_lborres * 10,
      TRUE ~ daily_potassium_lborres),
    daily_crp_lborres = case_when(
      daily_crp_lborres < 0 ~ abs(daily_crp_lborres),
      daily_crp_lborres > 750 ~ 750,
      TRUE ~ daily_crp_lborres),
    daily_lactate_lborres= case_when(
      daily_lactate_lborres < 0 ~ abs(daily_lactate_lborres),
      TRUE ~ daily_lactate_lborres)
  ) %>% 
  ff_relabel_df(topline)

# Calculate NEWS and qsofa

topline = topline %>% 
  mutate_if(is.factor, fct_recode, NULL = "Unknown") %>% 
  mutate_if(is.factor, fct_relevel, 'NO') %>% 
  mutate(
    # For NEWS
    alt_conscious = ifelse(daily_gcs_vsorres < 15, "Yes", "No"),
    hypoxic_target = "no", # This and next line for COPD patients only so not used.
    o2_rx = "no",
  ) %>% 
  news(rr = rr_vsorres, spo2 = oxy_vsorres, o2_rx = o2_rx,
       hypoxic_target = hypoxic_target, sbp = sysbp_vsorres, hr = hr_vsorres, 
       temp = temp_vsorres, alt_conscious = alt_conscious, 
       output = "df_vector",  na_to_zeros = FALSE) %>%
  
  qsofa(rr = rr_vsorres, sbp = sysbp_vsorres, gcs = daily_gcs_vsorres,
        na_to_zeros = FALSE, output = c("df_vector"))

surv_data = surv_data %>% 
  mutate(
    # For NEWS
    alt_conscious = ifelse(daily_gcs_vsorres < 15, "Yes", "No"),
    hypoxic_target = "no", # This and next line for COPD patients only so not used.
    o2_rx = "no",
  ) %>% 
  news(rr = rr_vsorres, spo2 = oxy_vsorres, o2_rx = o2_rx,
       hypoxic_target = hypoxic_target, sbp = sysbp_vsorres, hr = hr_vsorres, 
       temp = temp_vsorres, alt_conscious = alt_conscious, 
       output = "df_vector",  na_to_zeros = FALSE) %>%
  
  qsofa(rr = rr_vsorres, sbp = sysbp_vsorres, gcs = daily_gcs_vsorres,
        na_to_zeros = FALSE, output = c("df_vector")) %>% 
  ff_relabel_df(ccp_data_labels) 

# Now generate the NSAIDs variables
# Briefly, this looks to see if the patient is on topical NSAIDs and will then require another NSAID (likely systemic) in order to be classed as an NSAID user
# Aspirin, for this, is not an NSAID as its mainly used in UK for antiplatelet indications

nsaid_cream_gel = nsaid_cream_gel %>%  ungroup() %>%  
  filter(subjid %in% topline$subjid) 

#define variables
nsaids_data = topline %>% 
  select(subjid, propionic_acid_derivatives_antiinflammatory_and_antirheumatic_products, 
         ibuprofen, diclofenac, diclofenac_combinations, ketorolac, antiinflammatory_agents_non_steroidal_ophthalmologic,
         naproxen, celecoxib, etoricoxib, oxicams_antiinflammatory_and_antirheumatic_drugs, ibuprofen_combinations, parecoxib, chronic_nsaid_cmoccur) %>% 
  #mutate_all(as.character) %>% 
  replace(is.na(.), 'No') %>% 
  left_join(nsaid_cream_gel, by = 'subjid') %>% 
  mutate(  nsaids = ifelse(is_gel_cream == 'gel/cream', 'No NSAIDs', NA),
           nsaids = ifelse((propionic_acid_derivatives_antiinflammatory_and_antirheumatic_products == 'Yes' |  ibuprofen_combinations == 'Yes' |
                              ibuprofen == 'Yes' |
                              diclofenac == 'Yes' | 
                              diclofenac_combinations == 'Yes' | 
                              ketorolac == 'Yes' |
                              antiinflammatory_agents_non_steroidal_ophthalmologic == 'Yes' |
                              naproxen == 'Yes' |
                              celecoxib == 'Yes' |
                              parecoxib == 'Yes' |
                              etoricoxib == 'Yes' |
                              oxicams_antiinflammatory_and_antirheumatic_drugs == 'Yes') & is.na(is_gel_cream), 'Mixed NSAIDs', nsaids),
           nsaids = ifelse(((propionic_acid_derivatives_antiinflammatory_and_antirheumatic_products == 'Yes' &  ibuprofen_combinations == 'No' )|
                              ibuprofen == 'Yes') &
                             (diclofenac == 'No' & 
                                diclofenac_combinations == 'No' & 
                                ketorolac == 'No' &
                                antiinflammatory_agents_non_steroidal_ophthalmologic == 'No' &
                                naproxen == 'No' &
                                celecoxib == 'No' &
                                parecoxib == 'No' &
                                etoricoxib == 'No' &
                                oxicams_antiinflammatory_and_antirheumatic_drugs == 'No') & is.na(is_gel_cream), 'Ibuprofen (or similar propionic acid derivative)', nsaids),
           nsaids = ifelse((propionic_acid_derivatives_antiinflammatory_and_antirheumatic_products == 'No' &
                              ibuprofen == 'No' &
                              (diclofenac == 'Yes' |
                                 diclofenac_combinations == 'Yes' |
                                 ketorolac == 'Yes' |
                                 antiinflammatory_agents_non_steroidal_ophthalmologic == 'Yes' |
                                 naproxen == 'Yes') &
                              celecoxib == 'No' &
                              parecoxib == 'No' &
                              etoricoxib == 'No' &
                              oxicams_antiinflammatory_and_antirheumatic_drugs == 'No') & is.na(is_gel_cream), 'Diclofenac/ Ketorolac /Naproxen', nsaids),
           nsaids = ifelse(propionic_acid_derivatives_antiinflammatory_and_antirheumatic_products == 'No' &
                             ibuprofen == 'No' &
                             diclofenac == 'No' & 
                             diclofenac_combinations == 'No' & 
                             ketorolac == 'No' &
                             antiinflammatory_agents_non_steroidal_ophthalmologic == 'No' &
                             naproxen == 'No' &
                             (celecoxib == 'Yes' |
                                parecoxib == 'Yes' |
                                etoricoxib == 'Yes') &
                             oxicams_antiinflammatory_and_antirheumatic_drugs == 'No', 'COX-2 Inhibitor', nsaids),
           nsaids = ifelse(propionic_acid_derivatives_antiinflammatory_and_antirheumatic_products == 'No' &
                             ibuprofen == 'No' &
                             diclofenac == 'No' & 
                             diclofenac_combinations == 'No' & 
                             ketorolac == 'No' &
                             antiinflammatory_agents_non_steroidal_ophthalmologic == 'No' &
                             naproxen == 'No' &
                             celecoxib == 'No' &
                             parecoxib == 'No' &
                             etoricoxib == 'No' &
                             oxicams_antiinflammatory_and_antirheumatic_drugs == 'Yes', 'Oxicam', nsaids),
           nsaids = ifelse(propionic_acid_derivatives_antiinflammatory_and_antirheumatic_products == 'No' &
                             ibuprofen == 'No' &
                             diclofenac == 'No' & 
                             diclofenac_combinations == 'No' & 
                             ketorolac == 'No' &
                             antiinflammatory_agents_non_steroidal_ophthalmologic == 'No' &
                             naproxen == 'No' &
                             celecoxib == 'No' &
                             parecoxib == 'No' &
                             etoricoxib == 'No' &
                             oxicams_antiinflammatory_and_antirheumatic_drugs == 'No', 'No NSAIDs', nsaids) %>% fct_relevel('No NSAIDs', 'Ibuprofen (or similar propionic acid derivative)', 'Diclofenac/ Ketorolac /Naproxen', 
                                                                                                                            'COX-2 Inhibitor', 'Oxicam')) %>% 
  mutate(nsaids_yn = case_when(nsaids == 'Ibuprofen (or similar propionic acid derivative)' |
                                 nsaids == 'Diclofenac/ Ketorolac /Naproxen' |
                                 nsaids == 'COX-2 Inhibitor' |
                                 nsaids == 'Oxicam'|
                                 nsaids == 'Mixed NSAIDs' |
                                 chronic_nsaid_cmoccur == 'Yes' ~ 'NSAIDs', 
                               nsaids == 'No NSAIDs' |  chronic_nsaid_cmoccur == 'No'|   chronic_nsaid_cmoccur == '(Missing)'~ 'No NSAIDs') %>% factor() %>% ff_label('NSAIDs')) %>% 
  select(subjid, nsaids, nsaids_yn)

#Now tidy up some of the other variables, relevel the factors etc.
topline = topline %>% 
  ungroup() %>% 
  mutate_if(is.factor, fct_recode, 'No' = 'NO') %>% 
  mutate_if(is.factor, fct_recode, 'Yes' = 'YES') %>% 
  mutate_if(is.factor, fct_recode, NULL = 'Unknown') %>% 
  mutate_if(is.factor, fct_recode, NULL = 'N/K') %>% 
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  mutate_if(is.factor, fct_recode, NULL = 'Not specified') %>% 
  left_join(nsaids_data, by = 'subjid') %>% 
  distinct(subjid, .keep_all = T) %>% 
  mutate(diabetes_fct = case_when(
    diabetes_mhyn == 'Yes' ~ 'Diabetes without complications',
    diabetescom_mhyn == 'Yes' ~ 'Diabetes with complications',
    diabetes_mhyn == 'No' & diabetescom_mhyn =='No'~ 'No Diabetes',
    is.na(diabetes_mhyn) & is.na(diabetescom_mhyn) ~ NA_character_,
    is.na(diabetes_mhyn) | is.na(diabetescom_mhyn) ~ 'No Diabetes',
    TRUE  ~ NA_character_)   %>% factor() %>% ff_label('Diabetes')) %>% ff_relabel_df(topline)

surv_data = surv_data %>% 
  ungroup() %>% 
  mutate_if(is.factor, fct_recode, 'No' = 'NO') %>% 
  mutate_if(is.factor, fct_recode, 'Yes' = 'YES') %>% 
  mutate_if(is.factor, fct_recode, NULL = 'Unknown') %>% 
  mutate_if(is.factor, fct_recode, NULL = 'N/K') %>% 
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  mutate_if(is.factor, fct_recode, NULL = 'Not specified') %>% 
  left_join(nsaids_data, by = 'subjid') %>% 
  distinct(subjid, .keep_all = T) %>% 
  mutate(diabetes_fct = case_when(
    diabetes_mhyn == 'Yes' ~ 'Diabetes without complications',
    diabetescom_mhyn == 'Yes' ~ 'Diabetes with complications',
    diabetes_mhyn == 'No' & diabetescom_mhyn =='No'~ 'No Diabetes',
    is.na(diabetes_mhyn) & is.na(diabetescom_mhyn) ~ NA_character_,
    is.na(diabetes_mhyn) | is.na(diabetescom_mhyn) ~ 'No Diabetes',
    TRUE  ~ NA_character_)   %>% factor() %>% ff_label('Diabetes')) %>% ff_relabel_df(topline)


#ISARIC SCORES
# Adds the 4C score 
source('add_isaric_score.R')
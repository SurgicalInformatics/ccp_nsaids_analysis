library(magrittr)

## ----Variable changes, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------
cs = topline %>% 
  mutate_at(vars(chrincard, renal_mhyn, malignantneo_mhyn, modliv, obesity_mhyn, chronicpul_mhyn, malignantneo_mhyn, diabetes_mhyn, diabetescom_mhyn), fct_recode, 
            NULL = "Unknown",
            "No" = "NO",
            "Yes" = "YES") %>%
  mutate_at(vars(chrincard, renal_mhyn, malignantneo_mhyn, modliv, infiltrates_faorres, obesity_mhyn, 
                 diabetes_mhyn, diabetescom_mhyn, chronicpul_mhyn), fct_relevel, "No") %>% 
  ff_relabel_df(topline)

cs = cs %>% 
  mutate(sex = fct_recode(sex, NULL= "Not specified") %>% 
           ff_label("Sex"),
         infiltrates_faorres= fct_recode(infiltrates_faorres, 
                                         NULL = "N/A",
                                         "No" = "NO",
                                         "Yes" = "YES") %>% 
           fct_relevel("No") %>% 
           ff_label("CXR infiltrates"),
         daily_ldh_lborres = ff_label(daily_ldh_lborres, "LDH"),
         mort = ff_label(death, "Mortality"),
         
         # For NEWS
         alt_conscious = ifelse(daily_gcs_vsorres < 15, "Yes", "No"),
         hypoxic_target = "no", # This and next line for COPD patients only so not used.
         o2_rx = "no",
         
         # Neutrophil / Lymphocyte Ration
         NLR = daily_neutro_lborres / daily_lymp_lborres,
  )

### Co-morbidity count

# First sort out diabetes
cs = cs %>% 
  mutate(diabetes_combined = case_when(
    diabetes_mhyn == "Yes" | diabetescom_mhyn == "Yes" ~ "Yes",
    is.na(diabetes_mhyn) & is.na(diabetescom_mhyn) ~ NA_character_,
    TRUE ~ "No"))

# Then count comorbidity (all)
comorbid = cs %>% 
  select(subjid, chrincard, renal_mhyn, modliv, obesity_mhyn, chronicpul_mhyn, 
         malignantneo_mhyn, diabetes_combined)

comorbid$no_comorbid <- rowSums(comorbid == "Yes")


# Then count comorbidity (not including heart failure) for lasso model
comorbid_noheart = cs %>% 
  select(subjid, renal_mhyn, modliv, obesity_mhyn, chronicpul_mhyn, 
         malignantneo_mhyn, diabetes_combined)

comorbid_noheart$no_comorbid_h <- rowSums(comorbid_noheart == "Yes")

removals = tibble()
# Join with cs tibble
cs = cs %>% 
  left_join(comorbid) %T>% 
  {removals <<- bind_rows(removals, tibble(distinct = n_distinct(paste0(.$subjid, .$arm_n)),
                                           nrow = nrow(.),
                                           label = "cs-after left_join(comorbid) "))} %>% 
  left_join(comorbid_noheart) %>% 
  ff_relabel_df(cs) %T>% 
  {removals <<- bind_rows(removals, tibble(distinct = n_distinct(paste0(.$subjid, .$arm_n)),
                                           nrow = nrow(.),
                                           label = "cs-after left_join(comorbid_noheart)"))} 



## ----Create variables for surgisphere score, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------

# Comorbid score
comorbid_surgisphere = cs %>% 
  select(subjid, renal_mhyn, modliv, chronicpul_mhyn, 
         diabetes_combined, hypertension_mhyn, asthma_mhyn, smoking_mhyn_2levels) %>% 
  # YES/NO in smoking
  mutate(smoking = case_when(
    smoking_mhyn_2levels == "YES" ~ "Yes",
    smoking_mhyn_2levels == "NO" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # Hypertension
  mutate(hypertension = case_when(
    hypertension_mhyn == "1" ~ "Yes",
    hypertension_mhyn == "2" ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  # Asthma YES/NO
  mutate(asthma = case_when(
    asthma_mhyn == "YES" ~ "Yes",
    asthma_mhyn == "NO" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # renal
  mutate(renal = case_when(
    renal_mhyn == "Yes" ~ "Yes",
    renal_mhyn == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # Liver
  mutate(liver = case_when(
    modliv == "Yes" ~ "Yes",
    modliv == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # lung
  mutate(lung = case_when(
    chronicpul_mhyn == "Yes" ~ "Yes",
    chronicpul_mhyn == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # diabetes
  mutate(diabetes = case_when(
    diabetes_combined == "Yes" ~ "Yes",
    diabetes_combined == "No" ~ "No",
    TRUE ~ "No"
  ))  %>% 
  select(-hypertension_mhyn, -asthma_mhyn, -smoking_mhyn_2levels, -renal_mhyn,
         -modliv, -chronicpul_mhyn, -diabetes_combined)

# Add up "yes" response
comorbid_surgisphere$no_comorbid_SS <- rowSums(comorbid_surgisphere == "Yes")


# Immunosuppression
immunosuppession_surgisphere = cs %>% 
  select(subjid, malignantneo_mhyn, aidshiv_mhyn, malnutrition_mhyn, vulnerable_transplant,
         immno_cmtrt) %>% 
  # AIDS YES/NO
  mutate(aids = case_when(
    aidshiv_mhyn == "YES" ~ "Yes",
    aidshiv_mhyn == "NO" ~ "No",
    TRUE ~ "No")) %>% 
  # Malnutrition YES/NO
  mutate(malnurition = case_when(
    malnutrition_mhyn == "YES" ~ "Yes",
    malnutrition_mhyn == "NO" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # Cancer
  mutate(cancer = case_when(
    malignantneo_mhyn == "Yes" ~ "Yes",
    malignantneo_mhyn == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # Immuno drugs
  mutate(immuno_d = case_when(
    immno_cmtrt == "Yes" ~ "Yes",
    immno_cmtrt == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # Vulnerable Tx
  mutate(tx = case_when(
    vulnerable_transplant == "Yes" ~ "Yes",
    vulnerable_transplant == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  select(-aidshiv_mhyn, -malnutrition_mhyn, -malignantneo_mhyn, -immno_cmtrt, -vulnerable_transplant)

# Add up & mutate new column
immunosuppession_surgisphere$immuno_SS <- rowSums(immunosuppession_surgisphere == "Yes")

immunosuppession_surgisphere = immunosuppession_surgisphere %>% 
  mutate(immunoSS = case_when(
    immuno_SS >= 1 ~ 1,
    TRUE ~ 0
  )) %>% 
  select(subjid, immuno_SS)


# Surgisphere score
surgisphere_score = cs %>% 
  left_join(comorbid_surgisphere) %>% 
  left_join(immunosuppession_surgisphere) %>% 
  select(-smoking, -hypertension, -asthma, -renal, -liver, -lung, -diabetes) %>% 
  select(subjid, no_comorbid_SS, chrincard, immuno_SS) %>% 
  mutate(chrincard = case_when(
    chrincard == "Yes" ~ "Yes",
    chrincard == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # Score variable
  mutate(comorbid_SS = case_when(
    no_comorbid_SS >=2 | chrincard == "Yes" | immuno_SS == 1 ~ 1,
    TRUE ~ 0
  )) %>% 
  select(subjid, comorbid_SS)

cs = cs %>% 
  left_join(surgisphere_score) %>% 
  ff_relabel_df(cs) %T>% 
  {removals <<- bind_rows(removals, tibble(distinct = n_distinct(paste0(.$subjid, .$arm_n)),
                                           nrow = nrow(.),
                                           label = "cs-after left_join(surgisphere_score)"))}


## ----Add in Covid_gram score, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------

comorbid_gram = cs %>% 
  select(subjid, chrincard, renal_mhyn, modliv, chronicpul_mhyn, 
         diabetes_combined, hypertension_mhyn
  ) %>% 
  # Heart disease
  mutate(heart_disease = case_when(
    chrincard == "Yes" ~ "Yes",
    chrincard == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # Hypertension
  mutate(hypertension = case_when(
    hypertension_mhyn == "1" ~ "Yes",
    hypertension_mhyn == "2" ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  # renal
  mutate(renal = case_when(
    renal_mhyn == "Yes" ~ "Yes",
    renal_mhyn == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # Liver
  mutate(liver = case_when(
    modliv == "Yes" ~ "Yes",
    modliv == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # lung
  mutate(lung = case_when(
    chronicpul_mhyn == "Yes" ~ "Yes",
    chronicpul_mhyn == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # diabetes
  mutate(diabetes = case_when(
    diabetes_combined == "Yes" ~ "Yes",
    diabetes_combined == "No" ~ "No",
    TRUE ~ "No"
  ))  %>% 
  select(-chrincard, -hypertension_mhyn, -renal_mhyn,
         -modliv, -chronicpul_mhyn, -diabetes_combined)


# Add up "yes" response
comorbid_gram$comorbid_gram <- rowSums(comorbid_gram == "Yes")


immunosuppession_covid_gram = cs %>% 
  select(subjid, aidshiv_mhyn, malnutrition_mhyn, vulnerable_transplant,
         immno_cmtrt) %>% 
  # AIDS YES/NO
  mutate(aids = case_when(
    aidshiv_mhyn == "YES" ~ "Yes",
    aidshiv_mhyn == "NO" ~ "No",
    TRUE ~ "No")) %>% 
  # Malnutrition YES/NO
  mutate(malnurition = case_when(
    malnutrition_mhyn == "YES" ~ "Yes",
    malnutrition_mhyn == "NO" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # Immuno drugs
  mutate(immuno_d = case_when(
    immno_cmtrt == "Yes" ~ "Yes",
    immno_cmtrt == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  # Vulnerable Tx
  mutate(tx = case_when(
    vulnerable_transplant == "Yes" ~ "Yes",
    vulnerable_transplant == "No" ~ "No",
    TRUE ~ "No"
  )) %>% 
  select(-aidshiv_mhyn, -malnutrition_mhyn, -immno_cmtrt, -vulnerable_transplant)

immunosuppession_covid_gram$immuno_gram <- rowSums(immunosuppession_covid_gram == "Yes")

immunosuppession_covid_gram = immunosuppession_covid_gram %>% 
  mutate(immuno_gram = case_when(
    immuno_gram >= 1 ~ 1,
    TRUE ~ 0
  )) %>% 
  select(subjid, immuno_gram)


covid_gram_score = cs %>% 
  left_join(comorbid_gram) %>% 
  left_join(immunosuppession_covid_gram) %>% 
  select(-hypertension, -heart_disease, -renal, -liver, -lung, -diabetes) %>% 
  select(subjid, comorbid_gram, immuno_gram) %>% 
  rename(comorbid_gram2 = comorbid_gram)

covid_gram_score$comorbid_gram <- covid_gram_score$comorbid_gram2 + covid_gram_score$immuno_gram

cs = cs %>% 
  left_join(covid_gram_score) %>% 
  select(-comorbid_gram2, -immuno_gram) %>% 
  ff_relabel_df(cs)


## ----Add in combined_outcome, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------
# Create new tibble to avoid recoding renal_mhyn in original dataset
# cs2 = cs
# 
# cs2$renal_mhyn[is.na(cs2$renal_mhyn)] = "No"
# 
# cs2 = cs2 %>% 
#   select(subjid, mort, any_escalation, any_rrt, renal_mhyn, any_icu, any_invasive,
#          any_noninvasive, any_inotrope) %>% 
#   mutate(combined_outcome = case_when(
#     mort == "1" | any_escalation == "Yes" ~ "1",
#     is.na(mort) & is.na(any_escalation) ~ NA_character_,
#     TRUE ~ "0") %>% 
#       factor(),
#     
#     on_rrt = case_when(
#       any_rrt == "Yes" & renal_mhyn == "Yes" & any_icu == "No" & any_invasive == "NO" &
#         any_noninvasive == "No" & any_inotrope == "No" & mort != "1" ~ "Yes",
#       TRUE ~ "No")
#   ) %>% 
#   filter(on_rrt != "Yes") %>% 
#   select(-on_rrt) %>% 
#   # Excludes 49 patients likely RRT dependent before admission
#   
#   # Recode combined_outcome factor
#   # mutate(combined_outcome = fct_recode(combined_outcome, 
#   #                                      "No intervention or death" = "0",
#   #                                      "Intervention and/or death" = "1")) %>% 
#   select(subjid, combined_outcome)
# 
# 
# 
# # Join with 
# cs = cs %>% 
#   left_join(cs2, by = "subjid")
# 
# rm(cs2)

na_decision = FALSE

## ----Apply scores, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------
cs = cs %>% 
  
  # Isaric score (clinical parameters plus urea and crp)
  isaric_min(age = age, sex = sex, comorbid = no_comorbid, rr = rr_vsorres,
             spo2 = oxy_vsorres, gcs = daily_gcs_vsorres, 
             bun = daily_bun_lborres, crp = daily_crp_lborres,
             output = c("df_vector"), na_to_zeros = na_decision) %>%
  
  # Isaric score (clinical parameters only)
  isaric_clin(age = age, sex = sex, comorbid = no_comorbid, rr = rr_vsorres,
              spo2 = oxy_vsorres, gcs = daily_gcs_vsorres,
              output = c("df_vector"), na_to_zeros = na_decision) %>% 
  
  # Isaric Cox score
  isaric_cox(age = age, sex = sex, comorbid = no_comorbid, rr = rr_vsorres,
             spo2 = oxy_vsorres, gcs = daily_gcs_vsorres, 
             bun = daily_bun_lborres, crp = daily_crp_lborres,
             # plts = daily_plt_lborres, creat = daily_creat_lborres, nlr = NLR,
             output = c("df_vector"), na_to_zeros = na_decision) %>% 
  
  # Isaric Cox clinical score
  isaric_cox_clinical(age = age, sex = sex, comorbid = no_comorbid, rr = rr_vsorres,
                      spo2 = oxy_vsorres, gcs = daily_gcs_vsorres, 
                      # bun = daily_bun_lborres, crp = daily_crp_lborres,
                      # plts = daily_plt_lborres, creat = daily_creat_lborres, nlr = NLR,
                      output = c("df_vector"), na_to_zeros = na_decision) %>% 
  
  # Isaric Cox bootstrap average
  isaric_cox_average(age = age, sex = sex, comorbid = no_comorbid, rr = rr_vsorres,
                     spo2 = oxy_vsorres, gcs = daily_gcs_vsorres, 
                     bun = daily_bun_lborres, crp = daily_crp_lborres, #plts = daily_plt_lborres,
                     # plts = daily_plt_lborres, creat = daily_creat_lborres, nlr = NLR,
                     output = c("df_vector"), na_to_zeros = FALSE)


## ----Filter recent 4 weeks, message=FALSE, warning=FALSE------------------------------------------------------------------------------------------------------

# Create mortality cohort
cs = cs %T>% 
  {removals <<- bind_rows(removals, tibble(distinct = n_distinct(.$subjid),
                                           nrow = nrow(.),
                                           label = "cs-before filter(subjid %in% keep_14)"))} %>% 
  filter(subjid %in% keep_14) %T>% 
  {removals <<- bind_rows(removals, tibble(distinct = n_distinct(.$subjid),
                                           nrow = nrow(.),
                                           label = "cs-after filter(subjid %in% keep_14)"))} %>% 
  # filter(! subjid %in% keep_14_28) %T>% 
  # {removals <<- bind_rows(removals, tibble(distinct = n_distinct(.$subjid),
  #                                          nrow = nrow(.),
  #                                          label = "cs-after filter(! subjid %in% keep_14_28)"))} %>% 
  filter(age >= 18) %T>% 
  {removals <<- bind_rows(removals, tibble(distinct = n_distinct(.$subjid),
                                           nrow = nrow(.),
                                           label = "cs-after filter(age >= 18)"))} %>% 
  filter(!is.na(mort)) %T>% 
  {removals <<- bind_rows(removals, tibble(distinct = n_distinct(.$subjid),
                                           nrow = nrow(.),
                                           label = "cs-after filter(!is.na(mort))"))} %>% 
  ff_relabel_df(cs) %>% 
  mutate(NLR = ff_label(NLR, "NLR"))


## ----Date information, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------

# Create new tibble name for cs here to provide patient count for those who meet inclusion criteria
cs_number = cs

# Create date where patients admitted on this date or after are excluded
inclusion_date <- cs %>%
  select(subjid, hostdat) %>% 
  group_by(subjid) %>% 
  slice(which.max(hostdat))

date_cutoff <- max(inclusion_date$hostdat) + days(1)
# 
# 
# ## ----Remove outliers, message=FALSE, warning=FALSE------------------------------------------------------------------------------------------------------------
# cs = cs %>% 
#   filter(oxy_vsorres > 50 & oxy_vsorres < 101 | is.na(oxy_vsorres)) %T>% 
#   {removals <<- bind_rows(removals, tibble(distinct = n_distinct(paste0(.$subjid, .$arm_n)),
#                                            nrow = nrow(.),
#                                            label = "cs-after filter(oxy_vsorres > 50 & oxy_vsorres < 101)"))} %>% 
#   filter(daily_plt_lborres < 2000 | is.na(daily_plt_lborres)) %T>% 
#   {removals <<- bind_rows(removals, tibble(distinct = n_distinct(paste0(.$subjid, .$arm_n)),
#                                            nrow = nrow(.),
#                                            label = "cs-after filter(daily_plt_lborres < 2000)"))} %>% 
#   filter(NLR >= 0 & NLR < 200 | is.na(NLR)) %T>% 
#   {removals <<- bind_rows(removals, tibble(distinct = n_distinct(paste0(.$subjid, .$arm_n)),
#                                            nrow = nrow(.),
#                                            label = "cs-after filter(NLR >= 0 & NLR < 200) "))} %>% 
#   filter(admission_diabp_vsorres >20 & admission_diabp_vsorres <200 | is.na(admission_diabp_vsorres)) %T>% 
#   {removals <<- bind_rows(removals, tibble(distinct = n_distinct(paste0(.$subjid, .$arm_n)),
#                                            nrow = nrow(.),
#                                            label = "cs-after filter(admission_diabp_vsorres >20 & admission_diabp_vsorres <200)"))} %>% 
#   filter(sysbp_vsorres >40 & sysbp_vsorres < 250 | is.na(sysbp_vsorres)) %>% 
#   ff_relabel_df(cs) %T>% 
#   {removals <<- bind_rows(removals, tibble(distinct = n_distinct(paste0(.$subjid, .$arm_n)),
#                                            nrow = nrow(.),
#                                            label = "cs-after filter(sysbp_vsorres >40 & sysbp_vsorres < 250)"))}
# 
# 
# 
# ## ----Upload derivation and test datasets, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------
# # This has previously been created and merging allocation with cs using file 'patient_allocation_list.rda'
# # 
# # cs = cs %>% 
# #   left_join(patient_allocation_list)


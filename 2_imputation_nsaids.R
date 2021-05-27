#Load in libraries
library(tidyverse)
library(MatchIt)
library(cobalt)
library(survival)
library(finalfit)
library(MatchThem)
library(mice)
options(scipen = 999) # for ease of reading

#Remove any additional objects as this script requires a lot of memory to run
rm(comorbid, comorbid_gram, comorbid_noheart, comorbid_surgisphere, cs, cs_number, cs_selected, immunosuppession_covid_gram,
   immunosuppession_surgisphere, inclusion_date, medication_out, surgisphere_score, removals, isaric_clin, isaric_cox,
   isaric_cox_average, isaric_cox_clinical, isaric_min, lasso_score_full, lasso_score_full2)

#source preparation script
source('1_prep_nsaids.R')
source('mimi_functions.R') #contains functions for handling mimi objects - note this may need editing as broom / tidyverse seems to change its outputs recently

#Here we add in the aki outcome data - here for completeness
aki_data # = read_rds()

#Define variables to match on, the latter don't include age so these can be included as exact matches
explanatory_to_match_on = c('age', 'sex', 'chrincard', 'chronicpul_mhyn', 'diabetes_fct', 'obesity_mhyn', 'renal_mhyn', 'rheumatologic_mhyn', 'dementia_mhyn')
explanatory_to_match_on2 = c('sex', 'chrincard', 'chronicpul_mhyn', 'diabetes_fct', 'obesity_mhyn', 'renal_mhyn', 'rheumatologic_mhyn', 'dementia_mhyn')
explanatory_to_match_on3 = c('sex', 'chrincard', 'chronicpul_mhyn', 'diabetes_fct', 'obesity_mhyn', 'renal_mhyn', 'dementia_mhyn')

#Make sure ages are whole numbers
topline_data = topline_data %>% 
  mutate(age = as.integer(age)) 

topline = topline %>% 
  mutate(age = as.integer(age)) 

# First generate complete case matches
matched_data = topline %>% 
  filter(!is.na(death)) %>% 
  mutate(news = as.factor(news)) %>% 
  select(subjid, nsaids_yn, explan_severity_2, explanatory_to_match_on) %>% 
  mutate(nsaids_yn = factor(nsaids_yn)) %>% 
  mutate(nsaids_0_1 = ifelse(nsaids_yn == 'NSAIDs', 0,1)) %>% 
  drop_na(explan_severity_2, explanatory_to_match_on)

matchit_data = matchit(as.formula(ff_formula('nsaids_0_1', explanatory_to_match_on)), method = 'nearest', ratio = 1, calliper = 0.2, data = matched_data)

matched_data = match.data(matchit_data)


# Now impute as otherwise throw lots of patients away from match
topline_m_impute = topline %>% 
  filter(!is.na(death)) %>% #Remove those with missing outcomes 
  filter(!is.na(nsaids_yn)) %>% #Remove those with missing NSAID exposure - i.e. can't impute this as doesn't make sense to
  mutate(nsaids_yn = factor(nsaids_yn)) %>% 
  ff_relabel_df(topline_data) %>%
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  select(nsaids_yn, explan_severity_2, explanatory_to_match_on, death) %>% 
  mutate(nsaids_0_1 = nsaids_yn %>% factor() %>% fct_relevel('No NSAIDs')) %>% 
  mutate(nsaids_0_1 = ifelse(nsaids_0_1 == 'No NSAIDs', 0, 1)) %>%
  # distinct(subjid, .keep_all = T) %>% 
  mice(m = 1, maxit = 1)#Increase this to at least 5 and 5 - set 1 and 1 currently just for exploration

matched_imputed_surv = matchthem(formula = ff_formula('nsaids_0_1', explanatory_to_match_on), method = 'nearest', datasets = topline_m_impute, exact_vars = 'age') #Nearest neighbout, with exact age matching, as this is the most strong predictor of covid death

#Make minimal and full models
fit_imputed = with(matched_imputed_surv, glm(death ~ nsaids_yn, family = 'binomial'))
fit_imputed = with(matched_imputed_surv, glm(death ~ nsaids_yn, family = 'binomial'))
fit_imputed2 = with(matched_imputed_surv, glm(death ~ nsaids_yn, family = 'binomial' + age.factor + sex))
fit_imputed3 = with(matched_imputed_surv, glm(death ~ nsaids_yn, family = 'binomial'+ age.factor + sex + chrincard + obesity_mhyn + chronicpul_mhyn + diabetes_fct + renal_mhyn + rheumatologic_mhyn + dementia_mhyn))



bal.plot(matched_imputed_surv)
bal.tab(matched_imputed_surv)

#Ibuprofen
topline_m_impute_ibu = topline_data %>% 
  filter(nsaids == 'No NSAIDs' | nsaids == 'Ibuprofen (or similar propionic acid derivative)') %>% 
  #mutate_at(vars(explanatory_to_match_on), factor) %>% f
  filter(!is.na(death)) %>% 
  filter(!is.na(nsaids_yn)) %>% 
  mutate(nsaids_yn = factor(nsaids_yn)) %>% 
  ff_relabel_df(topline_data) %>%
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  select(nsaids_yn, nsaids, explan_severity_2, explanatory_to_match_on, death) %>% 
  mutate(nsaids_0_1 = nsaids_yn %>% factor() %>% fct_relevel('No NSAIDs')) %>% 
  mutate(nsaids_0_1 = ifelse(nsaids_0_1 == 'No NSAIDs', 0,1)) %>%
  # distinct(subjid, .keep_all = T) %>% 
  mice(m = 1, maxit = 1)

matched_imputed_topline_ibu = matchthem(formula = ff_formula('nsaids_0_1', explanatory_to_match_on), method = 'nearest', datasets = topline_m_impute_ibu, exact_vars = 'age')

fit_imputed_ibu = with(matched_imputed_topline_ibu, glm(death ~ nsaids_yn, family = 'binomial'))

#COX-2
topline_m_impute_cox2 = topline_data %>% 
  filter(nsaids == 'No NSAIDs' | nsaids == 'COX-2 Inhibitor') %>% 
  #mutate_at(vars(explanatory_to_match_on), factor) %>% f
  filter(!is.na(death)) %>% 
  filter(!is.na(nsaids_yn)) %>% 
  mutate(nsaids_yn = factor(nsaids_yn)) %>% 
  ff_relabel_df(topline_data) %>%
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  select(nsaids_yn, nsaids, explan_severity_2, explanatory_to_match_on, death) %>% 
  mutate(nsaids_0_1 = nsaids_yn %>% factor() %>% fct_relevel('No NSAIDs')) %>% 
  mutate(nsaids_0_1 = ifelse(nsaids_0_1 == 'No NSAIDs', 0,1)) %>%
  # distinct(subjid, .keep_all = T) %>% 
  mice(m = 5, maxit = 5)

matched_imputed_topline_cox2  = matchthem(formula = ff_formula('nsaids_0_1', explanatory_to_match_on), method = 'nearest', datasets = topline_m_impute_cox2, exact_vars = 'age' )

fit_imputed_cox2  = with(matched_imputed_topline_cox2 , glm(death ~ nsaids_yn, family = 'binomial'))

topline_m_impute_r = topline_data %>% 
  filter(rheumatologic_mhyn == 'Yes') %>% 
  #mutate_at(vars(explanatory_to_match_on), factor) %>% 
  filter(!is.na(death)) %>% 
  filter(!is.na(nsaids_yn)) %>% 
  mutate(nsaids_yn = factor(nsaids_yn)) %>% 
  ff_relabel_df(topline_data) %>%
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>%
  select(nsaids_yn, explan_severity_2, explanatory_to_match_on2, death) %>% 
  mutate(nsaids_0_1 = nsaids_yn %>% factor() %>% fct_relevel('No NSAIDs')) %>% 
  mutate(nsaids_0_1 = ifelse(nsaids_0_1 == 'No NSAIDs', 0,1)) %>%
  mice(m = 5, maxit = 5)

matched_imputed_topline_r = matchthem(formula = ff_formula('nsaids_0_1', explanatory_to_match_on3), method = 'nearest', datasets = topline_m_impute_r, exact_vars = 'age' )

fit_imputed_r = with(matched_imputed_topline_r, glm(death ~ nsaids_yn, family = 'binomial'))
summary(pool(fit_imputed_r), exponentiate = T, conf.int = T, conf.level = 0.95) 


#For severity
#qSOFA
qsofa_imp_m = topline %>% 
  filter(!is.na(nsaids_yn)) %>% 
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  mutate(qsofa = as.numeric(as.character(qsofa))) %>% 
  mice_lm_permute(dependent_regression ='qsofa', list_to_match_on = explanatory_to_match_on, reference_level = 'No NSAIDs', treat_var = 'nsaids_yn', exacts = 'age' )

qsofa_imp_pool = with(qsofa_imp_m, lm(as.formula(ff_formula('qsofa', 'nsaids_yn'))))

bal.tab(qsofa_imp_m)

#isaric
isaric_imp_m = topline %>% 
  filter(!is.na(nsaids_yn)) %>% 
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  mutate(isaric_clin = as.numeric(as.character(isaric_clin))) %>% 
  mice_lm_permute(dependent_regression ='isaric_clin', list_to_match_on = explanatory_to_match_on, reference_level = 'No NSAIDs', treat_var = 'nsaids_yn', exacts = 'age' )

isaric_imp_pool = with(isaric_imp_m, lm(as.formula(ff_formula('isaric_clin', 'nsaids_yn'))))

pool(isaric_imp_pool)

#NEWS
news_imp_m = topline %>% 
  filter(!is.na(nsaids_yn)) %>% 
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  mutate(qsofa = as.numeric(as.character(news))) %>% 
  mice_lm_permute(dependent_regression ='news', list_to_match_on = explanatory_to_match_on, reference_level = 'No NSAIDs', treat_var = 'nsaids_yn', exacts = 'age' )

news_imp_pool = with(news_imp_m, lm(as.formula(ff_formula('news', 'nsaids_yn'))))


#HR
hr_vsorres_imp_m = topline %>% 
  filter(!is.na(nsaids_yn)) %>% 
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  mice_lm_permute(dependent_regression ='hr_vsorres', list_to_match_on = explanatory_to_match_on, reference_level = 'No NSAIDs', treat_var = 'nsaids_yn', exacts = 'age' )

hr_vsorres_imp_pool = with(hr_vsorres_imp_m, lm(as.formula(ff_formula('hr_vsorres', 'nsaids_yn'))))


#RR
rr_vsorres_imp_m = topline %>% 
  filter(!is.na(nsaids_yn)) %>% 
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  mice_lm_permute(dependent_regression ='rr_vsorres', list_to_match_on = explanatory_to_match_on, reference_level = 'No NSAIDs', treat_var = 'nsaids_yn', exacts = 'age' )

rr_vsorres_imp_pool = with(rr_vsorres_imp_m, lm(as.formula(ff_formula('rr_vsorres', 'nsaids_yn'))))

#SpO2
oxy_vsorres_imp_m = topline %>% 
  filter(!is.na(nsaids_yn)) %>% 
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  mice_lm_permute(dependent_regression ='oxy_vsorres', list_to_match_on = explanatory_to_match_on, reference_level = 'No NSAIDs', treat_var = 'nsaids_yn', exacts = 'age' )

oxy_vsorres_imp_pool = with(oxy_vsorres_imp_m, lm(as.formula(ff_formula('oxy_vsorres', 'nsaids_yn'))))

#sbp
sysbp_vsorres_imp_m = topline %>% 
  filter(!is.na(nsaids_yn)) %>% 
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  mice_lm_permute(dependent_regression ='sysbp_vsorres', list_to_match_on = explanatory_to_match_on, reference_level = 'No NSAIDs', treat_var = 'nsaids_yn', exacts = 'age' )

sysbp_vsorres_imp_pool = with(sysbp_vsorres_imp_m, lm(as.formula(ff_formula('sysbp_vsorres', 'nsaids_yn'))))

#dbp
admission_diabp_vsorres_imp_m = topline %>% 
  filter(!is.na(nsaids_yn)) %>% 
  mutate_if(is.factor, fct_recode, NULL = '(Missing)') %>% 
  mice_lm_permute(dependent_regression ='admission_diabp_vsorres', list_to_match_on = explanatory_to_match_on, reference_level = 'No NSAIDs', treat_var = 'nsaids_yn', exacts = 'age' )

admission_diabp_vsorres_imp_pool = with(admission_diabp_vsorres_imp_m, lm(as.formula(ff_formula('admission_diabp_vsorres', 'nsaids_yn'))))

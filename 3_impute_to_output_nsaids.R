library(MatchThem)
library(tidyverse)
library(finalfit)
library(survival)
options(scipen = 999) # for readability

source('mimi_functions.R') # contains functions for handling mimi objects - note this may need editing as broom / tidyverse seems to change its outputs recently


#Main models for mortality
nsaids_main = exp_mimi(fit_imputed)
nsaids_main2 = exp_mimi(fit_imputed2)
nsaids_main3 = exp_mimi(fit_imputed3)

#Main models for ibuprofen / cox2 vs no NSAID
ibu_main = exp_mimi(fit_imputed_ibu)
cox2_main = exp_mimi(fit_imputed_cox2)

#Physiological parameters
matched_severity = bind_rows(mimi_lm(qsofa_imp_pool, 'qSOFA'),
                             mimi_lm(news_imp_pool, 'NEWS'),
                             mimi_lm(hr_vsorres_imp_pool, 'Heart rate'),
                             mimi_lm(rr_vsorres_imp_pool, 'Respiratory rate'),
                             mimi_lm(oxy_vsorres_imp_pool, 'SpO2'),
                             mimi_lm(sysbp_vsorres_imp_pool, 'Systolic blood pressure'),
                             mimi_lm(admission_diabp_vsorres_imp_pool, 'Diastolic blood pressure'))

#Now format all together
matched_severity %>% 
  tibble() %>% 
  mutate_if(is.double, round, 3) %>% 
  mutate(effect = paste0(estimate, ' (', l_95, ' to ', u_95, ', p = ', p_value,')')) %>% 
  select(variable, effect)


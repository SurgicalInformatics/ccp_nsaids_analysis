library(tidyverse)
library(MatchIt)
library(MatchThem)

#Imputation
mice_permute = function(.data, dependent_regression, treat_var, list_to_match_on, reference_level = 'No', method = 'nearest', exacts){
  require(MatchThem)
  require(tidyverse)
  match_data_in = .data
  
  match_data_impute = match_data_in %>% 
    mutate(!!dependent_regression := factor(!!sym(dependent_regression))) %>% 
    filter(!is.na(!!sym(dependent_regression))) %>% 
    select(dependent_regression,treat_var, list_to_match_on) %>% 
    mutate(treat_var_0_1 = !!sym(treat_var) %>% factor() %>% fct_relevel(reference_level)) %>% 
    mutate(treat_var_0_1 = ifelse(treat_var_0_1 == reference_level, 0,1)) %>%
    mice(m = 5, maxit = 5)
  
  matched_imputed = matchthem(formula = ff_formula(dependent = 'treat_var_0_1', explanatory = list_to_match_on), method = method, datasets = match_data_impute, calliper = 0.5, exact_vars = exacts)
  
  return(matched_imputed)
}


mice_lm_permute = function(.data, dependent_regression, treat_var, list_to_match_on, reference_level = 'No', method = 'nearest', exacts){
  require(MatchThem)
  require(tidyverse)
  match_data_in = .data
  
  match_data_impute = match_data_in %>% 
    filter(!is.na(!!sym(dependent_regression))) %>% 
    select(dependent_regression,treat_var, list_to_match_on) %>% 
    mutate(treat_var_0_1 = !!sym(treat_var) %>% factor() %>% fct_relevel(reference_level)) %>% 
    mutate(treat_var_0_1 = ifelse(treat_var_0_1 == reference_level, 0,1)) %>%
    mice(m = 5, maxit = 5)
  
  matched_imputed = matchthem(formula = ff_formula(dependent = 'treat_var_0_1', explanatory = list_to_match_on), method = method, datasets = match_data_impute, calliper = 0.5, exact_vars = exacts)
  
  return(matched_imputed)
}

exp_mimi = function(x){ 
  sum_df = summary(pool(x), exponentiate = T, conf.int = T, conf.level = 0.95) %>% select(estimate, `2.5 %`, `97.5 %`, p.value) %>% rownames_to_column('variable') %>% 
    rename(l_95 = `2.5 %`,
           u_95 = `97.5 %`,
           p_value = p.value) %>% 
    filter(variable != '(Intercept)')
  
  bind_rows(tibble(variable = 'No NSAIDs', estimate = 1, l_95 = 1, u_95 = 1, p_value = 1) , 
            sum_df) %>% 
    mutate(variable = gsub('nsaids_yn', '', variable)) 
  
}


exp_mimi2 = function(x, label){ 
  sum_df = summary(pool(x), exponentiate = T, conf.int = T, conf.level = 0.95) %>% select(estimate, `2.5 %`, `97.5 %`, p.value) %>% rownames_to_column('variable') %>% 
    rename(l_95 = `2.5 %`,
           u_95 = `97.5 %`,
           p_value = p.value) %>% 
    filter(variable != '(Intercept)')
  
  bind_rows(tibble(variable = 'No NSAIDs', estimate = 1, l_95 = 1, u_95 = 1, p_value = 1) , 
            sum_df) %>% 
    mutate(variable = gsub('nsaids_yn', '', variable)) %>% 
    mutate(variable = paste0(label, ' - ', variable))
  
}



exp_mimi_explicit = function(x, var_1, var_level, label){ 
  sum_df = summary(pool(x), exponentiate = T, conf.int = T, conf.level = 0.95)%>% rownames_to_column('term') %>% select(term, estimate, `2.5 %`, `97.5 %`, p.value) %>% rename(variable = term) %>% 
    rename(l_95 = `2.5 %`,
           u_95 = `97.5 %`,
           p_value = p.value) %>% 
    filter(variable != '(Intercept)')
  
  bind_rows(tibble(variable = var_level, estimate = 1, l_95 = 1, u_95 = 1, p_value = 1) , 
            sum_df) %>% 
    mutate(variable = gsub(var_1, '', variable)) %>% 
    mutate(variable = paste0(label, ' - ', variable))
  
}


mimi_lm = function(x, label){ 
  summary(pool(x), conf.int = T, conf.level = 0.95) %>% select(estimate, `2.5 %`, `97.5 %`, p.value) %>% rownames_to_column('variable')%>% 
    rename(l_95 = `2.5 %`,
           u_95 = `97.5 %`,
           p_value = p.value) %>% 
    filter(variable != '(Intercept)') %>% 
    mutate(variable = label)
}

library(tidyverse)

# Mean age at TB incidence, globally

tb_burden_estimates <- read.csv("TB_burden_age_sex_2024-06-25.csv")
# https://www.who.int/teams/global-tuberculosis-programme/data
# WHO TB incidence estimates disaggregated by age group, sex and risk factor [>0.5Mb]



tb_ages <- tb_burden_estimates %>% filter(year == 2022, sex %in% c("m","f"), risk_factor == "all",
                                age_group %in% c("0-4", "0-14", "5-14", "15-24", "25-34", "35-44", "45-54", "55-64", "65plus")) %>%
                                group_by(age_group, country) %>% 
                                summarise(both_sexes = sum(best)) %>%
                                pivot_wider(values_from = "both_sexes", names_from = age_group) %>% 
                                mutate(`5-14` = `0-14` - `0-4`) %>%
                                select(-`0-14`) %>%
                                # move the `5-14` column to the second column position
                                select(country, `0-4`, `5-14`, everything()) %>%
                                summarise(across(where(is.numeric), sum))
age_medians <- c(2, 9.5, 19.5, 29.5, 39.5, 49.5, 59.5, 69.5)
sum(tb_ages[1:4])/sum(tb_ages)
sum(tb_ages[1:5])/sum(tb_ages)
# so the median is in the 5th age group, 35-44
sum(unlist(tb_ages %>% ungroup())*age_medians)/sum(unlist(tb_ages %>% ungroup()))
# and mean is ~39.
library(tidyverse)
library(boot)

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



# Healthy life expectancy at age 40, globally
lifetable <- read.csv("region_life_tables.csv")
lifetable %>% 
    filter(Period == 2019,
           Dim2ValueCode == "AGEGROUP_AGE40-44",
           Indicator == "ex - expectation of life at age x") %>%
    select(Location, Value)


# Estimated rate of TB incidence decline by country
tb_burden_estimates <- read.csv("TB_burden_countries_2024-06-25.csv")
trends <- tb_burden_estimates %>% filter(year > 2015) %>% 
       group_by(country) %>%
       filter(max(e_inc_100k)>100) %>%
       summarise(data = list(cur_data_all()),
                 models = list(lm(e_inc_100k ~ year, data = data[[1]])))
rates <- lapply(trends$models,
              # get the rate of change relative to the 2021 value
              function(model) coef(model)[2]/predict(model, data.frame(year = 2021)))
summary(unlist(rates))
# IQR -0-4% per year


# TB case fatality ratio by HIV status
incidence <- 10600000
hivpos_incidence <- 703000
hivneg_incidence <-  incidence - hivpos_incidence

hivneg_deaths <- 1380000
hivpos_deaths <- 187000

# https://www.who.int/teams/global-tuberculosis-programme/tb-reports/global-tuberculosis-report-2022/tb-disease-burden/2-2-tb-mortality
# https://www.who.int/teams/global-tuberculosis-programme/tb-reports/global-tuberculosis-report-2022/tb-disease-burden/2-1-tb-incidence


(cfr_hivpos <- hivpos_deaths/hivpos_incidence)
(cfr_hivneg <- hivneg_deaths/hivneg_incidence)

cfr_hivpos/cfr_hivneg

# Covariance, mortality and duration
# If HIV+ = 2x risk of death, and 1/2 the duration, with cvs of 1, ...

# What if immune status is a latent class, and HIV is just strongly associated with it?
immune_status <- rnorm(n=10000, 0, 5)
# rank order immune status
immune_status <- immune_status[order(immune_status)]
hiv_status <- rbinom(n=10000, size=1, prob= 0.4*inv.logit(-immune_status))
mean(hiv_status)
plot(inv.logit(-immune_status), hiv_status)
cor(immune_status, hiv_status)
# sample durations 
duration <- rgamma(10000, 1,1)
# order durations with some error
duration <- duration[order(duration + rnorm(10000, 0, 1))]
mean(duration[hiv_status == 1])/mean(duration[hiv_status == 0])
plot(immune_status, duration)
summary(duration)
summary(immune_status)

# assign mortality, strongly correlated with immune status
mortality <- rbinom(10000, 1, 0.4*inv.logit(-0.5*immune_status -0.5))
summary(mortality)
mean(mortality[hiv_status == 1])/mean(mortality[hiv_status == 0])

cor(duration, mortality) * sd(duration) * sd(mortality)
cor(duration/mean(duration), mortality/mean(mortality)) * sd(duration/mean(duration)) * sd(mortality/mean(mortality))

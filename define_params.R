params <- list()

# First block of params, needed for averages

## TB morbidity
params$tb_symptom_duration <- c(3/12, 1/12, 5/12) #retrospective symptom reporting of ~120 days, not all as severe as basis for IHME DW estimate #nolint
params$tb_symptom_dw <- c(0.35, 0.2, 0.6) # IHME 0.33 for HIV-, 0.41 for HIV+

## TB mortality
params$tb_cfr <- c(0.12, 0.05, 0.2) # WHO Global Report 2023 mortality:incidence
params$tb_death_yearslost <- c(15, 5, 25) # average age at TB death ~55, remaining life expectancy 65-70 (reduced vs population average; https://www.who.int/data/gho/data/themes/mortality-and-global-health-estimates/ghe-life-expectancy-and-healthy-life-expectancy), but left skew on remaining life expectancy at time of TB death. # nolint
params$discounting_rate <- c(0.03, 0.0, 0.07) # determines weight of future post-TB and transmission DALYs

## Sequelae
params$posttb_symptom_duration <- c(25, 10, 35) # total post-TB survival. Average age ~40 and life expectancy ~65?.
params$posttb_symptom_dw <- c(0.035, 0.01, 0.10) # disability weight averaged over all post-TB survival, 0.035 per Menzies et al 2021. Corresponds to 20% prevalence of mild (DW 0.02), 5% of moderate (DW 0.2), and 5% of severe (DW 0.4) COPD, for example. # nolint
params$posttb_cfr <- c(0.05, 0.02, 0.15) # menzies paper, average 1.14x RR of mortality ==> ~.14/1.14=12% of TB survivors die of post-TB sequelae. But this seems high. # nolint
params$posttb_death_yearslost <- c(10, 5, 20) #  occurs preferentially among older TB patients (like TB deaths) but on ~ 5 year delay
params$posttb_death_timing <- c(5,1,10) #mean years to death, for those who die early of post-TB sequelae. Relevant to discounting only. 


# Transmission
params$downstream_cases <- c(0.9, 0.4, 2) # Estimated from Shrestha et al Step change model with one symptomatic or asympatomic case removed at present day. 
params$downstream_timing <- c(8, 2, 20) # Median time to downstream cases. Will apply time series discounting to exponential decay. 


# Second block of params, needed for time course

## Accrual before detectability

params$predetection_mm <- c(0.1,0,0.3) # proportion of morbidity, mortality, and post-TB disease DALYs that accrue before becomes detectable by screening algorithm
params$predetection_transmission <- c(0.05,0,0.1) # proportion of transmission that accrues before becomes detectable by screening algorithm

## Accrual after routine diagnosis

params$postrx_mm <- c(0.2,0.1,0.5) # proportion of morbidity, mortality, and post-TB disease DALYs that accrue after routine detection (must be less than proportion that die or resolve before detection)
params$postrx_transmission <- c(0.1,0.05,0.3) # proportion of transmission that occurs after routine detection

## Timing of accrual within detectable period
# Could model a correlation between MM accrual rate during the detectable period and how close a person is to getting routinely diagnosed.
# I.e., r(t) = r0 + r1 * t , if linear, or r(t) = r0 + r1 * t^2, if quadratic.
# We could parametrize this as a proportion of the total accrual that occurs in the first half of the detectable period (f), and the power of the above relationship (p). 


params$accrual_first_half_mm <- c(0.25, 0, 0.5) # proportion of morbidity, mortality, and post-TB disease that accrue during 2nd half of detectable period
params$accrual_first_half_transmission <- c(0.33, 0.1, 0.5)
params$accrual_power_mm <- c(3, 2, 4)
params$accrual_power_transmission <- c(1.5, 1, 2) 


# Third block of params, needed for differences beween average and detected cases

params$covariance_mortality_duration <- c(-0.1, -0.5, 0) # E.g. -0.5 could correspond to CV of 1 for both, and correlation of -0.5. Or CV of 1 for duration and 2 for mortality, and correlation of -0.25. If HIV+ = 2x risk of death, and 1/2 the duration, with cvs of 1, ...
params$covariance_transmission_duration <- c(1, 0.5, 1.5)


paramdf <- as.data.frame(t(as.data.frame(params)))
colnames(paramdf) <- c("mid", "low", "high")
paramdf$param <- rownames(paramdf)

write.csv(paramdf, file="DALY_model_param_values.csv")

midpoint_estimates <- as.list(paramdf$mid); names(midpoint_estimates) <- rownames(paramdf)


# Make a lookup table for the parameters and nice text names

nicenames <- list(
    "tb_symptom_duration" = "TB symptom duration",
    "tb_symptom_dw" = "TB disability weight",
    "tb_cfr" = "TB case fatality ratio",
    "tb_death_yearslost" = "Years of life lost per TB death",
    "discounting_rate" = "Annual discounting rate",
    "posttb_symptom_duration" = "Post-TB duration",
    "posttb_symptom_dw" = "Post-TB disability weight",
    "posttb_cfr" = "Post-TB case fatality ratio",
    "posttb_death_yearslost" = "Years of life lost per post-TB death",
    "posttb_death_timing" = "Post-TB death timing",
    "downstream_cases" = "Downstream cases",
    "downstream_timing" = "Downstream case timing",
    "predetection_transmission" = "Pre-detectability transmission",
    "predetection_mm" = "Pre-detectability morbidity and mortality",
    "postrx_transmission" = "Post-diagnosis transmission",
    "postrx_mm" = "Post-diagnosis morbidity and mortality",
    "accrual_first_half_mm" = "Early-detectable-period morbidity and mortality",
    "accrual_first_half_transmission" = "Early-detectable-period transmission",
    "accrual_power_mm" = "Late detectable period power for morbidity and mortality",
    "accrual_power_transmission" = "Late detectable period power for transmission",
    # "duration_cv" = "Disease duration heterogeneity",
    # "mortality_cv" = "Mortality risk heterogeneity",
    # "transmission_cv" = "Transmission heterogeneity",
    "covariance_mortality_duration" = "Mortality-duration covariance",
    "covariance_transmission_duration" = "Transmission-duration covariance")

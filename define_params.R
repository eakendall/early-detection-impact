params <- list()

# First block of params, needed for averages

## TB morbidity
params$tb_symptom_duration <- c(2/12, 1/12, 3/12) #retrospective symptom reporting of ~90 days, not all as severe as basis for IHME DW estimate #nolint
params$tb_symptom_dw <- c(0.35, 0.2, 0.6) # IHME 0.33 for HIV-, 0.41 for HIV+

## TB mortality
params$tb_cfr <- c(0.12, 0.05, 0.2) # WHO Global Report 2023 mortality:incidence
params$tb_death_yearslost <- c(15, 5, 30) # average age at TB death ~55, remaining life expectancy 65-70 (reduced vs population average; https://www.who.int/data/gho/data/themes/mortality-and-global-health-estimates/ghe-life-expectancy-and-healthy-life-expectancy), but left skew on remaining life expectancy at time of TB death. # nolint
params$discounting_rate <- c(0.03, 0.0, 0.07) # determines weight of future post-TB and transmission DALYs

## Sequelae
params$posttb_symptom_duration <- c(25, 10, 40) # total post-TB survival. Average age ~40 and life expectancy ~65?.
params$posttb_symptom_dw <- c(0.035, 0.01, 0.10) # disability weight averaged over all post-TB survival, per Menzies et al 2021. Corresponds to 20% prevalence of mild (DW 0.02), 5% of moderate (DW 0.2), and 5% of severe (DW 0.4) COPD, for example. # nolint
params$posttb_cfr <- c(0.05, 0.02, 0.15) # menzies paper, average 1.14x RR of mortality ==> ~.14/1.14=12% of TB survivors die of post-TB sequelae. But this seems high. # nolint
params$posttb_death_yearslost <- c(15, 5, 30) # similar to YLL from TB deaths, beacuse occur on a delay but also in people who were younger when they developed TB?
params$posttb_death_timing <- c(5,1,10) #mean years to death, for those who die early of post-TB sequelae. Relevant to discounting only. 


# Transmission
params$downstream_cases <- c(0.8, 0.5, 3) # should already be adjusted for temporal discounting. Estimated from Shrestha et al Step change model with one symptomatic or asympatomic case removed at present day. 


# Second block of params, needed for time course

## Timing of accrual within detectable period

params$second_half_vs_first_mm <- c(0.67, 0.5, 0.75 ) # proportion of morbidity, mortality, and post-TB disease that accrue during 2nd half of detectable period
params$second_half_vs_first_transmission <- c(0.67, 0.5, 0.75)

## Accrual before detectability

params$predetection_mm <- c(0.1,0,0.3) # proportion of morbidity, mortality, and post-TB disease DALYs that accrue before becomes detectable by screening algorithm
params$predetection_transmission <- c(0.05,0,0.1) # proportion of transmission that accrues before becomes detectable by screening algorithm

## Accrual after routine diagnosis

params$postrx_mm <- c(0.25,0,0.8) # proportion of morbidity, mortality, and post-TB disease DALYs that accrue after routine detection
params$postrx_transmission <- c(0.1,0,0.5) # proportion of transmission that occurs after routine detection



# Third block of params, needed for differences beween average and detected cases

params$duration_cv <- c(0.5,0.25, 1) # we will simulate gamma distribution of durations, with this coefficient of variation

params$duration_tbdeath_covarying_cv <- c(-0.5, -1, 0) # in the direction of correlation, how spread is tbdeath risk relative to duration (and is the correlation positive or negative). A value of 1 or -1 means positive or negative corelations, respectively, with the average mortaity at a given duration equal to the duration. Less than 1 means less spread, >1 means more spread than duration. 
params$duration_transmission_covarying_cv <- c(1, 0.5, 1.5)


paramdf <- as.data.frame(t(as.data.frame(params)))
colnames(paramdf) <- c("mid","low","high")
paramdf$param <- rownames(paramdf)

write.csv(paramdf, file="DALY_model_param_values.csv")

midpoint_estimates <- as.list(paramdf$mid); names(midpoint_estimates) <- rownames(paramdf)

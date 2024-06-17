params <- list()

# First block of params, needed for averages

## TB morbidity
params$tb_symptom_duration <- c(0.25, 2/12, 1)
params$tb_symptom_dw <- c(0.4, 0.3, 0.6)

## TB mortality
params$tb_cfr <- c(0.12, 0.05,0.2)
params$tb_death_yearslost <- c(20,10,35)
params$discounting_rate <- c(0.05,0.03,0.07)

## Sequelae
params$posttb_symptom_duration <- c(25,10,35)
params$posttb_symptom_dw <- c(0.1, 0.05, 0.2)
params$posttb_cfr <- c(0.05, 0.02, 0.10) # look at menzies paper
params$posttb_death_timing <- c(5,1,10) #mean years to death, for those who die early of post-TB sequelae
params$posttb_death_yearslost <- c(15,10,30) 

# Transmission
params$downstream_cases <- c(0.8, 0.5, 3) # should already be adjusted for temporal discounting


# Second block of params, needed for time course

## Timing of accrual within detectable period

params$second_half_vs_first_mm <- c(0.67, 0.5, 0.75 ) # proportion of morbidity, mortality, and post-TB disease that accrue during 2nd half of detectable period
params$second_half_vs_first_transmission <- c(0.67, 0.5, 0.75 ) 
params$resolving_detectable <- c(0.2,0.1,0.3) # proportion that will end in spontaneous resolution, after becoming detectable

## Accrual before detectability

params$predetection_mm <- c(0.1,0,0.3) # proportion of morbidity, mortality, and post-TB disease DALYs that accrue before becomes detectable by screening algorithm
params$predetection_transmission <- c(0.05,0,0.1) # proportion of transmission that accrues before becomes detectable by screening algorithm

## Accrual after routine diagnosis

params$postrx_mm <- c(0.25,0,1) # proportion of morbidity, mortality, and post-TB disease DALYs that accrue after routine detection
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

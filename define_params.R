params <- list()

params$tb_symptom_duration <- c(0.25, 2/12, 1)
params$tb_symptom_dw <- c(0.4, 0.3, 0.6)
params$tb_cfr <- c(0.12, 0.05,0.2)
params$tb_death_yearslost <- c(20,10,35)

params$posttb_symptom_duration <- c(25,10,35)
params$posttb_symptom_dw <- c(0.1, 0.05, 0.2)
params$posttb_cfr <- c(0.05, 0.02, 0.10) # look at menzies paper
params$posttb_death_timing <- c(5,1,10) #mean years to death, for those who die early of post-TB sequelae
params$posttb_death_yearslost <- c(15,10,30) 
params$downstream_cases <- c(0.8, 0.5, 3) # should already be adjusted for temporal discounting

params$discounting_rate <- c(0.05,0.03,0.07)

params$second_half_vs_first_mm <- c(0.67, 0.5, 0.75 ) # proportion of morbidity, mortality, and post-TB disease that accrue during 2nd half of detectable period
params$second_half_vs_first_transmission <- c(0.67, 0.5, 0.75 ) 
params$resolving_detectable <- c(0.2,0.1,0.3) # proportion that will end in spontaneous resolution, after becoming detectable
params$predetection_mm <- c(0.1,0,0.3) # proportion of morbidity, mortality, and post-TB disease DALYs that accrue before becomes detectable by screening algorithm
params$predetection_transmission <- c(0.05,0,0.1) # proportion of transmission that accrues before becomes detectable by screening algorithm
params$postrx_mm <- c(0.25,0,1) # proportion of morbidity, mortality, and post-TB disease DALYs that accrue after routine detection
params$postrx_transmission <- c(0.1,0,0.5) # proportion of transmission that occurs after routine detection


params$duration_cv <- c(0.5,0.25, 1) # will simulate gamma distribution of durations
params$duration_tbdeath_multiplier <- c(-0.5, -1, 0) # If case A has (1+a) times the average duration (and thus 1+a times the average risk of detection) in absence of ACF, what is its relative risk of death vs the average, in terms of a? For positive a this works fine, e.g. a value of 1 would mean 2x duration -> 2x mortality,  value of 0.5 would mean 2x duration -> 1.5x mortality, and value of 2 would mean 2x duration -> 3x mortality. For negatve a, e.g. -0.5, it means 2x duration -> 0.5x mortality, and -1 means 2x duration -> 0 mortality (so there a lower bound of -1).


params$duration_transmission_multiplier <- c(1, 0.5, 1.5) # Relationship between of total duration (in absence of ACF) and total transmission = 1?

paramdf <- as.data.frame(t(as.data.frame(params)))
colnames(paramdf) <- c("mid","low","high")
paramdf$param <- rownames(paramdf)

write.csv(paramdf, file="DALY_model_param_values.csv")

midpoint_estimates <- as.list(paramdf$mid); names(midpoint_estimates) <- rownames(paramdf)

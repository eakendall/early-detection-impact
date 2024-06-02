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


params$duration_cv <- c(0.5,0.25,0.75) # will simulation log-normal distribution of durations
params$duration_tbdeath_multiplier <- c(-0.5, -1, 0) # Relationship between total duration (in absence of ACF) and mortality risk = -0.5?
params$duration_transmission_multiplier <- c(1, 0.5, 1.5) # Relationship between of total duration (in absence of ACF) and total transmission = 1?
# Will map duration to DALY multipliers (relative to average DALYs), as quantiles of log-normal distributions. 
# I.e. multiplier of 1 means distributions have same sd and direction, 2 means DALY distribution has twice the sd (cv) and same direction,
# and -0.5 means DALY distribution has half the sd and we invert the mapping, i.e. nth quantile maps to (1-n)th. 

paramdf <- as.data.frame(t(as.data.frame(params)))
colnames(paramdf) <- c("mid","low","high")
paramdf$param <- rownames(paramdf)

write.csv(paramdf, file="DALY_model_param_values.csv")

midpoint_estimates <- as.list(paramdf$mid); names(midpoint_estimates) <- rownames(paramdf)

library(tidyverse)
if(!exists("midpoint_estimates")) source("define_params.R")

#### Per average case ####

dalys_per_average_case <- function(estimates = midpoint_estimates, plot = FALSE)
{
  with(estimates, {

  # TB morbidity DALYs
    morbidity_average_case <- tb_symptom_duration*tb_symptom_dw
    mortality_average_case <- tb_cfr*sum((1-discounting_rate)^(0:(tb_death_yearslost-1)))
  # Post-TB DALYs
    sequelae_average_case <- sum((1-discounting_rate)^(0:(posttb_symptom_duration-1)))*posttb_symptom_dw + 
      posttb_cfr*sum((1-discounting_rate)^(posttb_death_timing:(posttb_death_timing + posttb_death_yearslost))) 
  
  # Downstream case DALYs
  transmission_average_case <- downstream_cases*(mortality_average_case + morbidity_average_case + sequelae_average_case)
  
  # output a dataframe with the different components, classified as "cumulative" (vs averted) dalys for hte "average" (vs detected) case
  averages = tibble(
    "transmission" = transmission_average_case,
    "morbidity" = morbidity_average_case,
    "mortality" = mortality_average_case,
    "sequelae" = sequelae_average_case) %>% 
    mutate(cumulative_or_averted = "cumulative",
           average_or_detected = "average") %>%
    pivot_longer(cols = c(transmission, morbidity, mortality, sequelae))
  
  if(plot)
    print(plot_averages(averages))
    
  return(averages)

  })
}


####  Time weighting of DALYs within a case --> proportion avertible through early detection ####
# Of the symptoms, mortality risk, and sequelae accrual of a given case, (and then the same for transmission,)
# what proportion can't be prevented through screening because it occurs before they're detectable or after they're routinly diagnoses,
# and how much of the rest is in the second half of the disease course vs the first? 
# Assuming a linearly increase over time for those who eventually get treated, and a symmetric rise and fall for those who spontaneously resolve. 

within_case_avertible <- function(estimates = midpoint_estimates)
{
  with(estimates, {
    
    average_proportion_avertible_mm <- (1-predetection_mm - postrx_mm)*
      (0.5*resolving_detectable + second_half_vs_first_mm*(1-resolving_detectable)) 

    average_proportion_avertible_transmission <- (1-predetection_transmission - postrx_transmission)*
      (0.5*resolving_detectable + second_half_vs_first_transmission*(1-resolving_detectable))

      return(list(
        "average_proportion_avertible_mm" = average_proportion_avertible_mm,
        "average_proportion_avertible_transmission" = average_proportion_avertible_transmission))
  })
}



within_case_cumulative_and_averted <- function(estimates = midpoint_estimates)
{
  averages <- dalys_per_average_case(estimates)
  proportions <- within_case_avertible(estimates)
  averages_and_avertibles <- rbind(averages, 
                                   averages %>%
             mutate(cumulative_or_averted = "averted",
             value = value * case_when(name=="transmission" ~ proportions$average_proportion_avertible_transmission,
                                      TRUE ~ proportions$average_proportion_avertible_mm)))
  
  return(averages_and_avertibles)
  
}

within_case_cumulative_and_averted()


  # Utility function because we'll need to get the mu and sigma parameters for lognormal distributions with a desired mean and sd:
  solve_for_log_normal_parameters <- function(mean, sd)
  {  sigma2 = log(sd^2/mean^2 + 1)
  mu = log(mean) - sigma2/2
  return(list("mu"=mu, "sigma2"=sigma2))
  }
  
  # Gamma distribution version of the above function: for specified sd and mean,  get the shape and scale parameters
  solve_for_gamma_parameters <- function(mean, sd)
  {
    shape = mean^2/sd^2
    scale = sd^2/mean
    return(list("shape"=shape, "scale"=scale))
  }

#### Differences between detected and not detected cases ####
between_case_differences <- function(estimates = midpoint_estimates,
                                      N = 1000000) # found I needed >100k to get stable estimates
{

  # Now simulate durations ~ probabilities of detection for N incident cases:
  qs <- runif(N, 0, 1)

  with(estimates, {

    relative_durations <- qgamma( p=qs, 
                                  shape=(solve_for_gamma_parameters(1, duration_cv))[['shape']], 
                                  scale=(solve_for_gamma_parameters(1, duration_cv))[['scale']])
  
    # and for each, get the corresponding relative-DALY multipliers for mortality and transmission:
    # (don't need noise here since we'll just be averaging)
    relative_dalys_mortality <- qgamma(p = `if`(duration_tbdeath_multiplier>0, qs, 1-qs),
                                       shape=(solve_for_gamma_parameters(1, duration_cv * abs(duration_tbdeath_multiplier)))[['shape']],
                                       scale=(solve_for_gamma_parameters(1, duration_cv * abs(duration_tbdeath_multiplier)))[['scale']])
    relative_dalys_transmission <- qgamma(p = `if`(duration_transmission_multiplier>0, qs, 1-qs),
                                           shape=(solve_for_gamma_parameters(1, duration_cv * abs(duration_transmission_multiplier)))[['shape']],
                                            scale=(solve_for_gamma_parameters(1, duration_cv * abs(duration_transmission_multiplier)))[['scale']])
    
    # And simulate an ACF sample and their relative mortality and transmission DALYs, relative to the averages
    ACF_cases <- sample(1:N, size=N, replace = T, prob = relative_durations)
    avertible_mortality_multiplier_detected <- mean(relative_dalys_mortality[ACF_cases])
    avertible_transmission_multiplier_detected <- mean(relative_dalys_transmission[ACF_cases])

    return(list(
      "avertible_mortality_multiplier_detected" = avertible_mortality_multiplier_detected,
      "avertible_transmission_multiplier_detected" = avertible_transmission_multiplier_detected))
    })
}

#### Put it all together #####
daly_estimator <- function(estimates = midpoint_estimates)
{
  
  cumulativerows <- within_case_cumulative_and_averted(estimates) 
  detectedrows <- within_case_cumulative_and_averted(estimates) %>% mutate(
    average_or_detected = "detected",
    value = value * case_when(name=="transmission" ~ between_case_differences(estimates)$avertible_transmission_multiplier_detected,
                              name=="mortality"~ between_case_differences(estimates)$avertible_mortality_multiplier_detected,
                              TRUE ~ 1))
  
  return(rbind(cumulativerows, detectedrows))
}

daly_estimator(midpoint_estimates)


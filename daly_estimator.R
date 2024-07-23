library(tidyverse)
if(!exists("midpoint_estimates")) source("define_params.R")

#### Per average case ####

dalys_per_average_case <- function(estimates = midpoint_estimates, plot = FALSE)
{
  with(estimates, {

  # TB morbidity DALYs
    morbidity_average_case <- tb_symptom_duration*tb_symptom_dw
    mortality_average_case <- tb_cfr * sum((1 - discounting_rate)^(0:(tb_death_yearslost - 1)))

  # Post-TB DALYs
    sequelae_average_case <- sum((1 - discounting_rate)^(0:(posttb_symptom_duration - 1))) * posttb_symptom_dw +
      posttb_cfr * sum((1 - discounting_rate)^(posttb_death_timing:(posttb_death_timing + posttb_death_yearslost)))
  
  # Downstream case DALYs
  transmission_average_case <-
    downstream_cases *
    (mortality_average_case + morbidity_average_case + sequelae_average_case) *
    ((1 - discounting_rate)^(downstream_timing))
  
  
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
      (second_half_vs_first_mm) 

    average_proportion_avertible_transmission <- (1-predetection_transmission - postrx_transmission)*
      (second_half_vs_first_transmission)

      return(list(
        "average_proportion_avertible_mm" = average_proportion_avertible_mm,
        "average_proportion_avertible_transmission" = average_proportion_avertible_transmission))
  })
}



within_case_cumulative_and_averted <- function(averages = NULL, proportions = NULL, estimates = midpoint_estimates)
{
  if(missing(averages)) averages <- dalys_per_average_case(estimates)
  if(missing(proportions)) proportions <- within_case_avertible(estimates)
  averages_and_avertibles <- rbind(averages,
                                   averages %>%
             mutate(cumulative_or_averted = "averted",
             value = value * case_when(name == "transmission" ~ proportions$average_proportion_avertible_transmission,
                                       TRUE ~ proportions$average_proportion_avertible_mm)))
  
  return(averages_and_avertibles)
  
}



#### Differences between detected and not detected cases ####
between_case_differences <- function(estimates = midpoint_estimates)
{
  # New plan: If we specify the covariance (which includes the correlation and the coefficients of variation), then things are simple. But users need to know how to interpret covariance, so we can provide an illustration for specified duration_cv and correlation. 
  # covariance = correlation * duration_cv * mortality_cv (or transmission_cv)

  # We're defining these on relative scales, so among all incident TB, E[D] = E[M] = E[T] = 1.
  # We have joint distributions f[D,T] and f[D,M]. 
  # The covariance between D and T is integral dD dT (D-E[D])(T-E[T]) f[D,T] = integral dD dT (f D T - f D - f T + f) = 
  # integral dD dT f D T - integral dT (integral dD D f) - integral dD (integral dT T f) + integral dD dT f  = 
  # integral dD dT f D T - integral dT (E[D]) - integral dD (E[T]) + 1  = 
  # integral dD dT f D T - integral dT 1 - integral dD 1 + 1  = 
  # integral dD dT f D T - 1 - 1 + 1 

  # And the quantity we want is the expected value of T when weighted by D, i.e. 
  # integral dD dT f D T = cov(D,T) + 1. 

  # So we can just use the covariance plus 1 as the expected value of T when sampling by D.

  return(list(
    "avertible_mortality_multiplier_detected" = covariance_mortality_duration + 1,
    "avertible_transmission_multiplier_detected" = covariance_transmission_duration + 1))
}




#### Put it all together #####
daly_estimator <- function(within_case = NULL,
                           between_case = NULL,
                           estimates = midpoint_estimates)
{
  if (missing(within_case)) within_case <- within_case_cumulative_and_averted(estimates = estimates)
  if (missing(between_case)) between_case <- between_case_differences(estimates = estimates)
  
  cumulativerows <- within_case
  detectedrows <- within_case %>% mutate(
    average_or_detected = "detected",
    value = value * case_when(name == "transmission" ~ between_case$avertible_transmission_multiplier_detected,
                              name == "mortality"~ between_case$avertible_mortality_multiplier_detected,
                              TRUE ~ 1))
  
  return(rbind(cumulativerows, detectedrows))
}

# daly_estimator()


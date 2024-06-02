
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


#### Differences between detected and not detected cases ####
between_case_differences <- function(estimates = midpoint_estimates)
{
  # Utility function because we'll need to get the mu and sigma parameters for lognormal distributions with a desired mean and sd:
  solve_for_log_normal_parameters <- function(mean, sd)
  {  sigma2 = log(sd^2/mean^2 + 1)
  mu = log(mean) - sigma2/2
  return(list("mu"=mu, "sigma2"=sigma2))
  }
  
  # Now simulate durations ~ probabilities of detection for N incident cases:
  N <- 1000000 # found I needed >100k to get stable estimates
  qs <- runif(N, 0, 1)

  with(estimates, {

    relative_durations <- qlnorm(qs, unlist(solve_for_log_normal_parameters(1, duration_cv)))
    
    # and for each, get the corresponding relative-DALY multipliers for mortality and transmission:
    relative_dalys_mortality <- qlnorm(p = `if`(duration_tbdeath_multiplier>0, qs, 1-qs),
                                       unlist(solve_for_log_normal_parameters(1, duration_cv*abs(duration_tbdeath_multiplier))))
    relative_dalys_transmission <- qlnorm(p = `if`(duration_transmission_multiplier>0, qs, 1-qs),
                                          unlist(solve_for_log_normal_parameters(1, duration_cv*abs(duration_transmission_multiplier))))
    
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
  with(dalys_per_average_case(estimates), 
       with(within_case_avertible(estimates),
            with(between_case_differences(estimates), {
                                     
                transmission_averted <- average_proportion_avertible_transmission*transmission_average_case*avertible_transmission_multiplier_detected
                morbidity_averted <- average_proportion_avertible_mm * morbidity_average_case
                mortality_averted <- average_proportion_avertible_mm * mortality_average_case * avertible_mortality_multiplier_detected
                sequelae_averted <-  average_proportion_avertible_mm * sequelae_average_case
                
                total_averted <- transmission_averted + morbidity_averted + mortality_averted + sequelae_averted  
              
            
                return(list("total_averted" = total_averted, 
                            "transmission_averted" = transmission_averted,
                            "morbidity_averted" = morbidity_averted, 
                            "mortality_averted" = mortality_averted, 
                            "sequelae_averted" = sequelae_averted 
                            # "total_average_case" = downstream_dalys_average + tb_morbidity_dalys_average + tb_mortality_dalys_average + posttb_dalys_average,
                            # "transmission_average_case" = downstream_dalys_average,
                            # "morbidity_average_case" = tb_morbidity_dalys_average,
                            # "mortality_average_case" = tb_mortality_dalys_average,
                            # "sequelae_average_case" = posttb_dalys_average,
                            # "average_proportion_avertible_mm" = average_proportion_avertible_mm,
                            # "average_proportion_avertible_transmission" = average_proportion_avertible_transmission,
                            # "avertible_mortality_multiplier_detected" = case_dalys_mortality_multiplier,
                            # "avertible_transmission_multiplier_detected" = case_dalys_transmission_multiplier
                            ))
  })))
}

daly_estimator(midpoint_estimates)

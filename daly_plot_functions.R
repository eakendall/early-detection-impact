
# setwd("~/Google Drive/My Drive/DALY impact of ACF 2024/ACF-impact-Rshiny/")
source("daly_estimator.R")

library(gridExtra)
library(kableExtra)
library(ggpattern)
library(magick)
library(MASS)
library(conflicted)
conflicted::conflict_prefer(name = "select", winner = "dplyr")
conflicted::conflict_prefer(name = "filter", winner = "dplyr")

plot_averages <- function(output_dalys_per_average_case = NULL, 
                          input_first_block = midpoint_estimates, 
                          number_labels = TRUE, ymax = NULL)
{
 
  if(missing(output_dalys_per_average_case)) 
    output_dalys_per_average_case <- dalys_per_average_case()
   
  plot <- ggplot(output_dalys_per_average_case %>% 
                   filter(cumulative_or_averted == "cumulative",
                                   average_or_detected == "average")  %>% 
                   mutate(component = case_when(str_detect(name, "mortality") ~ "TB Mortality",
                                                str_detect(name, "morbidity") ~ "TB Morbidity",
                                                str_detect(name, "sequelae") ~ "Post-TB Sequelae",
                                                str_detect(name, "transmission") ~ "Transmission")) %>%
                   mutate(ordered_component = fct_relevel(component, "Transmission", "Post-TB Sequelae", "TB Mortality", "TB Morbidity")),
                 aes(x=average_or_detected, y=value, fill=ordered_component)) +  
    geom_col(position = "stack", width = 1) +  
    theme_minimal() + xlab("") + ylab("DALYs") + ggtitle("Total DALYS generated per average case") + 
    guides(fill="none") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(size = 16))
    
  
  if (number_labels) plot <- plot + geom_text(aes(label = paste0  (ordered_component, ", ",round(value,2))),
                                              position = position_stack(vjust = .5), size=5) else
                                                plot <- plot + geom_text(aes(label = ordered_component),  
                                                                         position = position_stack(vjust = .5),
                                                                         size = 5)

  if (!is.null(ymax)) plot <- plot + ylim(0,ymax)
                                              
 return(plot)
}

plot_detectable_proportion <- function(averages_plot, estimates=midpoint_estimates)
{
  if(missing(averages_plot)) averages_plot <- plot_averages(dalys_per_average_case(estimates), number_labels = F)
  
  with(estimates, {
    detectable_period_plot <- 
      averages_plot + geom_rect(aes(xmin=0.5, xmax=0.5 + predetection_mm, ymin=0, 
                                  ymax= sum(averages_plot$data %>% filter(name!="transmission") %>% select(value))), 
                              fill="gray", alpha=0.3) + 
      geom_rect(aes(xmin=0.5, xmax= 0.5 + predetection_transmission,
                    ymin= sum(averages_plot$data %>% filter(name!="transmission") %>% select(value))), 
                ymax= sum(averages_plot$data %>% select(value)), 
                fill="gray", alpha=0.3) + 
      annotate(geom="text", x=0.5 + predetection_mm/2, y=0.1, 
               label="Accrues before detectability", angle=90, hjust=0) +
      geom_rect(aes(xmin=1.5 - postrx_mm, xmax=1.5, ymin=0, 
                    ymax= sum(averages_plot$data %>% filter(name!="transmission") %>% select(value))), 
                fill="gray", alpha=0.3) + 
      geom_rect(aes(xmin=1.5 - postrx_transmission, xmax= 1.5,
                    ymin= sum(averages_plot$data %>% filter(name!="transmission") %>% select(value))), 
                ymax= sum(averages_plot$data %>% select(value)), 
                fill="gray", alpha=0.3) + 
      annotate(geom="text", x=1.5 - postrx_mm/2, y=0.1, 
               label="Accrues after routine diagnosis", angle=90, hjust=0) + 
      ggtitle("Avertible DALYs per average case") +
      ylab("Cumulative DALYs")
    return(detectable_period_plot)
  } )
}

# plot_detectable_proportion()


plot_time_course <- function(within_case = NULL, estimates = midpoint_estimates)
{
  if(missing(within_case)) within_case <- within_case_cumulative_and_averted(estimates = estimates)

  plotdata <- rbind(# total DALYs
                    within_case,
                    # during detectabnle period
                    within_case %>% filter(cumulative_or_averted=="cumulative") %>%
                      mutate(cumulative_or_averted="detectable",
                             value = value*case_when(name=="transmission" ~ (1-estimates$predetection_transmission -    estimates$postrx_transmission),
                                                     TRUE ~ (1-estimates$predetection_mm - estimates$postrx_mm))),
                    # before detectable period
                    within_case %>% filter(cumulative_or_averted=="cumulative") %>%
                      mutate(cumulative_or_averted="pre",
                             value = value*case_when(name=="transmission" ~ estimates$predetection_transmission,
                                                     TRUE ~ estimates$predetection_mm)),
                    # after detectable period
                    within_case %>% filter(cumulative_or_averted=="cumulative") %>%
                      mutate(cumulative_or_averted="post",
                             value = value*case_when(name=="transmission" ~ estimates$postrx_transmission,
                                                     TRUE ~ estimates$postrx_mm)))
                   
  
  # arbitrary detectable period width 1 (from 0 to +1) and 
  # transmission and mm DALYs accrual rates each r0 + r1*t^p, 
  # with f, p, and {r0, r1}(f,p) specific to mm or transmission.
  
  accrual_rate_params_transmission <- 
      accrual_rate_params(f = estimates$accrual_first_half_transmission, 
                                                          p = estimates$accrual_power_transmission)
  accrual_rate_params_mm <- 
      accrual_rate_params(f = estimates$accrual_first_half_mm,
                          p = estimates$accrual_power_mm)
  

  # plot accrual rates in the detectable periods as stacked polygons
  plotdata <- data.frame(x = seq(0,1,0.01)) %>% 
    mutate( y_transmission =   (accrual_rate_params_transmission[1] + 
                        accrual_rate_params_transmission[2] * x^estimates$accrual_power_transmission) * 
                        unlist((within_case %>% filter(cumulative_or_averted=="cumulative", 
                                                  name == "transmission") %>% 
                                           select(value))),
          y_morbidity = (accrual_rate_params_mm[1] +
                    accrual_rate_params_mm[2] * x^estimates$accrual_power_mm) * 
                    unlist((within_case %>% filter(cumulative_or_averted=="cumulative", 
                                                  name == "morbidity") %>% 
                                           select(value))),
          y_mortality = (accrual_rate_params_mm[1] + 
                  accrual_rate_params_mm[2] * x^estimates$accrual_power_mm) * 
                  unlist((within_case %>% filter(cumulative_or_averted=="cumulative", 
                                                name == "mortality") %>% 
                                      select(value))),
          y_sequelae = (accrual_rate_params_mm[1] +
                  accrual_rate_params_mm[2] * x^estimates$accrual_power_mm) * 
                  unlist((within_case %>% filter(cumulative_or_averted=="cumulative", 
                                                name == "sequelae") %>% 
                                      select(value)))) %>%
    # remove "y_" when creating component column of long format
    pivot_longer(-x, names_to = "component", names_pattern = "y_(.*)", values_to = "y")

  # plot the accrual rates as stacked polygons
  detectable_plot_period <- ggplot(data = plotdata %>% 
                                  # reverse order of levels so that transmission is on top
                                  mutate(component = fct_rev(fct_relevel(component, "morbidity", "mortality", "sequelae", "transmission"))),
                                   ) +
    geom_area(aes(x = x, y = y, group = component, fill = component), 
          position = position_stack()) +
    scale_fill_discrete(breaks = c("transmission", "sequelae", "mortality", "morbidity")) +
    xlab("Time during detectable period") +
    ylab("DALY accrual rate") +
    scale_x_continuous(labels = NULL, breaks = NULL, limits=c(-0.25, 1.25)) +
    scale_y_discrete(labels = NULL, breaks = NULL) +
    theme_minimal() +
    geom_vline(xintercept = 0.5) +
    ggtitle("Timing of DALY accrual") +
    theme(axis.text = element_text(size = 16),
          plot.title = element_text(size = 16))

  # Add gradient rectangle
  n <- 100
  ymax <- max(ggplot_build(detectable_plot_period)$data[[1]]$y)
  x_steps <- seq(from = -0.25, to = 0, length.out = n + 1)
  alpha_steps <- seq(from = 0.3, to = 0, length.out = n)
  rect_grad <- data.frame(xmin = x_steps[-1], 
                          xmax = x_steps[-n], 
                          alpha = rev(alpha_steps),
                          ymin = 0,
                          ymax = ymax)

  time_course_plot <- 
    detectable_plot_period +
    geom_rect(data=rect_grad, 
              aes(xmin=xmin, xmax=xmax,
                  ymin=ymin, ymax=ymax, 
                  alpha=alpha), fill="gray") + 
    guides(alpha = "none", fill = guide_legend(reverse = FALSE)) + 
    geom_vline(aes(xintercept = 1)) +
    annotate(geom = "text", x = -0.05, y = ymax/2,
             label = "before detectability", angle = 90, size=4, fontface = "italic") +
    annotate(geom = "text", x = 1.05, y = ymax / 2,
             label = "after routine detection", angle = 90, size=4, fontface = "italic") +
    ggtitle("Timing of DALY accrual") +
    theme(axis.text = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.position = "inside",
          legend.position.inside = c(.35,.7)) + 
    # at top of plot, add horizontal arrows from o to 1 and from 0 to -1
     annotate("segment", x = 0.52, y = ymax, xend = 0.95, yend = ymax, 
         linejoin = "mitre", linewidth = 5, color = "gray40",
         arrow = arrow(type = "closed", length = unit(0.01, "npc"))) +
    annotate("text", x = 0.54, y = ymax, label = "More likely to avert", color = "white", 
         hjust = 0, size = 3) + 
    annotate("segment", x = 0.48, y = ymax, xend = 0.05, yend = ymax, 
         linejoin = "mitre", linewidth = 5, color = "gray40",
         arrow = arrow(type = "closed", length = unit(0.01, "npc"))) +
    annotate("text", x = 0.46, y = ymax, label = "Less likely to avert", color = "white", 
         hjust = 1, size = 3)

      
    return(time_course_plot)
}

# # # For manuscript, add lines showing potential timing of detection
# fig3 <- plot_time_course() + 
#   geom_vline(xintercept =0.8, linetype = "dashed") + 
#   geom_vline(xintercept = 0.2, linetype = "dotted") + 
#   geom_vline(xintercept = 0, linetype = "dotdash") +
#   theme(legend.background = element_rect(fill = "white"),
#           axis.title.y = element_text(vjust=0), 
#           axis.title =  element_text(face = 'bold')) + 
#   ggtitle("")
  

# As three vertically arranged panels, 
# Plot the distubion of probabilities of detection during cross-section screening 
# (corresponding to the duration of the detectable period),
# and then scatted plots showign the relationship between this duration and the 
# corresponding relative contributions to transmission and mortality. 

# Utility function: Gamma distribution version of the above function: for specified sd and mean,  get the shape and scale parameters.
solve_for_gamma_parameters <- function(mean, sd)
{
  shape <- mean^2 / sd^2
  scale <- sd^2 / mean
  return(list("shape" = shape, "scale" = scale))
}

# Illustrate covariance, assuming lognormal distirubtions with similar coefficient of variation for duration and mortality (and duration and transmission), and specified covariances.
plot_heterogeneity <- function(estimates = midpoint_estimates,
                               N = 500) # just for visualization, too few for stable estimates
{
  
  #  Y = exp(X) where X ~ N(mu, sigma)
  #  E[Y]_1 = 1 = exp(mu_1 + 1/2 sigma_11) ==> mu_1 + 1/2 sigma_11 = 0 ==> mu_1 = -1/2 sigma_11
  #  E[Y]_2 = 1 = exp(mu_2 + 1/2 sigma_22) ==> mu_2 = -1/2 sigma_22
  # cov_12^2 = cov_21^2 = exp(mu_1 + mu_2 + 1/2 sigma_11 + 1/2 sigma_jj) (exp(sigma_12) - 1) = 
  #  1*1*(exp(sigma_12) - 1) = exp(sigma_12) - 1 ==> 
  #  sigma_12 = log(cov_12^2 + 1)
  
  #  and we can choose nearly any sigma_ii's we want here. 
  # Suppose we choose sigma_11 = sigma_22 = 1.
  #  then m_1 = m_2 = -1/2, and sigma_12 = log(cov_12^2 + 1)
  # If sigma_11 = sigma_22 = 1. then correlation = covariance, and max cov is 1. 
  #  So let's set sigma_ii to ceiling(sqrt(covariance)) when cov is positive.

  simulate_correlated_variables <- function(covarianceA, covarianceB, sigma11 = NULL)
  { 
    if (covarianceA < -1/exp(1)) stop("Covariance must be at least -1/e this illustration using log-normals to work, although our DALY model is valid for covariances from -1 to infinity.")

    if(missing(sigma11)) sigma11 <- max(ceiling(sqrt(abs(covarianceA))),
                                        ceiling(sqrt(abs(covarianceB))),
                                        0.2)
    
    sigma22A <- max(ceiling(abs(covarianceA)/sigma11),0.2)
    sigma22B <- max(ceiling(abs(covarianceB)/sigma11),0.2)
    sigma12A <- log(covarianceA^2 + 1)
    sigma12B <- log(covarianceB^2 + 1)
    mu1 <- -1/2 * sigma11
    mu2A <- -1/2 * sigma22A 
    mu2B <- -1/2 * sigma22B
    sigma <- matrix(c(sigma11, sigma12A, sigma12B,  
                      sigma12A, sigma22A, 0,
                      sigma12B, 0, sigma22B), nrow=3)
    mu <- c(mu1, mu2A, mu2B)
    X <- mvrnorm(N, mu, sigma)
    Y <- exp(X)
    colnames(Y) <- c("Y1", "Y2", "Y3")
          
    return(Y)
  }

  # choose a sigma_11 for duration that works for both covariances:
  sigma_duration <- max(ceiling(sqrt(abs(estimates$covariance_mortality_duration))),
                        ceiling(sqrt(abs(estimates$covariance_transmission_duration))),
                        0.2)

  # change any out-of-range covariances to 0 and generate error flag
  covarianceA <- estimates$covariance_mortality_duration
  if (covarianceA < -1/exp(1))
  {
    covarianceA <- -1/exp(1)
    error_flag_mortality <- 1
  } else error_flag_mortality <- 0
  
  covarianceB <- estimates$covariance_transmission_duration
  if (covarianceB < -1/exp(1))
  {
    covarianceB <- -1/exp(1)
    error_flag_transmission <- 1
  } else error_flag_transmission <- 0

  MVN <- simulate_correlated_variables(covarianceA, covarianceB, sigma_duration)

  scatter1 <- 
    ggplot(data=MVN, aes(x=Y1, y=Y2)) + 
    geom_point(alpha=0.3, shape=16)  + 
    geom_rug(col=rgb(.5,0,0,alpha=.2)) + 
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    annotate(geom = "text", x = 1, y = 0.1, 
      label = "Mean duration of detectable period", angle = 90, hjust = 0, vjust = -1) +
    geom_vline(xintercept = 1, linetype = "dashed") + 
    annotate(geom = "text", y = 1, x = 0.9, 
      label = "Average mortality risk per TB episode", angle = 0, hjust = -1, vjust = -1) +
    xlab("Relative duration \n(= relative probability of detection during ACF))") +
    ylab("Relative DALYs\nfrom TB mortality") +
    xlim(0, quantile(MVN[,"Y1"], 0.99)) +
    ylim(0, quantile(MVN[,"Y2"], 0.99))

  if(error_flag_mortality) scatter1 <- scatter1 +
  # add text on top of plot
    annotate(geom = "text", x = 1.1, y = 1.1,
      label = "Covariance too low for\nillustration with lognormals\n(but still valid for DALY model)",
      angle = 0, hjust = 0, vjust = 0, size=6, fontface="bold")
    
  scatter2 <- 
    ggplot(data=MVN, aes(x=Y1, y=Y3)) + 
    geom_point(alpha=0.3, shape=16)  + 
    geom_rug(col=rgb(.5,0,0,alpha=.2)) + 
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    # annotate(geom = "text", x = 1, y = 0.1,
    #   label = "Mean duration of detectable period", angle = 90, hjust = 0, vjust = -1) +
    geom_vline(xintercept = 1, linetype = "dashed") + 
    annotate(geom = "text", y = 1, x = 0.9, 
      label = "Average transmission per TB episode", angle = 0, hjust = -1, vjust = -1) +
    xlab("Relative duration \n(= relative probability of detection during ACF))") +
    ylab("Relative DALYs\nfrom transmission") + 
    annotate(geom = "text", y = 1, x = 0.9, label = "Average transmission per TB episode", 
      angle = 0, hjust = -1, vjust = -1) + 
    xlim(0, quantile(MVN[,"Y1"], 0.99)) +
    ylim(0, quantile(MVN[,"Y3"], 0.99))

  if(error_flag_transmission) scatter2 <- scatter2 + 
  # add text on top of plot
    annotate(geom = "text", x = 1.1, y = 1.1, 
      label = "Covariance too low for\nillustration with lognormals\n(but still valid for DALY model)",
      angle = 0, hjust = 0, vjust = 0, size=6, fontface="bold")

  # Arrange the two figures in a column, with transmission first

  return(grid.arrange(scatter2, scatter1, ncol=1))

}

# plot_heterogeneity()


# Display a table of numerical estimates

output_table <- function(output, forsummary = 0)
{
  if (missing(output)) output <- daly_estimator()

  useoutput <- output %>% 
    pivot_wider(names_from = cumulative_or_averted, values_from = value) %>%
    pivot_wider(names_from = average_or_detected, values_from = c(cumulative, averted)) %>%
    mutate(name = str_to_title(name)) %>%
    bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total")) 

    if (forsummary==1) 
    outputtable <- useoutput %>% 
      select(name, cumulative_average, averted_detected) %>%
      kable(., format = "html",
          digits = 2,
          col.names=c("", "Total cumulative DALYs per case (average case)", "Averted by early detection (detected case)")) %>%
      kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>% 
      row_spec(5, bold = T, hline_after = T) else

    if (forsummary==2) 
    outputtable <- useoutput %>% 
      select(name, cumulative_average, averted_average, averted_detected) %>%
      kable(., format = "html",
          digits = 2,
          col.names=c("", "Total cumulative DALYs per case (average case)", "Averted by early detection (average case)", "Averted by early detection (detected case)")) %>%
      kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>% 
      row_spec(5, bold = T, hline_after = T) else

    if (forsummary==3) 
    outputtable <- useoutput %>% 
      kable(., format = "html",
          digits = 2,
          col.names=c("", "Total cumulative DALYs per case (average case)", "Total cumulative DALYs per case (detected case)", "Averted by early detection (average case)", "Averted by early detection (detected case)")) %>%
      kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>% 
      row_spec(5, bold = T, hline_after = T) else

    outputtable <- useoutput %>% 
    kable(., format = "html",
          col.names=c("Source", rep(c("Average incident case", "Average detected case"), times = 2)),
          digits = 2) %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
    add_header_above(., header =
      c(" " = 1, "Total cumulative DALYs per case" = 2, "Averted by early detection" = 2)) %>% 
    row_spec(5, bold = T, hline_after = T)

return(outputtable)
}

# output_table()

# Now a figure similar to plot_averages() but showing DALYs averted per case detected, using the "Averted by early detection, Average detected case" from table above. 


plot_averted <- function(output, number_labels = TRUE, ymax = NULL)
{
 
  if(missing(output)) output <- daly_estimator()
   
  plot <- ggplot(output %>% filter(cumulative_or_averted == "averted",
                                   average_or_detected == "detected")  %>% 
                   mutate(component = case_when(str_detect(name, "mortality") ~ "TB Mortality",
                                                str_detect(name, "morbidity") ~ "TB Morbidity",
                                                str_detect(name, "sequelae") ~ "Post-TB Sequelae",
                                                str_detect(name, "transmission") ~ "Transmission")) %>%
                   mutate(ordered_component = fct_relevel(component, "Transmission", "Post-TB Sequelae", "TB Mortality", "TB Morbidity")),
                 aes(x=average_or_detected, y=value, fill=ordered_component)) +  
    geom_col(position = "stack", width = 1) +  
    theme_minimal() + xlab("") + ylab("DALYs") + ggtitle("DALYS averted per case detected") + 
    guides(fill="none") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(size = 16))
    
  
  if (number_labels) plot <- plot + geom_text(aes(label = paste0  (ordered_component, ", ",round(value,2))),
                                              position = position_stack(vjust = .5)) else
                                                plot <- plot + geom_text(aes(label = ordered_component),
                                                                         position = position_stack(vjust = .5))

  if (!is.null(ymax)) plot <- plot + ylim(0,ymax)
                                              
 return(plot)
}

# plot_averted()


plot_averted_portion <- function(output, ymax = NULL, base = "detected")
{
 
  if(missing(output)) output <- daly_estimator()
   
  plot <- ggplot(output %>% filter(average_or_detected == base)  %>% 
                   mutate(component = case_when(str_detect(name, "mortality") ~ "TB Mortality",
                                                str_detect(name, "morbidity") ~ "TB Morbidity",
                                                str_detect(name, "sequelae") ~ "Post-TB Sequelae",
                                                str_detect(name, "transmission") ~ "Transmission")) %>%
                  pivot_wider(names_from = "cumulative_or_averted", values_from = "value") %>%
                  mutate(difference = cumulative - averted) %>%
                  pivot_longer(cols = c("cumulative", "averted", "difference"), names_to = "cumulative_or_averted", values_to = "value") %>%
                  filter(cumulative_or_averted != "cumulative") %>%
                  mutate(ordered_component = fct_relevel(component, "Transmission", "Post-TB Sequelae", "TB Mortality", "TB Morbidity")),
                  
                   # make a stacked bar plot of "cumulative", and shade the "averted" portion of each in gray
                aes(x=average_or_detected, y=value, fill=ordered_component, pattern = cumulative_or_averted)) +
          geom_col_pattern(position = "stack", width = 1,
            pattern_fill = 'black',  pattern_spacing = 0.015, pattern_size = 0.02, pattern_density = 0.1) +
          scale_pattern_manual(name = "Averted by early detection?", 
                               values = c("stripe", "none"), 
                               breaks = c("averted", "difference"),
                               labels = c("Averted", "Not averted")) +
         theme_minimal() + xlab("") + ylab("DALYs") + 
        guides(fill = guide_legend(title = "DALY component")) +
        theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(size = 16),
          legend.position = "bottom",
          # stack the legends vertically for fill and pattern
          legend.box = "vertical")

    
  if (base == "detected") plot <- plot + ggtitle("DALYs per *detected* case") else
    plot <- plot + ggtitle("DALYS per *average* case")
  
  if (!is.null(ymax)) plot <- plot + ylim(0,ymax)
                                              
 return(plot)
}

# plot_averted_portion()

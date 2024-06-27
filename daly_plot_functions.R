
# setwd("~/Google Drive/My Drive/DALY impact of ACF 2024/ACF-impact-Rshiny/")
source("daly_estimator.R")

library(gridExtra)
library(kableExtra)

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

  plotdata <- rbind(within_case,
                    within_case %>% filter(cumulative_or_averted=="cumulative") %>%
                      mutate(cumulative_or_averted="detectable",
                             value = value*case_when(name=="transmission" ~ (1-estimates$predetection_transmission -    estimates$postrx_transmission),
                                                     TRUE ~ (1-estimates$predetection_mm - estimates$postrx_mm))),
                    within_case %>% filter(cumulative_or_averted=="cumulative") %>%
                      mutate(cumulative_or_averted="pre",
                             value = value*case_when(name=="transmission" ~ estimates$predetection_transmission,
                                                     TRUE ~ estimates$predetection_mm)),
                    within_case %>% filter(cumulative_or_averted=="cumulative") %>%
                      mutate(cumulative_or_averted="post",
                             value = value*case_when(name=="transmission" ~ estimates$postrx_transmission,
                                                     TRUE ~ estimates$postrx_mm)))
                   
  
  # If second half area is s fraction of total, with arbitrary detectable period width 2 (from -1 to +1) and midpoint height 1, 
  # then overall detectable area is 2;second half area is 2s;
  # ending height h2 is such that (h2+1)/2*1= 2s -> h2 = 4s - 1;
  # and initial height h1 of the first half is such that (h1+1)/2*1= 2-2s) --> h1 = 3 - 4s 
  
  h1_transmission <- 3 - 4 * estimates$second_half_vs_first_transmission
  h2_transmission <- 4 * estimates$second_half_vs_first_transmission - 1
  h1_mm <- 3 - 4 * estimates$second_half_vs_first_mm
  h2_mm <- 4 * estimates$second_half_vs_first_mm - 1
  
  # we're going to go around the polygon, first across the top. left to right, then bottom right to left. 
  # xs are defined by start of predetect for transmission and mm, then start/end of detectable, then end of post rx for transmission and mm. 
  # we'll define heights at each point (with bases at 0), then stack. 
  rect_points <- as_tibble(t(c("x1" = -1,
                               "x2" = -1,
                               "x3" = 1,
                               "x4" = 1)))
  
  rect_points <- rbind(rect_points %>% mutate(component = "morbidity"),
                       rect_points %>% mutate(component = "mortality"),
                       rect_points %>% mutate(component = "sequelae"),
                       rect_points %>% mutate(component = "transmission"))
                       
  
  rect_points <- rect_points %>% mutate(y1=0,
                                        y2=unlist(c(h1_mm, h1_mm, h1_mm, h1_transmission)*
                                                    (within_case %>% filter(cumulative_or_averted=="cumulative") %>% 
                                                     arrange(match(name, rect_points$component)) %>% select(value))),
                                        y3= unlist(c(h2_mm, h2_mm, h2_mm, h2_transmission)*
                                                     (within_case %>% filter(cumulative_or_averted=="cumulative") %>% 
                                                        arrange(match(name, rect_points$component)) %>% select(value))),
                                        y4=0)
  # now need to stack them
  rect_points_stacked <- rect_points %>% mutate(y2 = cumsum(y2), y3 = cumsum(y3)) %>%
                                          mutate(y1 = c(0, y2[1:3]), y4 = c(0, y3[1:3]))

  toplot <- rect_points_stacked %>%
    pivot_longer(-component,
                 names_to = c(".value", "id"),
                 names_pattern = "(\\D)(\\d+)") %>%
    mutate(component = factor(component, levels = rev(c("morbidity","mortality", "sequelae", "transmission"))))

  time_course_plot <- ggplot(data = toplot) +
    geom_polygon(aes(x = x, y = y, group = component, fill = component)) +
    scale_fill_discrete(breaks = rev(levels(toplot$component))) +
    xlab("Time (arbitrary scale)") +
    ylab("DALY accrual rate (arbitrary scale)") +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    scale_y_discrete(labels = NULL, breaks = NULL) +
    theme_minimal() +
    geom_rect(data = rect_points_stacked, aes(xmin = min(x1), xmax = 2*min(x1), ymin = 0, ymax = max(y3)),
              fill = "gray", alpha = 0.3) +
    geom_vline(data = rect_points_stacked, aes(xintercept = min(x4))) +
    annotate(geom = "text", x = 1.5*max(rect_points_stacked$x1), y = max(rect_points_stacked$y4) / 2,
             label = "before detectability", angle = 90) +
    annotate(geom = "text", x = 1.1*min(rect_points_stacked$x4), y = max(rect_points_stacked$y4) / 2,
             label = "after routine detection", angle = 90) +
    guides(fill = guide_legend(reverse = TRUE)) +
    ggtitle("Timing of DALY accrual") +
    theme(axis.text = element_text(size = 14),
          plot.title = element_text(size = 16))


      
    return(time_course_plot)
}



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

plot_heterogeneity <- function(estimates = midpoint_estimates,
                               N = 500) # just for visualization, too few for stable estimates
{
  qs <- runif(N, 0, 1)
  qs_mort <- qs_trans <- qs
  if (estimates$duration_tbdeath_covarying_cv < 0) qs_mort <- 1 - qs
  if (estimates$duration_transmission_covarying_cv < 0) qs_trans <- 1 - qs

  relative_durations <-  qgamma(p = qs,
                                shape = (solve_for_gamma_parameters(mean = 1, sd = estimates$duration_cv))[['shape']],
                                scale = (solve_for_gamma_parameters(mean = 1, sd = estimates$duration_cv))[['scale']])
    
  relative_dalys_mortality <-
    qgamma(p = qs_mort,
           shape=(solve_for_gamma_parameters(
                    mean = 1,
                    sd = estimates$duration_cv * abs(estimates$duration_tbdeath_covarying_cv)))[['shape']],
           scale=(solve_for_gamma_parameters(
                    mean = 1,
                    sd = estimates$duration_cv * abs(estimates$duration_tbdeath_covarying_cv)))[['scale']])
  relative_dalys_mortality_withnoise <-
      rnorm(N, mean = relative_dalys_mortality, sd = 0.2 * relative_dalys_mortality)
  
  relative_dalys_transmission <- qgamma(p = qs_trans,
                                        shape=(solve_for_gamma_parameters(mean=1, sd=estimates$duration_cv*abs(estimates$duration_transmission_covarying_cv)))[['shape']],
                                        scale=(solve_for_gamma_parameters(mean=1, sd=estimates$duration_cv*abs
                                        (estimates$duration_transmission_covarying_cv)))[['scale']]) 
  relative_dalys_transmission_withnoise <- 
      rnorm(N, mean = relative_dalys_transmission, sd = 0.2 * relative_dalys_transmission)
    
# Start first of three ggplot figures to be arranged in a column of panels
figure1 <- ggplot(as_tibble(relative_durations)) + 
  geom_density(aes(x=value), fill="blue", alpha=0.5) +
  xlim(0, quantile(relative_durations, 0.999)) +
  theme_minimal() + 
  labs(x = NULL, y = "Proportion of incident TB") +
    theme(axis.text.y = element_blank(), plot.margin = margin(10, 10, 10, 10),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 14))
  
# Now scatter plot of relative durations vs relative DALYs for mortality
figure2 <- ggplot(as_tibble(cbind(relative_durations, relative_dalys_mortality_withnoise))) + 
  geom_point(aes(x=relative_durations, y=relative_dalys_mortality_withnoise), alpha=0.5) +
  xlim(0, quantile(relative_durations, 0.999)) +
  theme_minimal() + 
  labs(x = NULL, y = "Relative DALYs\nfrom TB mortality") +
  theme(plot.margin = margin(10, 10, 10, 10),
        axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 14))

# and vs relative DALYs for transmission
figure3 <- ggplot(as_tibble(cbind(relative_durations, relative_dalys_transmission_withnoise))) + 
  geom_point(aes(x=relative_durations, y=relative_dalys_transmission_withnoise), alpha=0.5) +
  xlim(0, quantile(relative_durations, 0.999)) +
  theme_minimal() +
  labs(x = "Relative duration \n(= relative probability of detection during ACF))", y = "Relative DALYs\nfrom transmission") +
  theme(plot.margin = margin(10, 10, 10, 10),
        # increase axis label size
        axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 14))

# Arrange the three figures in a column

return(grid.arrange(figure1, figure2, figure3, ncol=1))

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



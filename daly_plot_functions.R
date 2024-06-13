source("daly_estimator.R")

plot_averages <- function(output, number_labels = TRUE)
{
 
  if(missing(output)) output <- dalys_per_average_case()
   
  plot <- ggplot(output %>% filter(cumulative_or_averted == "cumulative",
                                   average_or_detected == "average")  %>% 
                   mutate(component = case_when(str_detect(name, "mortality") ~ "TB Mortality",
                                                str_detect(name, "morbidity") ~ "TB Morbidity",
                                                str_detect(name, "sequelae") ~ "Post-TB Sequelae",
                                                str_detect(name, "transmission") ~ "Transmission")) %>%
                   mutate(ordered_component = fct_relevel(component, "Transmission", "Post-TB Sequelae", "TB Mortality", "TB Morbidity")),
                 aes(x=average_or_detected, y=value, fill=ordered_component)) +  
    geom_col(position = "stack", width = 1) +  
    theme_minimal() + xlab("") + ylab("DALYs") + ggtitle("Total DALYS per average case") + 
    guides(fill="none") 
  
  if (number_labels) plot <- plot + geom_text(aes(label = paste0(ordered_component, ", ",round(value,1))),
                                              position = position_stack(vjust = .5)) else
                                                plot <- plot + geom_text(aes(label = ordered_component),
                                                                         position = position_stack(vjust = .5))
                                              
                                              return(plot)
}

plot_averages()

plot_detectable_proportion <- function(totals_plot, estimates=midpoint_estimates)
{
  if(missing(totals_plot)) totals_plot <- plot_averages(dalys_per_average_case(estimates), number_labels = F)
  
  with(estimates, {
    detectable_period_plot <- 
      totals_plot + geom_rect(aes(xmin=0.5, xmax=0.5 + predetection_mm, ymin=0, 
                                  ymax= sum(totals_plot$data %>% filter(name!="transmission") %>% select(value))), 
                              fill="gray", alpha=0.4) + 
      geom_rect(aes(xmin=0.5, xmax= 0.5 + predetection_transmission,
                    ymin= sum(totals_plot$data %>% filter(name!="transmission") %>% select(value))), 
                ymax= sum(totals_plot$data %>% select(value)), 
                fill="gray", alpha=0.4) + 
      annotate(geom="text", x=0.5 + predetection_mm/2, y=0.1, 
               label="Accrues before detectability", angle=90, hjust=0) +
      geom_rect(aes(xmin=1.5 - postrx_mm, xmax=1.5, ymin=0, 
                    ymax= sum(totals_plot$data %>% filter(name!="transmission") %>% select(value))), 
                fill="gray", alpha=0.4) + 
      geom_rect(aes(xmin=1.5 - postrx_transmission, xmax= 1.5,
                    ymin= sum(totals_plot$data %>% filter(name!="transmission") %>% select(value))), 
                ymax= sum(totals_plot$data %>% select(value)), 
                fill="gray", alpha=0.4) + 
      annotate(geom="text", x=1.5 - postrx_mm/2, y=0.1, 
               label="Accrues after routine diagnosis", angle=90, hjust=0) + 
      ggtitle("Avertible DALYs per average case")
    return(detectable_period_plot)
  } )
}

plot_detectable_proportion()

plot_time_course <- function(estimates=midpoint_estimates)
{
  proportions <- within_case_avertible(estimates) 
  averages_and_avertibles <- within_case_cumulative_and_averted(estimates)

  plotdata <- rbind(averages_and_avertibles,
                    averages_and_avertibles %>% filter(cumulative_or_averted=="cumulative") %>%
                      mutate(cumulative_or_averted="detectable",
                             value = value*case_when(name=="transmission" ~ (1-estimates$predetection_transmission - estimates$postrx_transmission),
                                                     TRUE ~ (1-estimates$predetection_mm - estimates$postrx_mm)),
                    averages_and_avertibles %>% filter(cumulative_or_averted=="cumulative")) %>%
                      mutate(cumulative_or_averted="pre",
                             value = value*case_when(name=="transmission" ~ estimates$predetection_transmission,
                                                     TRUE ~ estimates$predetection_mm)),
                    averages_and_avertibles %>% filter(cumulative_or_averted=="cumulative") %>%
                      mutate(cumulative_or_averted="post",
                             value = value*case_when(name=="transmission" ~ estimates$postrx_transmission,
                                                     TRUE ~ estimates$postrx_mm)))
                   
  
  # If second half area is s fraction of total, with arbitrary detectable period width 2 (from -1 to +1) and midpoint height 1, 
  # then overall detectable area is 2;second half area is 2s;
  # ending height h2 is such that (h2+1)/2*1= 2s -> h2 = 4s - 1;
  # and initial height h1 of the first half is such that (h1+1)/2*1= 2-2s) --> h1 = 3 - 4s 
  
  # If pre period decays linearly to zero, and ratio of detectable to pre area is r,
  # then width of pre period w should be selected such that w*h1/2 = 2*r --> w =4r/h1 = 4r/(3-4s)
  # Similarly, if ratio of dertectable to post area is r2, then w2*h2/2 = 2*r2 --> w2 =4 r2/h2 = 4r2/(4s-1)
  
  # So now we can plot this.
  # Will need to also adjust for pre and post paramertrs beig fractions of total, not of detectable, 
  #  i.e. r = pre/detectable = pre/(1-pre-post) 

  r1_transmission <- estimates$predetection_transmission/(1-estimates$predetection_transmission-estimates$postrx_transmission)
  r2_transmission <- estimates$postrx_transmission/(1-estimates$predetection_transmission-estimates$postrx_transmission)
  r1_mm <- estimates$predetection_mm/(1-estimates$predetection_mm-estimates$postrx_mm)
  r2_mm <- estimates$postrx_mm/(1-estimates$predetection_mm-estimates$postrx_mm)
  
  h1_transmission <- 3 - 4*estimates$second_half_vs_first_transmission
  h2_transmission <- 4*estimates$second_half_vs_first_transmission - 1
  h1_mm <- 3 - 4*estimates$second_half_vs_first_mm 
  h2_mm <- 4*estimates$second_half_vs_first_mm - 1
  
  # we're going to go around the polygon, first across the top. left to right, then bottom right to left. 
  # xs are defined by start of predetect for transmission and mm, then start/end of detectable, then end of post rx for transmission and mm. 
  # we'll define heights at each point (with bases at 0), then stack. 
  transmission_predetection_wider <- estimates$predetection_transmission > estimates$predetection_mm
  transmission_postrx_wider <- estimates$postrx_transmission > estimates$postrx_mm
  rect_points <- as_tibble(t(c("x1"=ifelse(transmission_predetection_wider, -1-(4*r1_transmission/h1_transmission), -1-(4*r1_mm/h1_mm)), 
                               "x2"=ifelse(transmission_predetection_wider, -1-(4*r1_mm/h1_mm), -1-(4*r1_transmission/h1_transmission)),
                               "x3"=-1,
                               "x4"=1,
                               "x5"=ifelse(transmission_postrx_wider, 1+(4*r2_mm/h2_mm), 1+(4*r2_transmission/h2_transmission)),
                               "x6"=ifelse(transmission_postrx_wider, 1+(4*r2_transmission/h2_transmission), 1+(4*r2_mm/h2_mm))))) %>% 
    mutate(x7=x5, x8=x4, x9=x3, x10=x2)
  
  rect_points <- rbind(rect_points %>% mutate(component = "morbidity"),
                       rect_points %>% mutate(component = "mortality"),
                       rect_points %>% mutate(component = "sequelae"),
                       rect_points %>% mutate(component = "transmission"))
                       
  
  rect_points <- rect_points %>% mutate(y1=0,
                                        y3=unlist(c(h1_transmission, h1_mm, h1_mm, h1_mm)*
                                                    (averages_and_avertibles %>% filter(cumulative_or_averted=="cumulative") %>% 
                                                     arrange(match(name, rect_points$component)) %>% select(value))),
                                        y4= unlist(c(h2_transmission, h2_mm, h2_mm, h2_mm)*
                                                     (averages_and_avertibles %>% filter(cumulative_or_averted=="cumulative") %>% 
                                                        arrange(match(name, rect_points$component)) %>% select(value))),
                                        y6=0,
                                        y7=0, y8=0, y9=0, y10=0) %>%
                                mutate(
                                        y2=y3*case_when(component=="transmission" ~ 
                                                       ifelse(transmission_predetection_wider, (x2+1)/(x1+1), 0),
                                                     component!="transmission" ~ 
                                                       ifelse(transmission_predetection_wider, 0, (x2+1)/(x1+1))),
                                        y5=y4*case_when(component=="transmission" ~ 
                                                       ifelse(transmission_postrx_wider, (x5-1)/(x6-1), 0),
                                                     component!="transmission" ~ 
                                                       ifelse(transmission_postrx_wider, 0, (x5-1)/(x6-1))))


  # now need to stack them
  rect_points_stacked <- rect_points %>% mutate(y2 = cumsum(y2), y3 = cumsum(y3), y4 = cumsum(y4), y5 = cumsum(y5)) %>%
                                          mutate(y7 = c(0, y5[1:3]), y8 = c(0, y4[1:3]), y9 = c(0, y3[1:3]), y10 = c(0, y2[1:3]))

    
  toplot <- rect_points_stacked %>% pivot_longer(-component, 
               names_to = c(".value","id"), 
               names_pattern = "(\\D)(\\d+)" )
  

  ggplot(data = toplot) + 
    geom_polygon(aes(x=x, y=y, group=component, fill=component)) + 
  xlab("Time (arbitrary scale)") + ylab ("DALY accrual rate (arbitraty scale)") + 
    scale_x_discrete(labels = NULL, breaks = NULL) + scale_y_discrete(labels = NULL, breaks = NULL) + 
    theme_minimal() + 
    geom_rect(data=rect_points_stacked, aes(xmin=min(x1), xmax=max(x3), ymin=0, ymax=max(y4)), col="gray", alpha=0.1) + 
    geom_rect(data=rect_points_stacked, aes(xmin=min(x4), xmax=max(x6), ymin=0, ymax=max(y4)), col="gray", alpha=0.1) +
    geom_vline(xintercept = -1, linetype=2) + geom_vline(xintercept = 1, linetype=2) + 
    annotate(geom = "text", x= (min(rect_points_stacked$x1)-1)/2, y=max(rect_points_stacked$y4)/2, label="before detectability", angle=90) + 
    annotate(geom = "text", x= (max(rect_points_stacked$x6)+1)/2, y=max(rect_points_stacked$y4)/2, label="after routine detection", angle=90)

    
    # reshape the above so that it's on a time scale. 
    # Mark start and end of detectable periods. 
    # - Same width for mm and transmission (and adjust the heights).
    # - Plot increase within detectable period, considering that mean of first half coccurs at 1/4 t, and of second half at 3/4 t. 
      # --     0.5*resolving_detectable + second_half_vs_first_mm*(1-resolving_detectable)
    # - Plot the before-detectability and after-diagnosis areas as decays to arbitrary time points, each ~1/3 of detectable period?. 
    
    

    
    
    # from prior funciton:
        detectable_period_plot <- 
      totals_plot + geom_rect(aes(xmin=0.5, xmax=0.5 + predetection_mm, ymin=0, 
                                  ymax= sum(totals_plot$data %>% filter(name!="transmission") %>% select(value))), 
                              fill="gray", alpha=0.4) + 
      geom_rect(aes(xmin=0.5, xmax= 0.5 + predetection_transmission,
                    ymin= sum(totals_plot$data %>% filter(name!="transmission") %>% select(value))), 
                ymax= sum(totals_plot$data %>% select(value)), 
                fill="gray", alpha=0.4) + 
      annotate(geom="text", x=0.5 + predetection_mm/2, y=0.1, 
               label="Accrues before detectability", angle=90, hjust=0) +
      geom_rect(aes(xmin=1.5 - postrx_mm, xmax=1.5, ymin=0, 
                    ymax= sum(totals_plot$data %>% filter(name!="transmission") %>% select(value))), 
                fill="gray", alpha=0.4) + 
      geom_rect(aes(xmin=1.5 - postrx_transmission, xmax= 1.5,
                    ymin= sum(totals_plot$data %>% filter(name!="transmission") %>% select(value))), 
                ymax= sum(totals_plot$data %>% select(value)), 
                fill="gray", alpha=0.4) + 
      annotate(geom="text", x=1.5 - postrx_mm/2, y=0.1, 
               label="Accrues after routine diagnosis", angle=90, hjust=0) + 
      ggtitle("Avertible DALYs per average case")
    return(detectable_period_plot)
  } )
}

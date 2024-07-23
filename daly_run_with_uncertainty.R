library(gtools) # for logit functions
library(tidyverse)
library(sensitivity) #prcc
library(clipr)
library(lhs)
library(abind)
library(kableExtra)
library(gridExtra)

source("daly_estimator.R")
daly_estimator()


paramdf <- read.csv("DALY_model_param_values.csv", header = T, row.names = 1)
midpoint_estimates <- as.list(paramdf$mid); names(midpoint_estimates) <- rownames(paramdf)

  
#### Uncertainty and sensitivity analyses ####

# Will sample params from transformed beta distributions with desired mean (mid) and support (low to high)

# First need to estimate the correct beta parameters for each model parameter: 
get_beta_pars <- function(paramdf, 
                          paramname,
                          precision = 2)
{
  mid = paramdf[paramname, "mid"]
  low = paramdf[paramname, "low"]
  high = paramdf[paramname, "high"]
  beta_mu <- (mid-low)/(high-low) # also = alpha/(alpha+beta)
  beta_alphapar <- beta_mu * precision
  beta_betapar <- precision - beta_alphapar
  return(list("alpha"=beta_alphapar, "beta"=beta_betapar))
}

for (paramname in paramdf$param)
{
  pars <- get_beta_pars(paramdf, paramname)
  paramdf$beta_alphapar[paramdf$param==paramname] <- pars$alpha
  paramdf$beta_betapar[paramdf$param==paramname] <- pars$beta
}

# will transform these after sampling to get actual param vals

# check low and high values of all
e <- as.list(paramdf$low); names(e) <- rownames(paramdf)
daly_estimator(estimates = e)  
e <- as.list(paramdf$high); names(e) <- rownames(paramdf)
daly_estimator(estimates = e)  

# sample all parameters independently from beta uncertainty distributions, using latin hypcercube sampling
N_paramsets <- 1000
paramsets <- array(dim = c(N_paramsets, nrow(paramdf)), dimnames = list("N"=1:N_paramsets, "paramname"=paramdf$param))
# sample a latin hypercube with dimensions N_paramsets x nrow(paramdf)
LHS <- randomLHS(N_paramsets, nrow(paramdf))
dimnames(LHS) <- list("N" = 1:N_paramsets, "paramname"=paramdf$param)

for (paramname in paramdf$param)
{
  paramsets[, paramname] <-
    paramdf[paramname, "low"] +
    (paramdf[paramname, "high"] - paramdf[paramname, "low"]) *
  
      qbeta(p = LHS[, paramname],
            shape1 = paramdf[paramname, "beta_alphapar"],
            shape2 = paramdf[paramname, "beta_betapar"])
}

# Generate results for all sampled parameter sets
sampled_dalys <- apply(paramsets, 1, function(x) daly_estimator(estimates = as.list(x)))
# add list name as column named "paramset" to each element of sampled_dalys
sampled_dalys_withparamsetnumbers <- 
  lapply(seq_along(sampled_dalys), 
         function(i) sampled_dalys[[i]] %>% 
          mutate(paramset = names(sampled_dalys)[[i]]))

sampled_dalys_df <- do.call(rbind, sampled_dalys_withparamsetnumbers)

# Show table of results
  fortable <- sampled_dalys_df %>% 
    pivot_wider(names_from = cumulative_or_averted, values_from = value) %>%
    pivot_wider(names_from = average_or_detected, values_from = c(cumulative, averted)) %>%
    mutate(name = str_to_title(name)) %>%
    # summarize total by paramset group 
    group_by(paramset) %>%
    bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))

# For use in tables etc:
summarize_uncertainty <- function(outcome_vector)
{return(sprintf(
  "%0.2f [%0.2f - %0.2f]",
  median(outcome_vector), quantile(outcome_vector, 0.025), quantile(outcome_vector, 0.975)))}


  # For each "name" and non-paramset column, calculate median and 95% inner quantile range
 fortable %>%
    group_by(name) %>%
    select(-paramset) %>%
    summarise_all(~if(is.numeric(.)) summarize_uncertainty(.)
                  else name) %>%
    mutate(name = str_to_title(name)) %>% 
    # move "Total" row to last
    arrange(name == "Total") %>%
    kable(., format = "html",
          col.names=c("Source", rep(c("Average incident case", "Average detected case"), times = 2)),
          digits = 2) %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
    add_header_above(., header =
      c(" " = 1, "Total cumulative DALYs per case" = 2, "Averted by early detection" = 2)) %>% 
    row_spec(5, bold = T, hline_after = T)



# Boxplot of results
result_boxplot <- ggplot(data = fortable %>% mutate(name = factor(name, 
                              levels = c("Morbidity", "Mortality", "Sequelae", "Transmission", "Total"))), 
       aes(x = name, y = averted_detected)) +
  geom_boxplot() +
  ylab("DALYs averted per case detected") + xlab("DALY Component") +
  coord_flip() +
  theme_minimal()



# # plot all the components of averted DALYs
# plotall_data <- as.data.frame(t(sampled_dalys_alloutcomes)) %>% 
#   mutate(
#     morbidity_detected_case = morbidity_average_case,
#     mortality_detected_case = mortality_average_case * avertible_mortality_multiplier_detected,
#     sequelae_detected_case = sequelae_average_case,
#     transmission_detected_case = transmission_average_case * avertible_transmission_multiplier_detected,
#     morbidity_avertible_average_case = morbidity_average_case * average_proportion_avertible_mm,
#     mortality_avertible_average_case = mortality_average_case * average_proportion_avertible_mm,
#     sequelae_avertible_average_case = sequelae_average_case * average_proportion_avertible_mm,
#     transmission_avertible_average_case = transmission_average_case * average_proportion_avertible_transmission
#   ) %>% 
#   mutate(
#     total_detected_case = morbidity_detected_case + mortality_detected_case + sequelae_detected_case + transmission_detected_case,
#     total_avertible_average_case = morbidity_avertible_average_case + mortality_avertible_average_case + sequelae_avertible_average_case + transmission_avertible_average_case
#   ) %>%
#   select(
#     total_average_case, transmission_average_case, morbidity_average_case, mortality_average_case, sequelae_average_case,
#     total_detected_case, transmission_detected_case, morbidity_detected_case, mortality_detected_case, sequelae_detected_case,
#     total_avertible_average_case, transmission_avertible_average_case, morbidity_avertible_average_case, mortality_avertible_average_case, sequelae_avertible_average_case,
#     total_averted, transmission_averted, morbidity_averted, mortality_averted, sequelae_averted) %>%
#   pivot_longer(cols = everything()) %>% 
#   mutate(total_or_component = case_when(str_starts(name, "total_") ~ "total", TRUE ~ "component"),
#          average_or_detected = case_when(str_ends(name, "_average_case") ~ "average", TRUE ~ "detected"),
#          cumulative_or_averted = case_when(str_detect(name, "_avert") ~ "avertible", TRUE ~ "cumulative"))

# # show 2x2 grid of total and averted, for average case and detected case?
# ggplot(plotall_data) + 
#   geom_boxplot(aes(x=name, y=value, fill = total_or_component)) + 
#   facet_wrap(average_or_detected ~ cumulative_or_averted)

# # separate total from components, and show components as stacked bars. don't wrap.
# ggplot(plotall_data %>% filter(total_or_component=="total") %>%
#          mutate(nicename = fct_recode(as.factor(name),
#                                   "Average case,\ncumulative" = "total_average_case", 
#                                   "Average case,\navertible" = "total_avertible_average_case",
#                                   "Detected case,\ncumulative" = "total_detected_case", 
#                                   "Detected case,\navertible" = "total_averted")) %>% 
#          mutate(orderedname = fct_relevel(nicename, "Average case,\ncumulative", 
#                                           "Average case,\navertible", 
#                                          "Detected case,\ncumulative", 
#                                           "Detected case,\navertible"))) + 
#   geom_boxplot(aes(x=orderedname, y=value)) + 
#   ylab("DALYs") + theme_minimal()

# ggplot(plotall_data %>% filter(total_or_component!="total") %>%
#          mutate(component = case_when(str_detect(name, "mortality") ~ "TB Mortality",
#                                       str_detect(name, "morbidity") ~ "TB Morbidity",
#                                       str_detect(name, "sequelae") ~ "Post-TB Sequelae",
#                                       str_detect(name, "transmission") ~ "Transmission")) %>%
#          group_by(average_or_detected, cumulative_or_averted, component) %>% summarise(dalys=median(value)) %>%
#          mutate(label = as.factor(str_to_sentence(paste0(average_or_detected, " case,\n",cumulative_or_averted)))) %>%
#          mutate(orderedlabel = fct_relevel(label, 
#                                            "Average case,\nCumulative", 
#                                             "Average case,\nAvertible", 
#                                             "Detected case,\nCumulative", 
#                                             "Detected case,\nAvertible"),
#                 orderedcomponent = fct_relevel(component, "Transmission", "Post-TB Sequelae", "TB Mortality", "TB Morbidity"))) + 
#   geom_col(aes(x=orderedlabel, y=dalys, fill=orderedcomponent)) + 
#   theme_minimal() + xlab("Case type and timeframe considered") + ylab("DALY estimate") + 
#   guides(fill=guide_legend(title="Component"))

# merge paramsets with nicenames
params_with_names <- data.frame(paramsets)
colnames(params_with_names) <- unlist(nicenames[colnames(paramsets)])

prccs <- pcc(X = data.frame(params_with_names),
             y = fortable %>%
                  filter(name=="Total") %>%
                  arrange(as.numeric(paramset)) %>%
                  select(averted_detected) %>%
                  pull())
plot(prccs)

prcc_plot <- ggplot(cbind(prccs$PCC, name=rownames(prccs$PCC)), 
    aes(x = fct_reorder(name, abs(original)), y=original)) + geom_col() + 
  coord_flip() + ylab("Partial rank correlation coefficient (PRCC) with DALYs averted per case detected") + xlab("Parameter") + theme_minimal() +
  # remove "."s from x axis labels
  scale_x_discrete(labels = function(x) gsub("\\.", " ", x))

# arrange result_boxplot and prcc_plot in a two-panel figure
# and move 
grid.arrange(result_boxplot + 
              ggtitle("A") + 
              theme(plot.title.position = "plot"), 
            prcc_plot + 
              ggtitle("B") +
              theme(plot.title.position = "plot"), 
            ncol=1) 



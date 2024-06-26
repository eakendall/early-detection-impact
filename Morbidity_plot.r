# Make an illustrative figure of the morbidity of TB over time
# X axis is "Time with TB (arbitrary scale)"
# Y axis is "Symptom severity (arbitrary scale)"
# There are horizontal lines indicating three different symptom thresholds: "Mild, nonspecific symptoms", "Recognized symptoms", and "Severe symptoms" 
# There is a monotonically increasing curve that starts at zero and ends at 1, representing the progression of symptoms over time. We'll make up data to define this. 
# The total area under this curve is filled in a color representing the total time with TB,
# and the areas under this curve but above each of the symptom thresholds are hatched in different directions to represent the time spent with each level of symptoms

# Make up data for how sypmtoms progress over time
time <- seq(0,1,0.01)
symptoms <- time^2
# combine into a dataframe
symptom_data <- data.frame(time, symptoms)

# Make up data for the symptom thresholds
thresholds <- c(0.2, 0.45, 0.7)

# Make the plot using ggplot2
library(ggplot2)
library(ggpattern)
ggplot(symptom_data) + geom_area(aes(x = time, y = symptoms), fill = "cadetblue1") +
  geom_hline(yintercept = thresholds, linetype = "dashed") +
    # fill in the area that's below the curve but above the lowest threshold, in red
    geom_ribbon(data = symptom_data %>% filter(symptoms>thresholds[1]), aes(x = time, ymax = symptoms, ymin=thresholds[1]), fill = "cadetblue2") +
   geom_ribbon(data = symptom_data %>% filter(symptoms>thresholds[2]), aes(x = time, ymax = symptoms, ymin=thresholds[2]), fill = "cadetblue3") +
      geom_ribbon(data = symptom_data %>% filter(symptoms>thresholds[3]), aes(x = time, ymax = symptoms, ymin=thresholds[3]), fill = "cadetblue4") +
  # label the thresholds using annotate
    annotate("text", x = 0.1, y = thresholds[1] + 0.01, label = "Mild, nonspecific symptoms", hjust = 0, vjust = 0) +
    annotate("text", x = 0.1, y = thresholds[2] + 0.01, label = "Recognized symptoms", hjust = 0, vjust = 0) +
    annotate("text", x = 0.1, y = thresholds[3] + 0.01, label = "Severe symptoms", hjust = 0, vjust = 0) +
  theme_minimal() +
  labs(x = "Time (arbitrary scale)", y = "Symptom severity (arbitrary scale)", title = "Illustrative TB disease course") + 
    theme(axis.text = element_blank())

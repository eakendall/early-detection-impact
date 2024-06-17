library(shiny)
library(bslib)

source("daly_estimator.R")
source("daly_plot_functions.R")

paramdf <- read.csv("DALY_model_param_values.csv", header = T)

sliderInput_fromfile <- function(id, label, paramtable=paramdf, step=NULL) {
  sliderInput(id, label, 
              min = round(paramtable[paramtable$param==id,"low"],2), 
              max = round(paramtable[paramtable$param==id,"high"],2), 
              value = round(paramtable[paramtable$param==id,"mid"],2),
              step=step)
}

ui <- fluidPage(
  
  # App title ----
  titlePanel("DALY impact of early case detection"),
  
  navlistPanel(
    nav_panel(
      bslib::card(
      title = "Total DALYs associated with an average TB case",

      shiny::tags$h5("TB morbidity"),
      sliderInput_fromfile("tb_symptom_duration", "Duration of symptomatic TB (years)"),
      sliderInput_fromfile("tb_symptom_dw", "Disability weight while symptomatic with TB"),

      shiny::tags$h5("TB mortality"),
      sliderInput_fromfile("tb_cfr", "TB case fatality ratio"),
      sliderInput_fromfile("tb_death_yearslost", "Years of life lost per TB death"),
      
      shiny::tags$h5("Sequelae"),
      sliderInput_fromfile("posttb_symptom_duration", "Time lived post-TB (years)"),
      sliderInput_fromfile("posttb_symptom_dw", "Disability weight after TB (average, including those without sequelae)", step=0.02),
      sliderInput_fromfile("posttb_cfr", "Post-TB fatality ratio (proportion who die early of TB sequelae)", step=0.01),
      sliderInput_fromfile("posttb_death_yearslost", "Years of life lost per death from TB sequelae"),
      sliderInput_fromfile("posttb_death_timing", "Mean time to death from TB sequelae (years)"),
      
      shiny::tags$h5("Temporal discounting"),
      sliderInput_fromfile("discounting_rate", "Annual discounting rate on health outcomes", step=0.005),
    
      shiny::tags$h5("Transmission"),
      sliderInput_fromfile("downstream_cases", "# of downstream cases attributable to average case [already discounted]", step=0.25),  #! consider changing this and adding discounting to this model instead
      
      tableOutput("testtable"),
      plotOutput("averages_plot") 
      )
    )

  #   bslib::card(
  #     title = "Timing of DALY accrual",
      
      
  #     sliderInput_fromfile("predetection_mm", "Proportion of morbidity and mortality that accrue before TB becomes detectable by screening algorithm"),
  #     sliderInput_fromfile("predetection_transmission", "Proportion of transmission that occurs before TB becomes detectable by screening algorithm"),
  #     sliderInput_fromfile("postrx_mm", "Proportion of morbidity and mortality that accrue after routine detection"),
  #     sliderInput_fromfile("postrx_transmission", "Proportion of transmission that occurs after routine detection"),

  #     sliderInput_fromfile("second_half_vs_first_mm", "Of morbidity, mortality, and post-TB disease that accrue during detectable period, proportion in 2nd half"),
  #     sliderInput_fromfile("second_half_vs_first_transmission", "Of transmission that occurs during detectable period, proportion in 2nd half"),
  #     sliderInput_fromfile("resolving_detectable", "Proportion of TB that will spontaneously resolve after becoming detectable"),

  #     plotOutput("time_course_plot")
  #   ),

  #   bslib::card(
  #     title = "How TB found through screening differs from the average case",

  #       sliderInput_fromfile("duration_cv", "Coefficient of variation in TB duration"),
  #       sliderInput_fromfile("duration_tbdeath_covarying_cv", "How mortality risk co-varies with total duration (in absence of ACF)"),
  #       sliderInput_fromfile("duration_transmission_covarying_cv", "How cumulative transmission covaries with total duration (in absence of ACF)"),

  #       plotOutput("distribution_plot")
  #   ),

  #   bslib::card(
  #     title = "Results",

  #     plotOutput("averted_plot"),
  #     tableOutput("results_table")
  # )
))

# Define server logic ----
server <- function(input, output) {

  # # Overall structure: 
  # # dalys_per_average_case depends on first block of slider values
  # averages = dalys_per_average_case(sliderValues1)
  # # within_case_avertible depends on second block of slider values
  # proportions = within_case_avertible(slidervalues2)
  # # within_case_cumulative_and_averted depends on both of the above, but not directly on slider inputs
  # within_case = within_case_cumulative_and_averted(averages = averages, proportions = proportions)
  # # between_case_differences depends on third block of slider values (duration_cv, duration_xx_covarying_cv)
  # between_case = between_case_differences(sliderValues3)
  # # daly_estimator detpends on within_case_cumulative_and_averted (and its dependencies) and between_case_differences
  # dalys_averted_per_case_detected = daly_estimator(within_case = within_case, between_case = between_case)

  # # plot functions depned on corresponding functions' outputs: 
  # # plot_averages depends on input_first_block or the resulting dalys_per_average_case
  # averages_plot = plot_averages(output_dalys_per_average_case = averages)
  # # plot_detectable_proportion depends on second input block and plot_averages output from first input block
  # proportions_plot = plot_detectable_proportion(averages_plot = averages_plot, estimates 
  #  slidervalues2)
  # # plot_time_course depends on within_case (with its dependencies on first and second input blocks) and on second input block directly
  # time_course_plot = plot_time_course(within_case = within_case, estimates = slidervalues2)
  # # plot_heterogeneity depends only on third input block parameters
  # heterogeneity_plot = plot_heterogeneity(estimates = sliderValues3)
  # # output_table depends on daly_estimator output, or on all 3 input blocks
  # output_table = output_table(dalys_averted_per_case_detected)
  # # as does plot_averted
  # averted_plot = plot_averted(dalys_averted_per_case_detected)



   # Create a reactive list of all slider input values for first block
  sliderValues1 <- reactive({
    list(
      tb_symptom_duration = input$tb_symptom_duration,
      tb_symptom_dw = input$tb_symptom_dw,
      tb_cfr = input$tb_cfr,
      tb_death_yearslost = input$tb_death_yearslost,
      posttb_symptom_duration = input$posttb_symptom_duration,
      posttb_symptom_dw = input$posttb_symptom_dw,
      posttb_cfr = input$posttb_cfr,
      posttb_death_yearslost = input$posttb_death_yearslost,
      posttb_death_timing = input$posttb_death_timing,
      discounting_rate = input$discounting_rate,
      downstream_cases = input$downstream_cases
    )
  })

  # then run the daly estimator function and plot results
  averages = dalys_per_average_case(sliderValues1())
  output$testtable = renderTable({
    as.data.frame(averages())
  })
  output$averages_plot = renderPlot({
    plot_averages(output_dalys_per_average_case = averages())
  })
  
     
  
  # # Show the result
  
  # output$average_plot <- renderPlot({
  #   plot_averages(dalyResult())
  # })
  # output$time_course_plot <- renderPlot({
  #   plot_time_course(dalyResult())
  # })
  # output$distribution_plot <- renderPlot({
  #   plot_distributions(dalyResult())
  # })
  # output$averted_plot <- renderPlot({
  #   plot_averted(dalyResult())
  # })
  # output$results_table <- renderTable({
  #   output_table(dalyResult())
  # })
}
  
# call the shiny app ----
shinyApp(ui = ui, server = server)



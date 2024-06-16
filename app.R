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
      bslib::card(
      title = "Total DALYs associated with an average TB case",
      sliderInput_fromfile("tb_symptom_duration", "Duration of symptomatic TB (years)"),
      sliderInput_fromfile("tb_symptom_dw", "Disability weight while symptomatic with TB"),
      sliderInput_fromfile("tb_cfr", "TB case fatality ratio"),
      sliderInput_fromfile("tb_death_yearslost", "Years of life lost per TB death"),
      
      sliderInput_fromfile("posttb_symptom_duration", "Time lived post-TB (years)"),
      sliderInput_fromfile("posttb_symptom_dw", "Disability weight after TB (average, including those without sequelae)", step=0.02),
      sliderInput_fromfile("posttb_cfr", "Post-TB fatality ratio (proportion who die early of TB sequelae)", step=0.01),
      sliderInput_fromfile("posttb_death_yearslost", "Years of life lost per death from TB sequelae"),
      sliderInput_fromfile("posttb_death_timing", "Mean time to death from TB sequelae (years)"),
      
      sliderInput_fromfile("discounting_rate", "Annual discounting rate on health outcomes", step=0.005),
  
      sliderInput_fromfile("downstream_cases", "# of downstream cases attributable to average case [already discounted]", step=0.25),
    # consider changing this and adding discounting to this model instead
      
      tableOutput("testtable")#,
      # plotOutput("average_plot") # need to make reactive and dependent on the above inputs
      
    ),

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
  #       sliderInput_fromfile("duration_tbdeath_multiplier", "Correlation between total duration (in absence of ACF) and mortality risk"),
  #       sliderInput_fromfile("duration_transmission_multiplier", "Correlation between total duration (in absence of ACF) and total transmission"),

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

   # Create a reactive list of all slider input values
  sliderValues <- reactiveValues()
  
  observe({
    sliderValues$tb_symptom_duration <- input$tb_symptom_duration
    sliderValues$tb_symptom_dw <- input$tb_symptom_dw
    sliderValues$tb_cfr <- input$tb_cfr
    sliderValues$tb_death_yearslost <- input$tb_death_yearslost
    sliderValues$posttb_symptom_duration <- input$posttb_symptom_duration
    sliderValues$posttb_symptom_dw <- input$posttb_symptom_dw
    sliderValues$posttb_cfr <- input$posttb_cfr
    sliderValues$posttb_death_yearslost <- input$posttb_death_yearslost
    sliderValues$posttb_death_timing <- input$posttb_death_timing
    sliderValues$discounting_rate <- input$discounting_rate
    sliderValues$downstream_cases <- input$downstream_cases

    print(sliderValues) 
  })

  # then run the daly estimator function 
  testResult <- dalys_per_average_case(sliderValues)
  testtable <- as.data.frame(testResult)
     
  
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



# to do: 
#compartmentalize sampling function so the rest runs faster, or figure out an analytical substitute
# formatting and subsetting of display, incl accordion of input subsections
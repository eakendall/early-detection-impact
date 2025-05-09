# setwd("~/Google Drive//My Drive//DALY impact of ACF 2024/ACF-impact-Rshiny/")
library(shiny)
library(bslib)
library(tidyverse)
library(plotly)

source("daly_estimator.R"); print("sourced estimator")
source("daly_plot_functions.R"); print("sourced plot functions")

paramdf <- read.csv("DALY_model_param_values.csv", header = T); print("read params")

slider_input_from_file <- function(id, label, paramtable = paramdf, step = NULL) {
  sliderInput(id, label, 
              min = round(paramtable[paramtable$param==id,"low"] ,2),
              max = round(paramtable[paramtable$param==id,"high"], 2),
              value = round(paramtable[paramtable$param==id,"mid"], 2),
              step = step)
}

.boot_dep <- "https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.0.0/css/bootstrap.min.css"

ui <- bslib::page_navbar(
  title = "Early TB detection DALY estimator",

  # change font sizes/formatting of all sliders
  tags$head(
    tags$style(HTML("
      .irs-bar {width: 100%; height: 25px; background: black; border-top: 1px solid black; border-bottom: 1px solid black;}
      .irs-bar-edge {background: black; border: 1px solid black; height: 25px; border-radius: 0px; width: 20px;}
      .irs-line {border: 1px solid black; height: 25px; border-radius: 0px;}
      .irs-grid-text {font-family: 'arial'; color: white; bottom: 17px; z-index: 1;}
      .irs-grid-pol {display: none;}
      .irs-max {font-family: 'arial'; color: black;}
      .irs-min {font-family: 'arial'; color: black;}
      .irs-single {color:black; background:#6666ff;}
      .irs-slider {width: 30px; height: 30px; top: 22px;}
    ")),
    tags$style('body {
      font-family: Arial; 
      font-size: 12px; 
      # font-style: italic; }'
    )
  ),

  collapsible = TRUE,

  tabPanel("Total DALYs per average TB episode",
    # organize as 4 navbarPages, for the three input sections plus the final results
    # within each tabPanel, include a results card
    layout_columns(
      col_widths = c(9,3),
      card(
        sidebarLayout(
          sidebarPanel(
            width = 5,
            style = "height: 90vh; overflow-y: auto;",
            accordion(
              id = "average case inputs",
              # multiple = TRUE,
              accordion_panel(
                title = "Inputs: TB morbidity",
                slider_input_from_file("tb_symptom_duration", "Duration of symptoms (years) per TB episode", 
                  step = 0.02),
                slider_input_from_file("tb_symptom_dw", "Disability weight while symptomatic")
              ),
              accordion_panel(
                title = "Inputs: TB mortality",
                slider_input_from_file("tb_cfr", "TB case fatality ratio"),
                slider_input_from_file("tb_death_yearslost", "Years of life lost per TB death")
              ),
              accordion_panel(
                title = "Inputs: Sequelae",
                slider_input_from_file("posttb_symptom_duration", "Time lived after TB (years)"),
                slider_input_from_file("posttb_symptom_dw", "Average disability weight after TB", step = 0.02),
                slider_input_from_file("posttb_cfr", "Post-TB fatality ratio 
                                      (proportion of all TB survivors who die early of TB sequelae)", step = 0.01),
                slider_input_from_file("posttb_death_yearslost",
                                      "Years of life lost per death from TB sequelae"),
                slider_input_from_file("posttb_death_timing", "Mean time to death from TB sequelae (years)")
              ),
              accordion_panel(
                title = "Inputs: Transmission",
                slider_input_from_file("downstream_cases",
                                       "# of attributable downstream cases", step = 0.25),
                slider_input_from_file("downstream_timing",
                                       "Mean years to downstream cases")
              ),
              accordion_panel(
              title = "Inputs: Temporal discounting",
              slider_input_from_file("discounting_rate", "Annual discounting rate",
                                     step = 0.005)
              ),
            )
          ),
          mainPanel(
            width=7,
            plotOutput("averages_plot", width = "100%")
          )
        ) # end sidebarLayout
      ), # end main card
      card(
        card_header("Results summary"),
        htmlOutput("results_table_forsummary1")
      )
    ) # end layout_columns
  ), # end tabPanel
  
  tabPanel("Average vs detected TB",
    card(
      min_height = 100,
      max_height = 250,
      full_screen = TRUE,
      card_header("Overview of heterogeneity portion of model"),
      p("Individuals whose TB is preferentially detected through screening may have disease courses that (even without screening) differ from the average."),
      p("Individuals whose TB is detectable for a longer time period (prior to routine diagnosis, death, or spontaneous resolution) are more likely to be detected through screening. These same individuals may generate more cumulative transmission, whereas those with more rapidly progressive disease of shorter duration may be at higher risk for mortality."),
      p("These relationships are parametrized as covariances, after normalizing duration, transmission contribution, and mortality risk to each have a mean of 1. For illustrative purposes, they are simulated here assuming log-normal distributions of duration, transmission, and mortality. The model used for DALY impact estimation is more general, however, and is valid even for covariance values not able to be illustrated in the figure.")
    ),
    # within each tabPanel, include a results card
    layout_columns(
      col_widths = c(9,3),
      card(
        sidebarLayout(
          sidebarPanel(
            width = 5,
            style = "height: 90vh; overflow-y: auto;",
            accordion(
              accordion_panel("Inputs: Duration and Transmission",
                slider_input_from_file("covariance_transmission_duration",
                                  "Transmission-duration covariance", step = 0.1)
              ),
              accordion_panel("Inputs: Duration and mortality",
                slider_input_from_file("covariance_mortality_duration",
                                  "Mortality-duration covariance")
              )
            )
          ), # end sidebarPanel
          mainPanel(
            width = 7,
            plotOutput("heterogeneity_plot", height = "900px")
          )
        ) # end sidebarLayout
      ), # end main card
      card(
        card_header("Results summary"),
        htmlOutput("results_table_forsummary3")
      ) # end results card
    ) # end layout_columns
  ),
  tabPanel("Timing of DALY accrual",
    card(
      max_height = 250,
      full_screen = TRUE,
      card_header("Overview of DALY accrual portion of model"),
      p("Disease in the present may increase the risk of future morbidity, mortality, and secondary cases. 'Accrual' refers to when future DALYs become inevitable, even if they haven't yet occurred."),
      p("Some DALYs may accrue before TB becomes detectable through screening, or after it would be diagnosed through routine care even without screening. These will not be affected by early detection."),
      p("Within the detectable period, the rate of DALY accrual may increase over time as disease becomes more severe. Cross-sectional screening will intercept TB at a random point in its disease course and thus preferentially averts the proportion of DALYs that accrue later in the detectable period.")
    ),
    layout_columns(
      col_widths = c(9,3),
      card(
        sidebarLayout(
          sidebarPanel(
            width = 5,
            style = "height: 90vh; overflow-y: auto;",
            accordion(
              accordion_panel("Inputs: Accrual before detectability",
                slider_input_from_file("predetection_mm",
                                      "Proportion of morbidity and mortality that accrue 
                                      before TB becomes detectable by screening algorithm"),
                slider_input_from_file("predetection_transmission",
                                      "Proportion of transmission that occurs 
                                      before TB becomes detectable by screening algorithm")
              ),
              accordion_panel("Inputs: Accrual after routine diagnosis",
                slider_input_from_file("postrx_mm",
                                      "Proportion of morbidity and mortality that accrue 
                                      after routine detection"),
                slider_input_from_file("postrx_transmission",
                                      "Proportion of transmission that occurs 
                                      after routine detection")
              ),
              accordion_panel("Inputs: Timing within detectable period",
                slider_input_from_file("accrual_first_half_mm",
                                      "Of personal DALYs that accrue during detectable period, 
                                      proportion in 1st half of that period"),
                slider_input_from_file("accrual_first_half_transmission",
                                      "Of transmission that occurs during detectable period, 
                                      proportion in 1st half of that period"),
                slider_input_from_file("accrual_power_mm",
                                      "Concentration of personal DALY accrual late in detectable period (power relationship)", step = 0.1),
                slider_input_from_file("accrual_power_transmission",
                                      "Concentration of transmission late in detectable period (power relationship)", step = 0.1)
              )
            ) # end accordion
          ), # end sidebarPanel
          mainPanel(
            style = "height: 90vh; overflow-y: auto;",
            width = 7,
            plotOutput("proportions_plot"),
            plotOutput("time_course_plot")
          )
        ) # end sidebarLayout
      ), # end main card
      card(
        card_header("Results summary"),
        htmlOutput("results_table_forsummary2")
      ) # end results card
    ) # end layout_columns
  ), # end tabPanel
  tabPanel("Results - DALYs averted per early detection",
    fluidRow(
      column(width = 4,
        plotOutput("averages_plot_with_averted")
      ),
      column(width = 4,
        htmlOutput("results_table") # kable (html) table
      ),
      column(width = 4,
        plotOutput("averted_plot")
      )
    )
  )
) # end ui


# Define server logic ----
server <- function(input, output) {

  # Create a reactive list of all slider input values for first block of inputs
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
      downstream_cases = input$downstream_cases,
      downstream_timing = input$downstream_timing
    )
  })

  # # Create a reactive list of all slider input values for second block
  sliderValues2 <- reactive({
    list(
      accrual_first_half_mm = input$accrual_first_half_mm,
      accrual_first_half_transmission = input$accrual_first_half_transmission,
      accrual_power_mm = input$accrual_power_mm,
      accrual_power_transmission = input$accrual_power_transmission,
      predetection_mm = input$predetection_mm,
      predetection_transmission = input$predetection_transmission,
      postrx_mm = input$postrx_mm,
      postrx_transmission = input$postrx_transmission
    )
  })

  # Create a reactive list of all slider input values for third block
  sliderValues3 <- reactive({
    list(
      covariance_transmission_duration = input$covariance_transmission_duration,
      covariance_mortality_duration = input$covariance_mortality_duration
    )
  })

  # Define the intermediate reactive functions
  averages <- reactive(dalys_per_average_case(sliderValues1()))
  proportions_reactive <- reactive(within_case_avertible(sliderValues2()))
  
  within_case <- reactive(within_case_cumulative_and_averted(
    averages = averages(),
    proportions = proportions_reactive()))
  
  between_case <- reactive(between_case_differences(estimates = sliderValues3()))

  dalys_averted_per_case_detected <- reactive(
    daly_estimator(within_case = within_case(),
                   between_case = between_case()))

  # plot functions depend on corresponding functions' outputs:
  averages_plot_reactive <- reactive(
    plot_averages(output_dalys_per_average_case = averages()))

  output$averages_plot <- renderPlot(averages_plot_reactive())

  output$proportions_plot = renderPlot(
    plot_detectable_proportion(averages_plot = averages_plot_reactive(),
                               estimates = sliderValues2()))

  output$time_course_plot <- renderPlot(
    plot_time_course(within_case = within_case(),
                     estimates = sliderValues2()))

  output$heterogeneity_plot <- renderPlot(
    {plot_heterogeneity(estimates = sliderValues3())}, height = 600)

  # both the results tab plots should have the same ymax

  output$results_table <- renderText({
    output_table(dalys_averted_per_case_detected(), forsummary = 0)})

    output$results_table_forsummary1 <- renderText({
    output_table(dalys_averted_per_case_detected(), forsummary = 0)})

    output$results_table_forsummary2 <- renderText({
    output_table(dalys_averted_per_case_detected(), forsummary = 0)})

    output$results_table_forsummary3 <- renderText({
    output_table(dalys_averted_per_case_detected(), forsummary = 0)})

  output$averages_plot_with_averted <- renderPlot(
    plot_averted_portion(dalys_averted_per_case_detected(),
    base = "average",
    ymax = max(dalys_averted_per_case_detected() %>%
    filter((cumulative_or_averted == "cumulative" &
      (average_or_detected == "average") |  average_or_detected == "detected")) %>%
     group_by(average_or_detected) %>%
     summarise(y = sum(value)) %>%
     pull(y))))

  output$averted_plot <- renderPlot(
    plot_averted_portion(dalys_averted_per_case_detected(),
    base = "detected",
    ymax = max(dalys_averted_per_case_detected() %>%
      filter((cumulative_or_averted == "cumulative" &
      (average_or_detected == "average") |  average_or_detected == "detected")) %>%
      group_by(average_or_detected) %>%
      summarise(y = sum(value)) %>%
      pull(y))))
}


print("read to call app")
# call the shiny app ----
shinyApp(ui = ui, server = server)

# setwd("~/Google Drive//My Drive//DALY impact of ACF 2024/ACF-impact-Rshiny/")
library(shiny)
library(bslib)
library(tidyverse)
library(plotly)

source("daly_estimator.R")
source("daly_plot_functions.R")

paramdf <- read.csv("DALY_model_param_values.csv", header = T)

slider_input_from_file <- function(id, label, paramtable = paramdf, step = NULL) {
  sliderInput(id, label, 
              min = round(paramtable[paramtable$param==id,"low"] ,2),
              max = round(paramtable[paramtable$param==id,"high"], 2),
              value = round(paramtable[paramtable$param==id,"mid"], 2),
              step = step)
}

.boot_dep <- "https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.0.0/css/bootstrap.min.css"

ui <- bslib::page_navbar( 
  title = "ACF Impact",

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
    )),

tabPanel("Total DALYs per average TB case",
  # organize as 4 navbarPages, for the three input sections plus the final results
  # within each tabPanel, include a results card
  layout_columns(
      col_widths = c(9,3),
      card(
        sidebarLayout(
          sidebarPanel(
            width = 5,
            # style = "height: 90vh; overflow-y: auto;",
            accordion(
              id = "average case inputs",
              # multiple = TRUE,
              accordion_panel(
                title = "TB morbidity",
                slider_input_from_file("tb_symptom_duration", "Duration of symptomatic TB (years)"),
                slider_input_from_file("tb_symptom_dw", "Disability weight while symptomatic")
              ),
              accordion_panel(
                title = "TB mortality",
                slider_input_from_file("tb_cfr", "TB case fatality ratio"),
                slider_input_from_file("tb_death_yearslost", "Years of life lost per TB death")
              ),
              accordion_panel(
                title = "Sequelae",
                slider_input_from_file("posttb_symptom_duration", "Time lived with TB sequelae (years)"),
                slider_input_from_file("posttb_symptom_dw", "Disability weight during TB sequelae", step = 0.02),
                slider_input_from_file("posttb_cfr", "Post-TB fatality ratio 
                                      (proportion of all TB survivors who die early of TB sequelae)", step = 0.01),
                slider_input_from_file("posttb_death_yearslost",
                                      "Years of life lost per death from TB sequelae"),
                slider_input_from_file("posttb_death_timing", "Mean time to death from TB sequelae (years)")
              ),
              accordion_panel(
                title = "Temporal discounting",
                slider_input_from_file("discounting_rate", "Annual discounting rate on health outcomes",
                                      step = 0.005)
              ),
              accordion_panel(
                title = "Transmission",
                slider_input_from_file("downstream_cases",
                                      "# of attributable downstream cases 
                                      [already discounted]", step = 0.25)
              ) #! consider changing this and adding discounting to this model instead
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
        htmlOutput("results_table_forsummary")
      )
    ) # end layout_columns
  ), # end tabPanel
  tabPanel("Timing of DALY accrual",
    sidebarLayout(
      sidebarPanel(
        style = "height: 90vh; overflow-y: auto;",
        accordion(
          accordion_panel("Accrual before detectability",
            slider_input_from_file("predetection_mm",
                                  "Proportion of morbidity and mortality that accrue 
                                  before TB becomes detectable by screening algorithm"),
            slider_input_from_file("predetection_transmission",
                                  "Proportion of transmission that occurs 
                                  before TB becomes detectable by screening algorithm")
          ),
          accordion_panel("Accrual after routine diagnosis",
            slider_input_from_file("postrx_mm",
                                  "Proportion of morbidity and mortality that accrue 
                                  after routine detection"),
            slider_input_from_file("postrx_transmission",
                                  "Proportion of transmission that occurs 
                                  after routine detection")
          ),
          accordion_panel("Timing within detectable period",
            slider_input_from_file("second_half_vs_first_mm",
                                  "Of personal DALYs that accrue during detectable period, 
                                  proportion in 2nd half"),
            slider_input_from_file("second_half_vs_first_transmission",
                                  "Of transmission that occurs during detectable period, 
                                  proportion in 2nd half"),
            slider_input_from_file("resolving_detectable",
                                  "Proportion of TB that will spontaneously resolve 
                                  after becoming detectable")
          )
        )
      ),
      mainPanel(
        style = "height: 90vh; overflow-y: auto;",
        plotOutput("proportions_plot"),
        plotOutput("time_course_plot")
      )
    )
  ),
  tabPanel("Average vs detected cases",
    sidebarLayout(
      sidebarPanel(
        accordion(
          accordion_panel("Variance in disease duration",
            slider_input_from_file("duration_cv", "Coefficient of variation in TB duration (in absence of ACF)")
          ),
          accordion_panel("Covariance of outcomes with duration",
            slider_input_from_file("duration_tbdeath_covarying_cv",
                              "How mortality risk co-varies with duration"),
            slider_input_from_file("duration_transmission_covarying_cv",
                              "How cumulative transmission co-varies with duration")
          )
        )
      ),
      mainPanel(
        plotOutput("heterogeneity_plot", height = "600px")
      )
    )
  ),
  tabPanel("Results - DALYs averted per case detected",
    fluidRow(
      column(width = 4,
        plotOutput("averages_plot_ymax")
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
      downstream_cases = input$downstream_cases
    )
  })

  # Create a reactive list of all slider input values for second block
  sliderValues2 <- reactive({
    list(
      second_half_vs_first_mm = input$second_half_vs_first_mm,
      second_half_vs_first_transmission = input$second_half_vs_first_transmission,
      resolving_detectable = input$resolving_detectable,
      predetection_mm = input$predetection_mm,
      predetection_transmission = input$predetection_transmission,
      postrx_mm = input$postrx_mm,
      postrx_transmission = input$postrx_transmission
    )
  })

  # Create a reactive list of all slider input values for third block
  sliderValues3 <- reactive({
    list(
      duration_cv = input$duration_cv,
      duration_tbdeath_covarying_cv = input$duration_tbdeath_covarying_cv,
      duration_transmission_covarying_cv = input$duration_transmission_covarying_cv
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
    output_table(dalys_averted_per_case_detected())})

    output$results_table_forsummary <- renderText({
    output_table(dalys_averted_per_case_detected(),forsummary = TRUE)})

  output$averages_plot_ymax <- renderPlot(
    plot_averages(output_dalys_per_average_case = averages(), 
    ymax = max(dalys_averted_per_case_detected() %>% 
    filter((cumulative_or_averted == "cumulative" & average_or_detected == "average") |
     (cumulative_or_averted == "averted" & average_or_detected == "detected") ) %>%
     group_by(cumulative_or_averted) %>%
     summarise(y = sum(value)) %>%
     pull(y)) ))

  output$averted_plot <- renderPlot(
    plot_averted(dalys_averted_per_case_detected(),
    ymax = max(dalys_averted_per_case_detected() %>% 
    filter((cumulative_or_averted == "cumulative" & average_or_detected == "average") |
     (cumulative_or_averted == "averted" & average_or_detected == "detected") ) %>%
     group_by(cumulative_or_averted) %>%
     summarise(y = sum(value)) %>%
     pull(y))))

  
}

# call the shiny app ----
shinyApp(ui = ui, server = server)

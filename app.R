library(shiny)
library(bslib)

source("daly_estimator.R")

paramdf <- read.csv("DALY_model_param_values.csv", header = T)

sliderInput_fromfile <- function(id, label, paramtable=paramdf, step=NULL) {
  sliderInput(id, label, 
              min = round(paramtable[paramtable$param==id,"low"],2), 
              max = round(paramtable[paramtable$param==id,"high"],2), 
              value = round(paramtable[paramtable$param==id,"mid"],2),
              step=step)
}

##*** TO DO:
# add differentiation into more and less important parameters, 
# with user option whether to see the "minor" parameters

# Show figures illustrating steps on each tab

# Results tab

# Table of sub-estimates, e.g. components of average-case DALYs, proportion averetible on average, and case_dalys_mortality_multiplier and case_dalys_transmission_multiplier


ui <- fluidPage(
  
  # App title ----
  titlePanel("DALY impact of early case detection"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    sidebarPanel(
      tabsetPanel(
        tabPanel("Total DALYs associated with the average TB case",
          sliderInput_fromfile("tb_symptom_duration", "Duration (years) of symptoms"),
          sliderInput_fromfile("tb_symptom_dw", "Disability weight while symptomatic with TB"),
          sliderInput_fromfile("tb_cfr", "TB case fatality ratio"),
          sliderInput_fromfile("tb_death_yearslost", "Years of life lost per TB death"),
          
          sliderInput_fromfile("posttb_symptom_duration", "Duration (years) lived post-TB"),
          sliderInput_fromfile("posttb_symptom_dw", "Disability weight after TB (average, including those without sequelae)", step=0.02),
          sliderInput_fromfile("posttb_cfr", "Post-TB fatality ratio (proportion who die early of TB sequelae)", step=0.01),
          sliderInput_fromfile("posttb_death_yearslost", "Years of life lost per death from TB sequelae"),
          sliderInput_fromfile("posttb_death_timing", "Mean time (years) to death from TB sequelae"),
          
          sliderInput_fromfile("discounting_rate", "Annual discounting rate on health outcomes", step=0.005),
      
          sliderInput_fromfile("downstream_cases", "Downstream cases attributable to average case [already discounted]", step=0.25),
        # consider changing this and adding discounting to this model instead
        
        ### may want to break into subgroups. 
      ),
    
        tabPanel("Timing of DALY accrual versus detection through screening",
          sliderInput_fromfile("second_half_vs_first_mm", "Of morbidity, mortality, and post-TB disease that accrue during detectable period, proportion in 2nd half"),
          sliderInput_fromfile("second_half_vs_first_transmission", "Of transmission that occurs during detectable period, proportion in 2nd half"),
          sliderInput_fromfile("resolving_detectable", "Proportion of TB that will spontaneously resolve after becoming detectable"),
    
          sliderInput_fromfile("predetection_mm", "Proportion of morbidity and mortality that accrue before TB becomes detectable by screening algorithm"),
          sliderInput_fromfile("predetection_transmission", "Proportion of transmission that occurs before TB becomes detectable by screening algorithm"),
          sliderInput_fromfile("postrx_mm", "Proportion of morbidity and mortality that accrue after routine detection"),
          sliderInput_fromfile("postrx_transmission", "Proportion of transmission that occurs after routine detection")
        ),
      
        tabPanel("How TB found through screening differs from the average case",
            sliderInput_fromfile("duration_cv", "Coefficient of variation in TB duration"),
            sliderInput_fromfile("duration_tbdeath_multiplier", "Correlation between total duration (in absence of ACF) and mortality risk"),
            sliderInput_fromfile("duration_transmission_multiplier", "Correlation between total duration (in absence of ACF) and total transmission")
        )
      )
    ),

    mainPanel(
    
      textOutput("result"),
      tableOutput("table") # will put intermediate results here
      ## Will add figure(s) here too??
      
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  
  # Reactive expression to create data frame of all input values ----
  dalyResult <- reactive({
    daly_estimator(input)
    })
  
  intermediateResults <- reactive({
    daly_estimator(input)
  })
   
  
  # Show the result
  ## Need to make these prettier
  output$result <- renderPrint(dalyResult())
  output$table <- renderTable(intermediateResults())
  
}
  
# call the shiny app ----
shinyApp(ui = ui, server = server)

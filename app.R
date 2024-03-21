library(shiny)
library(brms)
library(dplyr)
library(DT)
library(ggplot2)
library(gridExtra)

ui <- fluidPage(
  titlePanel("Bayesian Statistics for Seroepidemiology"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose CSV File", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      selectInput("numerator", "Numerator:", choices = NULL),
      selectInput("denominator", "Denominator:", choices = NULL),
      selectInput("randomEffect", "Random Effect:", choices = NULL),
      selectInput("locusSubset", "Locus Subset:", choices = NULL, multiple = TRUE),
      numericInput("chains", "Number of Chains:", value = 2, min = 1),
      numericInput("iterations", "Iterations per Chain:", value = 2000, min = 500),
      actionButton("runModel", "Run Model")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Model Results", downloadButton("downloadData", "Download Results"), dataTableOutput("modelSummary")),
        tabPanel("Density Plot", plotOutput("densityPlot"))
      )
    )
  ),
  div(style = "position: fixed; bottom: 0; width: 100%; text-align: center; background-color: lightgrey;",
      "Version 0.0.1, Last update 27th February 2024.")
)

server <- function(input, output, session) {
  
  results <- reactiveVal(data.frame())
  
  uploadedData <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    read.csv(inFile$datapath)
  })
  
  observe({
    df <- uploadedData()
    if (is.null(df)) return()
    updateSelectInput(session, "numerator", choices = names(df))
    updateSelectInput(session, "denominator", choices = names(df))
    updateSelectInput(session, "randomEffect", choices = names(df))
    updateSelectInput(session, "locusSubset", choices = unique(df$locus), selected = unique(df$locus)[1])
  })
  
  plotData <- reactiveVal(data.frame(locus = character(), value = numeric(), median = numeric(), lower = numeric(), upper = numeric()))
  
  observeEvent(input$runModel, {
    showModal(modalDialog(title = "Processing", "The models are now running. This may take some time.", easyClose = TRUE, footer = NULL))
    
    df <- uploadedData()
    if (is.null(df)) {
      removeModal()
      return()
    }
    
    tmpPlotData <- data.frame(locus = character(), value = numeric(), median = numeric(), lower = numeric(), upper = numeric())
    
    aggregateResults <- list()
    
    for (locusValue in input$locusSubset) {
      df_subset <- df %>% filter(locus == locusValue)
      
      if (nrow(df_subset) > 0 && input$numerator %in% names(df_subset) && input$denominator %in% names(df_subset) && input$randomEffect %in% names(df_subset)) {
        formula <- as.formula(paste(input$numerator, "| trials(", input$denominator, ") ~ (1 |", input$randomEffect, ")"))
        fit <- brm(formula, data = df_subset, family = binomial(), chains = input$chains, iter = input$iterations, warmup = floor(input$iterations / 4), control = list(adapt_delta = 0.95))
        
        post <- posterior_samples(fit)
        postTransformed <- exp(post[,1])/(1 + exp(post[,1]))
        
        medianVal <- quantile(postTransformed, probs = 0.5)
        lowerCI <- quantile(postTransformed, probs = 0.025)
        upperCI <- quantile(postTransformed, probs = 0.975)
        
        tmpPlotData <- rbind(tmpPlotData, data.frame(locus = locusValue, value = postTransformed, median = medianVal, lower = lowerCI, upper = upperCI))
        
        summary_tab <- data.frame(
          locus = locusValue,
          type = "Adjusted",
          geog = "global",
          est = medianVal,
          upper = upperCI,
          lower = lowerCI
        )
        
        aggregateResults[[length(aggregateResults) + 1]] <- summary_tab
      }
    }
    
    results(do.call(rbind, aggregateResults))
    plotData(tmpPlotData)
    
    removeModal()
  })
  
  output$modelSummary <- DT::renderDataTable({
    datatable(results(), options = list(pageLength = 5, autoWidth = TRUE), class = 'cell-border stripe')
  })
  
  output$densityPlot <- renderPlot({
    tmpPlotData <- plotData()
    
    plots <- lapply(unique(tmpPlotData$locus), function(locus) {
      df_subset <- tmpPlotData[tmpPlotData$locus == locus, ]
      p <- ggplot(df_subset, aes(x = value)) +
        geom_density(alpha = 0.75, fill = "#56B4E9") +
        geom_vline(aes(xintercept = median), color = "red", linetype = "dashed") +
        geom_vline(aes(xintercept = lower), color = "blue", linetype = "dashed") +
        geom_vline(aes(xintercept = upper), color = "blue", linetype = "dashed") +
        labs(title = paste("Posterior Distribution for", locus), x = "Transformed Posterior Samples") +
        theme_minimal() +
        xlim(0, 0.15)
      return(p)
    })
    
    do.call(gridExtra::grid.arrange, c(plots, ncol = 1))
  }, res = 96)
  
  output$downloadData <- downloadHandler(
    filename = function() { "model_results.csv" },
    content = function(file) {
      write.csv(results(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)

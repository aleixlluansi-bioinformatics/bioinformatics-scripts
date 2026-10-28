# ============================================================
# Shiny App: Data Distribution Assessment
# ============================================================
# Features:
#   1. Upload CSV files
#   2. Preview dataset
#   3. Descriptive statistics
#   4. Histogram and boxplot
#   5. Missing-value assessment
#   6. Shapiro-Wilk normality test
#   7. Normal Q-Q plot
#
# Run:
#   shiny::runApp()
# ============================================================


# 1. Load packages --------------------------------------------

library(shiny)
library(ggplot2)


# 2. User interface -------------------------------------------

ui <- fluidPage(
  
  titlePanel("Data Distribution Assessment"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      fileInput(
        "file",
        "Upload a CSV file",
        accept = ".csv"
      ),
      
      helpText(
        "Upload a dataset containing numeric and/or categorical variables."
      ),
      
      selectInput(
        "variable",
        "Select a numeric variable",
        choices = NULL
      ),
      
      checkboxInput(
        "remove_na",
        "Remove missing values from analysis",
        value = TRUE
      )
    ),
    
    mainPanel(
      
      h3("Dataset preview"),
      tableOutput("data_preview"),
      
      h3("Descriptive statistics"),
      tableOutput("descriptive_stats"),
      
      h3("Distribution histogram"),
      plotOutput("histogram"),
      
      h3("Boxplot"),
      plotOutput("boxplot"),
      
      h3("Normality assessment"),
      
      p(
        "The Shapiro-Wilk test evaluates the null hypothesis ",
        "that the observations come from a normal distribution."
      ),
      
      tableOutput("normality_test"),
      
      verbatimTextOutput("normality_interpretation"),
      
      h4("Normal Q-Q plot"),
      plotOutput("qq_plot"),
      
      h3("Missing values"),
      tableOutput("missing_values")
    )
  )
)


# 3. Server ----------------------------------------------------

server <- function(input, output, session) {
  
  # Read uploaded CSV file.
  uploaded_data <- reactive({
    
    req(input$file)
    
    read.csv(
      input$file$datapath,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  
  
  # Update numeric variable selector.
  observeEvent(uploaded_data(), {
    
    data <- uploaded_data()
    
    numeric_columns <- names(data)[
      vapply(data, is.numeric, logical(1))
    ]
    
    updateSelectInput(
      session,
      "variable",
      choices = numeric_columns,
      selected = if (length(numeric_columns) > 0) {
        numeric_columns[1]
      } else {
        character(0)
      }
    )
  })
  
  
  # Prepare the selected numeric variable.
  selected_values <- reactive({
    
    data <- uploaded_data()
    
    req(input$variable)
    
    validate(
      need(
        input$variable %in% names(data),
        "Select a valid numeric variable."
      )
    )
    
    x <- data[[input$variable]]
    
    validate(
      need(
        is.numeric(x),
        "The selected variable must be numeric."
      )
    )
    
    if (input$remove_na) {
      x <- x[!is.na(x)]
    }
    
    validate(
      need(
        length(x) > 0,
        "No observations available."
      )
    )
    
    x
  })
  
  
  # Dataset preview.
  output$data_preview <- renderTable({
    
    head(uploaded_data(), 10)
    
  }, rownames = TRUE)
  
  
  # Descriptive statistics.
  output$descriptive_stats <- renderTable({
    
    x <- selected_values()
    
    data.frame(
      Statistic = c(
        "Number of observations",
        "Missing values",
        "Mean",
        "Standard deviation",
        "Median",
        "Minimum",
        "Maximum",
        "First quartile",
        "Third quartile"
      ),
      
      Value = c(
        length(x),
        sum(is.na(uploaded_data()[[input$variable]])),
        mean(x),
        if (length(x) > 1) sd(x) else NA_real_,
        median(x),
        min(x),
        max(x),
        unname(quantile(x, 0.25)),
        unname(quantile(x, 0.75))
      )
    )
    
  }, digits = 3)
  
  
  # Histogram.
  output$histogram <- renderPlot({
    
    x <- selected_values()
    
    ggplot(
      data.frame(Value = x),
      aes(x = Value)
    ) +
      geom_histogram(
        bins = 20,
        color = "black"
      ) +
      labs(
        title = paste("Distribution of", input$variable),
        x = input$variable,
        y = "Frequency"
      ) +
      theme_minimal()
    
  })
  
  
  # Boxplot.
  output$boxplot <- renderPlot({
    
    x <- selected_values()
    
    ggplot(
      data.frame(Value = x),
      aes(y = Value)
    ) +
      geom_boxplot() +
      labs(
        title = paste("Boxplot of", input$variable),
        y = input$variable,
        x = NULL
      ) +
      theme_minimal()
    
  })
  
  
  # Shapiro-Wilk normality test.
  output$normality_test <- renderTable({
    
    x <- selected_values()
    
    # Shapiro-Wilk requires between 3 and 5000 observations.
    validate(
      need(
        length(x) >= 3,
        "At least 3 observations are required."
      ),
      need(
        length(x) <= 5000,
        paste(
          "Shapiro-Wilk accepts at most 5000 observations.",
          "Consider a representative sample or another assessment."
        )
      ),
      need(
        length(unique(x)) > 1,
        "The variable has no variation; normality cannot be assessed."
      )
    )
    
    test <- shapiro.test(x)
    
    data.frame(
      Test = "Shapiro-Wilk",
      N = length(x),
      W = unname(test$statistic),
      P_value = test$p.value
    )
    
  }, digits = 4)
  
  
  # Interpret the Shapiro-Wilk test.
  output$normality_interpretation <- renderText({
    
    x <- selected_values()
    
    validate(
      need(length(x) >= 3, "Insufficient observations."),
      need(length(x) <= 5000, "More than 5000 observations."),
      need(
        length(unique(x)) > 1,
        "The variable has no variation."
      )
    )
    
    test <- shapiro.test(x)
    
    if (test$p.value < 0.05) {
      
      paste0(
        "P = ", signif(test$p.value, 3),
        ". The null hypothesis of normality is rejected ",
        "at the 0.05 significance level."
      )
      
    } else {
      
      paste0(
        "P = ", signif(test$p.value, 3),
        ". There is insufficient evidence to reject ",
        "the null hypothesis of normality at the 0.05 level. ",
        "This does not prove that the data are normally distributed."
      )
      
    }
    
  })
  
  
  # Normal Q-Q plot.
  output$qq_plot <- renderPlot({
    
    x <- selected_values()
    
    ggplot(
      data.frame(Value = x),
      aes(sample = Value)
    ) +
      stat_qq() +
      stat_qq_line() +
      labs(
        title = paste("Normal Q-Q plot:", input$variable),
        x = "Theoretical quantiles",
        y = "Sample quantiles"
      ) +
      theme_minimal()
    
  })
  
  
  # Missing values by variable.
  output$missing_values <- renderTable({
    
    data <- uploaded_data()
    
    data.frame(
      Variable = names(data),
      
      Missing_values = vapply(
        data,
        function(x) sum(is.na(x)),
        integer(1)
      ),
      
      Missing_percentage = round(
        vapply(
          data,
          function(x) mean(is.na(x)) * 100,
          numeric(1)
        ),
        2
      )
    )
    
  })
  
}


# 4. Launch application ----------------------------------------

shinyApp(
  ui = ui,
  server = server
)

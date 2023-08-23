# Load required libraries
library(shiny)
library(ggplot2)
library(shinythemes)
library(plotly)

# Condition variables
cond_vars <- c("groupcutoffs", "subcutoffs", "rmsea.cuts", "srmr.cuts",
               "nnfi.cuts", "cfi.cuts", "n.excellent")
summary_vars <- c("diff_adj_sum_mean_i", "mean_diff_cent_i")

# Data for testing
comp_test <- comp_pers %>% 
  select(diff_adj_sum_mean_i, mean_diff_cent_i, all_of(cond_vars)) %>% 
  slice(1:50)
comp_test2 <- comp_pers %>% 
  select(diff_adj_sum_mean_i, mean_diff_cent_i, all_of(cond_vars)) %>% 
  slice(51:100)

saveRDS(comp_test, here("Network-Multiverse/data/comp_test_pers.RDS"))
saveRDS(comp_test2, here("Network-Multiverse/data/comp_test_emot.RDS"))

# Define the UI for the Shiny app
ui <- navbarPage(
  title = "Network Multiverse",
  theme = shinytheme("lumen"),
 
  # Dataset Selection
  tabPanel("Select Dataset",
           sidebarLayout(
             sidebarPanel(
               selectInput("dataset", 
                           label = "Choose a Dataset", 
                           choices = c("Personality",
                                       "Emotion"))
             )
           ,
           mainPanel(
             tableOutput("dataset_summary"),
             # tableOutput("table")
           )
           )
           ),
  
  # Compute summary statistics
  tabPanel("Summary Statistics",
           sidebarLayout(
             sidebarPanel(
               selectInput("groupby",
                           label = "Group by:",
                           choices = cond_vars,
                           selected = "None"),
               selectInput("filter",
                           label = "Filter:",
                           choices = cond_vars,
                           selected = "None")

               )
             ,
             mainPanel(
               verbatimTextOutput("summary"),
               # tableOutput("table")
             )
           )
           ),
  
  # Boxplots
  tabPanel("Boxplot",
           sidebarLayout(
             sidebarPanel(
               selectInput("boxplot_column",
                           label = "Select Column:",
                           choices = summary_vars,
                           selected = summary_vars[1])
             ),
             mainPanel(
               plotlyOutput("boxplot")
             )
           )
  )


)

# Define the server for the Shiny app
server <- function(input, output, session) {
  
  # Load the selected dataset based on user input
  dataset <- reactive({
    dataset_name <- input$dataset
    
    if(dataset_name == "Personality") {
      return(readRDS("data/comp_test_pers.RDS"))
    } else if (dataset_name == "Emotion") {
      return(readRDS("data/comp_test_emot.RDS"))
    }
  })
  
  output$dataset_summary <- renderTable({
    # Use a reactive expression by calling it like a function
    summary(dataset())
  })
  
  # Create ggplot boxplot
  output$boxplot <- renderPlotly({
    selected_column <- input$boxplot_column
    
    # Generate ggplot boxplot
    boxplot_plot <- ggplot(dataset(), aes(x = "", y = .data[[selected_column]])) +
      geom_boxplot() +
      labs(title = paste("Boxplot of", selected_column))
    
    # Convert ggplot plot to a Plotly plot
    ggplotly(boxplot_plot)
  })

 
}

# Run the Shiny app
shinyApp(ui, server)




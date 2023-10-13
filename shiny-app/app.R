#--- Load required libraries
library(shiny)
library(waiter)
library(ggplot2)
library(shinythemes)
library(plotly)
library(RColorBrewer)
library(tidyverse)
library(here)
library(fontawesome)#
library(knitr)
library(kableExtra)

#--- Variable names
cond_vars <- c("groupcutoffs", "subcutoffs", "rmsea.cuts", "srmr.cuts",
               "nnfi.cuts", "cfi.cuts", "n.excellent")
summary_vars <- c("Heterogeneity", 
                  # "VI", "ARI", "Modularity",
                  "Nonzero Edge Diff.", "Adjacency Diff.", "Centrality Diff.", 
                  "Density Temp. Abs. Diff.", "Density Cont. Abs. Diff.")

# --- Create actual data sets
# Personality
# comp_pers <- readRDS(here("output/mv_personality/comp_pers_upd.RDS"))
# 
# comp_pers_app <- comp_pers %>%
#     dplyr::rename(Heterogeneity = heterogeneity_g,
#                   VI = vi,
#                   ARI = ari,
#                   Modularity = modularity,
#                   "Nonzero Edge Diff." = mean_nonzero_diff_edge_i,
#                   "Adjacency Diff." = diff_adj_sum_mean_i,
#                   "Centrality Diff." = mean_diff_cent_i,
#                   "Density Temp. Abs. Diff." = mean_abs_diff_dens_temp_i,
#                   "Density Cont. Abs. Diff." = mean_abs_diff_dens_cont_i) %>%
#     tidyr::unnest(Heterogeneity) %>%
#     dplyr::select(all_of(summary_vars), all_of(cond_vars))
# 
# comp_emot <- readRDS(here("output/mv_emotion/comp_emot_upd.RDS"))
# 
# comp_emot_app <- comp_emot %>%
#   dplyr::rename(Heterogeneity = heterogeneity_g,
#                 VI = vi,
#                 ARI = ari,
#                 Modularity = modularity,
#                 "Nonzero Edge Diff." = mean_nonzero_diff_edge_i,
#                 "Adjacency Diff." = diff_adj_sum_mean_i,
#                 "Centrality Diff." = mean_diff_cent_i,
#                 "Density Temp. Abs. Diff." = mean_abs_diff_dens_temp_i,
#                 "Density Cont. Abs. Diff." = mean_abs_diff_dens_cont_i) %>%
#   tidyr::unnest(Heterogeneity) %>%
#   dplyr::select(all_of(summary_vars), all_of(cond_vars))
# 
# 
# saveRDS(comp_pers_app, here("shiny-app/data/comp_pers_app.RDS"))
# saveRDS(comp_emot_app, here("shiny-app/data/comp_emot_app.RDS"))


# # --- Data for testing
# comp_pers <- readRDS(here("output/mv_personality/comp_pers.RDS"))
# # nicer column names
# # TODO will have to do this before I save the datasets for the app
# # names(comp_pers)
# comp_pers <- comp_pers %>%
#   dplyr::rename(Heterogeneity = heterogeneity_g,
#                 VI = vi,
#                 ARI = ari,
#                 Modularity = modularity,
#                 "Nonzero Edge Diff." = mean_nonzero_diff_edge_i,
#                 "Adj. Diff." = diff_adj_sum_mean_i,
#                 "Cent. Diff." = mean_diff_cent_i,
#                 "Abs. Diff. Density Temp." = mean_abs_diff_dens_temp_i,
#                 "Abs. Diff. Density Cont." = mean_abs_diff_dens_cont_i) %>% 
#   unnest(Heterogeneity)
# 
# comp_test <- comp_pers %>%
#   select(all_of(summary_vars), all_of(cond_vars)) %>%
#   slice_sample(n = 1000)
# comp_test2 <- comp_pers %>%
#   select(all_of(summary_vars), all_of(cond_vars)) %>%
#   slice_sample(n = 1000)
# saveRDS(comp_test, here("shiny-app/data/comp_test_pers.RDS"))
# saveRDS(comp_test2, here("shiny-app/data/comp_test_emot.RDS"))

#--- Plots data set
# Personality
# mv_res_pers <- readRDS(here("output/mv_personality/mv_res_pers2.RDS"))
# 
# mv_res_plots <- lapply(mv_res_pers, function(x) {
#    x[c("group_plot_paths", "path_counts")]
#  })
# 
# saveRDS(mv_res_plots, here("shiny-app/data/mv_pers_plots.RDS"))
# rm(mv_res_pers, mv_res_plots)
# 
# # Emotion
# mv_res_emot <- readRDS(here("output/mv_emotion/mv_res_emot_new51.RDS"))
# mv_res_plots <- lapply(mv_res_emot, function(x) {
#   x[c("group_plot_paths", "path_counts")]
# })
# 
# saveRDS(mv_res_plots, here("shiny-app/data/mv_emot_plots.RDS"))
# rm(mv_res_emot, mv_res_plots)


#--- ggplot Theming
theme_multiverse <- function() {
  # add google font
  sysfonts::font_add_google("News Cycle", "news")
  # use showtext
  showtext::showtext_auto()
  # theme
  ggplot2::theme_minimal(base_family = "news") +
    ggplot2::theme(
      # remove minor grid
      panel.grid.minor = ggplot2::element_blank(),
      # Title and Axis Texts
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2), hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1.1), hjust = 0.5),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, b = 10)),
      
      # Faceting
      strip.text = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.1), hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      # Grid
      panel.grid = ggplot2::element_line(colour = "#F3F4F5"),
      # Legend
      legend.title = ggplot2::element_text(face = "plain"),
      legend.position = "top",
      legend.justification = 1,
      # Panel/Facets
      panel.spacing.y = ggplot2::unit(1.5, "lines")
    )
}


#--- SCA Plot
# For lineplot
plot_outcome <- function(mv_res, 
                         var,
                         specs = NULL, # hard-coded for now
                         y_label){   
  
  var <- enquo(var)
  
  mv_res %>% 
    dplyr::mutate(variable = as.numeric(!!var)) %>% 
    dplyr::arrange(variable) %>% 
    dplyr::mutate(iteration = dplyr::row_number()) %>%
    ggplot(aes(x = .data$iteration,
               y = variable)) + 
    geom_point(size = 0.8)+
    theme_multiverse()+
    labs(x = "")
}

# For Specification Plot
plot_specification <- function(mv_res,
                               var,
                               specs = NULL){    # hard-coded for now
  
  var <- enquo(var)
  
  mv_res %>% 
    dplyr::mutate(variable = as.numeric(!!var)) %>% 
    dplyr::arrange(variable) %>% 
    dplyr::mutate(iteration = dplyr::row_number()) %>%
    dplyr::mutate(across(c(groupcutoffs, subcutoffs,
                           rmsea.cuts, srmr.cuts,
                           cfi.cuts, nnfi.cuts,
                           n.excellent), ~ as.factor(.))) %>% 
    dplyr::mutate(groupcutoffs = fct_recode(groupcutoffs, !!!setNames(as.character(group_cuts), group_levels)),
                  subcutoffs = fct_recode(subcutoffs, !!!setNames(as.character(sub_cuts), sub_levels)),
                  rmsea.cuts = fct_recode(rmsea.cuts, !!!setNames(as.character(rmsea_cuts), rmsea_levels)),
                  srmr.cuts = fct_recode(srmr.cuts, !!!setNames(as.character(srmr_cuts), srmr_levels)),
                  cfi.cuts = fct_recode(cfi.cuts, !!!setNames(as.character(cfi_cuts), cfi_levels)),
                  nnfi.cuts = fct_recode(nnfi.cuts, !!!setNames(as.character(nnfi_cuts), nnfi_levels)),
                  n.excellent = fct_recode(n.excellent, !!!setNames(as.character(n_excels), n_excels_levels))) %>% 
    tidyr::pivot_longer(cols = c(groupcutoffs,
                                 subcutoffs,
                                 rmsea.cuts,
                                 srmr.cuts,
                                 cfi.cuts,
                                 nnfi.cuts,
                                 n.excellent),
                        values_to = "value", names_to = "specification") %>%
    dplyr::mutate(specification = dplyr::case_match(specification,
                                                    "groupcutoffs" ~ "Group",
                                                    "subcutoffs" ~ "Subgroup",
                                                    "rmsea.cuts" ~ "RMSEA",
                                                    "srmr.cuts" ~ "SRMR",
                                                    "cfi.cuts" ~ "CFI",
                                                    "nnfi.cuts" ~ "NNFI",
                                                    "n.excellent" ~ "N° excellent")) %>% 
    dplyr::mutate(specification = as.factor(specification)) %>% 
    dplyr::mutate(specification = forcats::fct_relevel(specification, 
                                                       "Group", 
                                                       "Subgroup",
                                                       "RMSEA", 
                                                       "SRMR",
                                                       "CFI",
                                                       "NNFI")) %>% 
    dplyr::mutate(value = forcats::fct_relevel(value, 
                                               "liberal",
                                               "medium-liberal",
                                               "medium",
                                               "medium-strict",
                                               "strict")) %>% 
    ggplot(aes(x = .data$iteration,
               y = 1,
               color = .data$value)) + 
    geom_point(shape = 124, size = 15
               #pch='.'   #for faster plotting
    )+
    theme_multiverse()+
    scale_y_continuous(limits = c(0.99, 1.01), expand = c(0,0))+
    scale_color_manual(values = palette_full)+
    facet_wrap(specification~., 
               ncol = 1, 
               strip.position = "left")+
    labs(y = "",
         x = "Iteration",
         color = "Specification")+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y.left = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          legend.text = element_text(size = rel(1.3)),
          legend.title = element_text(size = rel(1.3)),
          legend.spacing.x = unit(0.2, "cm"),
          legend.text.align = 0,
          strip.text = element_text(size = rel(1.4)),
          axis.text.x = element_text(size = rel(1.2)),
          axis.title.x = element_text(size = rel(1.3)))
}

# Color palette
palette_3 <- colorRampPalette(brewer.pal(9, "RdBu"))(3)
palette_4 <- colorRampPalette(brewer.pal(9, "RdBu"))(4)


# Set values explicitly 
palette_full <- c(palette_4[1:2], "lightblue", palette_4[3:4])
names(palette_full) <- c("liberal", "medium-liberal", "medium", "medium-strict", "strict")

# Define desired levels for each factor
group_cuts <- c(.50, .60, .75, .80)
sub_cuts <- c(.50, .60, .75, .80)
rmsea_cuts <- c(.03, .05, .08)
srmr_cuts <- c(.03, .05, .08)
nnfi_cuts <- c(.90, .95, .97)
cfi_cuts <- c(.90, .95, .97)
n_excels <- c(1, 2, 3)

group_levels <- c("liberal", "medium-liberal", "medium-strict", "strict")
sub_levels <- c("liberal", "medium-liberal", "medium-strict", "strict")
rmsea_levels <- c("strict", "medium", "liberal")
srmr_levels <- c("strict", "medium", "liberal")
nnfi_levels <- c("liberal", "medium", "strict")
cfi_levels <- c("liberal", "medium", "strict")
n_excels_levels <- c("liberal", "medium", "strict")

multiverse_grid <- expand.grid(
  group_cuts = group_cuts,
  sub_cuts = sub_cuts,
  rmsea_cuts = rmsea_cuts,
  srmr_cuts = srmr_cuts,
  nnfi_cuts = nnfi_cuts,
  cfi_cuts = cfi_cuts,
  n_excels = n_excels
) %>% 
  mutate(spec = dplyr::row_number())
# -------------------------------------------------------------------------
# UI ----------------------------------------------------------------------
# -------------------------------------------------------------------------



# Define the UI for the Shiny app
ui <- fluidPage(
  autoWaiter(),

  navbarPage(
  title = "Network Multiverse",
  theme = shinytheme("flatly"),

  
  # Introduction
  tabPanel("Introduction",
           mainPanel(
             tags$div(
               class = "well",
               tags$h2("How to use this app"),
               tags$p("Welcome to the Network Multiverse Shiny App! This interactive tool is a supplement to the analyses in our manuscript on Multiverse analysis in dynamic networks.."),
               tags$p("Follow these steps to use the app:"),
               tags$ol(
                 tags$li("Select a dataset from the dropdown menu below. Choose between 'Personality' and 'Emotion'."),
                 tags$li("Below, you will find an explanation of the variable abbrevations used in this app."),
                 tags$li("Explore the various tabs to analyze the data and visualize insights."),
               ),
               tags$h2("For more information"),
               tags$p("If you need more details about the methodology behind the research, please refer to the corresponding sections in our paper."),
               tags$p("For any inquiries or feedback, feel free to contact us at bjoern.siepe@uni-marburg.de.")

             ),

             selectInput("dataset", 
                         label = "Choose a Dataset", 
                         choices = c("Personality", "Emotion")),
             textOutput("dataset_info"),
             # Display the variable explanation table
             h2("Variable Name Explanations"),
             p("All differences here are calculated as reference fit - multiverse fit. For example, a negative density difference means that the multiverse fit has a higher density."),
             p("Also, note that the lowest subgroup cutoff for the emotion data set is 51% instead of 50%, which is explained in the manuscript."),
             tableOutput("variable_table")
           )
  ),

  
  # Compute summary statistics
  tabPanel("Summary Statistics",
           sidebarLayout(
             sidebarPanel(
               tags$div(
                 class = "well",
                 tags$h3("Grouped Summary Statistics"),
                 tags$p("Grouped summary statistics provide insights into how different variables interact and vary across various groups."),
                 tags$p("Choose one or more variables from the dropdown menu to group the summary statistics by those variables."),
                 tags$p("The string after the underscore(_) at the end of a variable name indicates the summary statistic for the grouping."),
                 tags$p("For example, Heterogeneity_sd is the standard deviation of the heterogeneity within a given grouping."),
                 selectInput("groupby_vars",
                             label = "Group by:",
                             choices = cond_vars,
                             selected = "None",
                             multiple = TRUE)
               )
             ),
             mainPanel(
               tableOutput("fg_summary")
             )
           )
  ),
  
  # Boxplots
  tabPanel("Boxplot",
           sidebarLayout(
             sidebarPanel(
               tags$div(
                 class = "well",
                 tags$h3("Boxplot Visualization"),
                 tags$p("A boxplot is a graphical representation of the distribution of data. It displays the median, quartiles, and potential outliers within a dataset."),
                 tags$p("Select a column from the dropdown menu to generate an interactive boxplot for that specific variable."),
                 selectInput("boxplot_column",
                             label = "Select Column:",
                             choices = summary_vars,
                             selected = summary_vars[1])
               )),
             mainPanel(
               plotlyOutput("boxplot")
             )
           )
  ),
  
  # SCA 
  tabPanel("Specification Curve",
           sidebarLayout(
             sidebarPanel(
               tags$div(
                 class = "well",
                 tags$h3("Specification Curve Analysis"),
                 tags$p("A Specification Curve Analysis (SCA) plot is a visual tool used to analyze the effects of a variable across different specifications."),
                 tags$p("Select a variable from the dropdown menu to explore its impact as ordered by its size across specifications."),
                 tags$p("Rendering of the plot can take some time. This plot is not interactive to speed up computation time."),
                 selectInput("sca_column",
                             label = "Select Column:",
                             choices = summary_vars,
                             selected = summary_vars[1])
               )
             ),
             mainPanel(
               plotOutput("sca_out"),
               plotOutput("sca_spec", height = 650)
             )
           )
  ),
  tabPanel("Network",
           sidebarLayout(
             sidebarPanel(
               tags$div(
                 class = "well",
                 tags$h3("Select Specification Values"),
                 p("Select one value for each specification:"),
                 selectInput("net_groupcutoffs", "groupcutoffs", choices = c(group_cuts), selected = group_cuts[[3]]),
                 selectInput("net_subcutoffs", "subcutoffs", choices = c(sub_cuts), selected = group_cuts[[3]]),
                 selectInput("net_rmsea.cuts", "rmsea.cuts", choices = c(rmsea_cuts), selected = rmsea_cuts[[2]]),
                 selectInput("net_srmr.cuts", "srmr.cuts", choices = c(srmr_cuts), selected = srmr_cuts[[2]]),
                 selectInput("net_nnfi.cuts", "nnfi.cuts", choices = c(nnfi_cuts), selected = nnfi_cuts[[2]]),
                 selectInput("net_cfi.cuts", "cfi.cuts", choices = c(cfi_cuts), selected = cfi_cuts[[2]]),
                 selectInput("net_n.excellent", "n.excellent", choices = c(n_excels), selected = n_excels[[2]]),
               )
             ),
             mainPanel(
               plotOutput("net_plot", width = 500, height = 500),
               tableOutput("net_matrix")
             )
           )
     
  )
  
  
  
  ),
  # Footer-like section
  
  # Footer-like section
  div(
    class = "footer",
    style = "background-color: #f8f9fa; padding: 10px; text-align: center;",
    HTML(paste("Shiny App for the paper \"Network Multiverse\" (Siepe, Schumacher & Heck, 2023). Find the source code on", 
               a(href = "https://github.com/bsiepe/network-multiverse", fa("github")),
               "GitHub."
    ))
  ),
  
  # Custom CSS for the footer-like section
  tags$head(
    tags$style(
      HTML("
        .footer {
          position: fixed;
          bottom: 0;
          width: 100%;
          color: #6c757d; /* Grey font color */
        }
      ")
    )
  )
)

# Define the server for the Shiny app
server <- function(input, output, session) {
  
  #--- Load the selected dataset based on user input
  sel_data <- reactive({
    dataset_name <- input$dataset
    
    if(dataset_name == "Personality") {
      return(readRDS("data/comp_pers_app.RDS"))
    } else if (dataset_name == "Emotion") {
      return(readRDS("data/comp_emot_app.RDS"))
    }
  })
  
  
  output$dataset_info <- renderText({
    dataset_name <- input$dataset
    
    if (dataset_name == "Personality") {
      "You have chosen the personality data set by Wright et al.(2019) containing data of 94 participants.\n
      It can be obtained at https://osf.io/95hyr/"
    } else if (dataset_name == "Emotion") {
      "You have chosen emotion data set by Kullar et al. (2023) containing data of 105 participants.\n
      It can be obtained at https://github.com/mkullar/DataDrivenEmotionDynamics/"
    } else {
      ""
    }
  })
  
  #--- Variables explained
  variable_data <- data.frame(
    Variable = c("Heterogeneity", "Nonzero Edge Diff.", "Adjacency Diff.", "Centrality Diff.",
                 "Density Temp. Abs. Diff.", "Density Cont. Abs. Diff."),
    Explanation = c(
      "Proportion of group-level edges to all edges.",
      "Difference between nonzero edges.",
      "Adjacency matrix difference.",
      "Centrality difference.",
      "Absolute difference of density for temporal network.",
      "Absolute difference of density for contemporaneous network."
    )
  )
  
  # Render the table
  output$variable_table <- renderTable({
    variable_data
  }, rownames = FALSE)

  

  #--- Grouped summary
  # Create reactive expressions for filter and group_by
  grouped_data <- reactive({
    dataset <- sel_data()
    

    group_by_cols <- syms(input$groupby_vars)
    if (!is.null(group_by_cols) & length(group_by_cols) > 0) {
      dataset <- dataset %>%
        group_by(!!!group_by_cols)
    }
    
    return(dataset)
  })
  
  # Compute and display summary statistics
  output$fg_summary <- renderTable({
    grouped_data() %>%
      summarise(across(summary_vars, 
                       .fns = list(mean = mean, sd = sd),
                       .names = "{.col}_{.fn}"))
    
    
  })
  
  
  
  #--- Create ggplot boxplot
  output$boxplot <- renderPlotly({
    selected_column <- input$boxplot_column
    
    # Generate ggplot boxplot
    boxplot_plot <- ggplot(sel_data(), aes(x = "", y = .data[[selected_column]])) +
      geom_boxplot() +
      labs(title = paste("Boxplot of", selected_column))+
      theme_multiverse()
    
    # Convert ggplot plot to a Plotly plot
    ggplotly(boxplot_plot) 
  })
  

  #--- Create SCA
  output$sca_out <- renderPlot({
    selected_column <- input$sca_column
    
    # Generate ggplot boxplot
    plot_outcome(mv_res = sel_data(), var = .data[[selected_column]])+
      theme(axis.title = element_text(size = rel(1.5)),
            axis.text = element_text(size = rel(1.5)))
    
    
    
  })
  
  output$sca_spec <- renderPlot({
    selected_column <- input$sca_column
    
    # Generate ggplot boxplot
    plot_specification(mv_res = sel_data(), var = .data[[selected_column]])+
      theme(axis.title = element_text(size = rel(1.5)),
            axis.text = element_text(size = rel(1.5)))
    
    # # Convert ggplot plot to a Plotly plot
    # ggplotly(sca_spec_plot) 
  })
  
  #--- Create networks
  # Load data based on input 
  sel_net_data <- reactive({
    dataset_name <- input$dataset
    
    if(dataset_name == "Personality") {
      return(readRDS("data/mv_pers_plots.RDS"))
    } else if (dataset_name == "Emotion") {
      return(readRDS("data/mv_emot_plots.RDS"))
    }
  })
  
  
  output$net_plot <- renderPlot({
    spec_ind_plot <- reactive({multiverse_grid %>% 
      filter(group_cuts == input$net_groupcutoffs,
             sub_cuts == input$net_subcutoffs,
             rmsea_cuts == input$net_rmsea.cuts, 
             srmr_cuts == input$net_srmr.cuts,
             nnfi_cuts == input$net_nnfi.cuts,
             cfi_cuts == input$net_cfi.cuts,
             n_excels == input$net_n.excellent) %>% 
        pull(spec)})
      
    plot(sel_net_data()[[spec_ind_plot()]]$group_plot_paths, margin = c(5,5,5,5))  
    
    
  })
  output$net_matrix <- function(){
    spec_ind_mat <- reactive({multiverse_grid %>% 
        filter(group_cuts == input$net_groupcutoffs,
               sub_cuts == input$net_subcutoffs,
               rmsea_cuts == input$net_rmsea.cuts, 
               srmr_cuts == input$net_srmr.cuts,
               nnfi_cuts == input$net_nnfi.cuts,
               cfi_cuts == input$net_cfi.cuts,
               n_excels == input$net_n.excellent) %>% 
        pull(spec)})
    knitr::kable(sel_net_data()[[spec_ind_mat()]]$path_counts) %>% 
      add_header_above(c(" " =  1, "Node of Origin" = dim(sel_net_data()[[spec_ind_mat()]]$path_counts)[2])) %>% 
      kable_styling()
    
  }
  
 
}

# Run the Shiny app
shinyApp(ui, server)




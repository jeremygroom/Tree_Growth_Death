# Forest Analysis Shiny Application
# Load required libraries
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(readxl)
library(here)

# Source required files
source("Global.R")
source(paste0(CODE.LOC, "Functions.R"))

# Load species names data
spp_names_full <- readxl::read_xlsx(paste0(DATA.LOC, "FullSppNames.xlsx"))

# Details for plotting
virid.use <- viridis_pal(option = "H", begin = 0.1, end = 0.9)(n_domain)  # Get colors for plotting

# Create species choice list with readable labels
species_choices <- spp_names_full %>%
  filter(SPCD %in% as.numeric(gsub("X", "", SEL.SPP))) %>%
  mutate(
    species_code = paste0("X", SPCD),
    display_name = paste0(GENUS, " ", SPECIES, " (", COMMON_NAME, ")")
  ) %>%
  arrange(GENUS, SPECIES) %>%
  pull(display_name, name = species_code)

# Load analysis arrays once at startup
load_analysis_data <- function() {
  # Load the main analysis arrays
  if (file.exists(paste0(DATA.LOC, "analysis.arrays.zip"))) {
    temp_file <- unzip(paste0(DATA.LOC, "analysis.arrays.zip"), "analysis.arrays.RDS")
    analysis_data <- readRDS(temp_file)
    file.remove(temp_file)  # Clean up
    return(analysis_data)
  } else {
    stop("analysis.arrays.zip file not found in Data/ directory")
  }
}

# Pre-load the analysis data
analysis_arrays <- load_analysis_data()

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Forest Climate Response Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Analysis", tabName = "analysis", icon = icon("tree"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .box {
          border-radius: 5px;
        }
        .selectize-input {
          font-size: 14px;
        }
      "))
    ),
    
    tabItems(
      tabItem(tabName = "analysis",
              fluidRow(
                box(
                  title = "Analysis Parameters", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 12,
                  
                  fluidRow(
                    column(6,
                           selectInput(
                             "species_selection",
                             "Select Tree Species:",
                             choices = species_choices,
                             selected = names(species_choices)[1],
                             width = "100%"
                           )
                    ),
                    column(6,
                           radioButtons(
                             "analysis_type",
                             "Analysis Type:",
                             choices = list(
                               "Growth Analysis" = "grow",
                               "Mortality Analysis" = "mort"
                             ),
                             selected = "grow",
                             inline = TRUE
                           )
                    )
                  ),
                  
                  fluidRow(
                    column(12,
                           actionButton(
                             "generate_plot",
                             "Generate Analysis Plot",
                             class = "btn-primary btn-lg",
                             style = "width: 200px; margin-top: 10px;"
                           )
                    )
                  )
                )
              ),
              
              fluidRow(
                box(
                  title = "Analysis Results", 
                  status = "success", 
                  solidHeader = TRUE,
                  width = 12,
                  height = "600px",
                  
                  conditionalPanel(
                    condition = "output.plot_ready",
                    div(
                      style = "text-align: center;",
                      plotOutput("analysis_plot", height = "550px")
                    )
                  ),
                  
                  conditionalPanel(
                    condition = "!output.plot_ready",
                    div(
                      style = "text-align: center; margin-top: 200px;",
                      h4("Select species and analysis type, then click 'Generate Analysis Plot'"),
                      icon("chart-line", "fa-3x", style = "color: #3c8dbc; margin-top: 20px;")
                    )
                  )
                )
              ),
              
              fluidRow(
                box(
                  title = "Analysis Information",
                  status = "info",
                  solidHeader = TRUE,
                  width = 12,
                  
                  h4("About This Analysis"),
                  p("This application analyzes tree growth and mortality responses to climate variables across different environmental domains."),
                  
                  h5("Climate Variable:"),
                  p(paste("Current analysis uses:", CLIM.VAR.USE, "data")),
                  
                  h5("Domains:"),
                  p(paste("Analysis includes", n_domain, "environmental domains:", paste(DOMAIN.LEVELS, collapse = ", "))),
                  
                  h5("Plot Components:"),
                  tags$ul(
                    tags$li("Left panel: Response rates by environmental domain"),
                    tags$li("Center panel: Distribution of plots across climate conditions"),
                    tags$li("Right panel: Geographic distribution of sample plots")
                  )
                )
              )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive value to track plot readiness
  plot_ready <- reactiveVal(FALSE)
  
  # Load domain analysis output based on analysis type
  load_domain_output <- reactive({
    k <- ifelse(input$analysis_type == "grow", 1, 2)
    output_file <- paste0(save.loc.fcn(k), "Domain_Analysis_Output.RDS")
    
    if (file.exists(output_file)) {
      return(readRDS(output_file))
    } else {
      return(NULL)
    }
  })
  
  # Generate plot when button is clicked
  observeEvent(input$generate_plot, {
    
    # Show loading message
    showNotification("Generating plot...", type = "message", duration = 2)
    
    # Set plot as not ready while generating
    plot_ready(FALSE)
    
    # Generate the plot
    tryCatch({
      
      # Determine k value (1 = growth, 2 = mortality)
      k <- ifelse(input$analysis_type == "grow", 1, 2)
      
      # Load the domain analysis output
      plt_dat <- load_domain_output()
      
      if (is.null(plt_dat)) {
        showNotification(
          paste("Domain analysis output not found for", input$analysis_type), 
          type = "error", 
          duration = 5
        )
        return()
      }
      
      plt_dat2 <- plt_dat$domain.summaries
      
      # Get the appropriate analysis arrays
      mort_grow_dat <- if (k == 1) analysis_arrays$arrays.grow else analysis_arrays$arrays.mort
      # Need climate names for files and axes.
      clim.names <- read_csv(paste0(DATA.LOC, "ClimateNames.csv"), show_col_types = FALSE) %>%
        filter(filename == CLIM.VAR.USE)
      
      var1 <<- clim.names$values[grep("pre", clim.names$values)]
      var.delt <<- clim.names$values[grep("d", clim.names$values)]
      
      
      var.label <<- clim.names$label[clim.names$values == var1]
      #var.filename <- clim.names$filename[clim.names$values == var1]
      var.delt.label <<- clim.names$label[clim.names$values == var.delt]
      
      
      domain_matrix <- mort_grow_dat$domain.matrix
      domain_n <- mort_grow_dat$domain.n 
      quant_lims <- mort_grow_dat$quant.lims
      quant.lims.delta <<- mort_grow_dat$quant.lims.delt
      #browser()
      # Generate the plot using your existing function
      pair.plts.fcn(
        sppnum.to.plot = names(which(species_choices == input$species_selection)),
        use.dat = plt_dat2,
        domain.matrix = domain_matrix,
        quant.lims = quant_lims,
        domain.n = domain_n,
        k = k,
        SHINYAPP.IN.USE = TRUE
      )
      
      # Set plot as ready
      plot_ready(TRUE)
      
      showNotification("Plot generated successfully!", type = "message", duration = 3)
      
      
    }, error = function(e) {
      showNotification(
        paste("Error generating plot:", e$message), 
        type = "error", 
        duration = 10
      )
      plot_ready(FALSE)
    })
  })
  
  # Render the plot
  output$analysis_plot <- renderPlot({
    if (plot_ready()) {
      k <- ifelse(input$analysis_type == "grow", 1, 2)
      plot_file <- paste0(save.loc.fcn(k), input$species_selection, "_plots.png")
      
      if (file.exists(plot_file)) {
        # Read and display the saved plot
        img <- png::readPNG(plot_file)
        grid::grid.raster(img)
      }
    }
  })
  
  # Output for conditional panel
  output$plot_ready <- reactive({
    plot_ready()
  })
  outputOptions(output, "plot_ready", suspendWhenHidden = FALSE)
  
  # Display selected species info
  observe({
    if (!is.null(input$species_selection)) {
      species_info <- spp_names_full %>%
        filter(paste0("X", SPCD) == input$species_selection)
      
      if (nrow(species_info) > 0) {
        updateSelectInput(
          session,
          "species_selection",
          label = paste0("Selected: ", species_info$GENUS, " ", species_info$SPECIES, " (", species_info$COMMON_NAME, ")")
        )
      }
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
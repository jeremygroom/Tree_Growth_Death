# Tree Growth and Mortality Analysis Shiny App

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)

# Source required files
source("Global.R")
source(paste0(CODE.LOC, "Functions.R"))

# UI
ui <- fluidPage(
  tags$head(  # Container for logo
    tags$style(HTML("
      .logo-container {
        position: absolute;
        top: 10px;
        right: 20px;
        z-index: 1000;
      }
      .logo-container img {
        transition: opacity 0.3s ease;
      }
      .logo-container img:hover {
        opacity: 0.8;
      }
    "))
  ),
  
  titlePanel("FIA Tree Growth and Mortality Analysis"),
  
  # Logo
  tags$div(
    class = "logo-container",
    tags$a(
      href = 'http://www.groomanalytics.com',
      target = "_blank",
      tags$img(src = 'GroomAnalyticsH.jpg', 
               height = "40px",
               alt = "Groom Analytics",
               title = "Groom Analytics Home")
    )
  ),
  
  # Control panel at top
  fluidRow(
    column(3,
           selectInput("species", 
                       "Select Species:",
                       choices = NULL,  # Will be populated in server
                       selected = NULL)
           ),    
    column(3,
           selectInput("analysis_type", 
                       "Analysis Type:",
                       choices = list("Growth" = 1, "Mortality" = 2),
                       selected = 1)
    )
  ),
  

  
  # Main plot area
  fluidRow(
    column(12,
           plotOutput("main_plot", height = "600px")
    )
  ),
  
  # Status text
  fluidRow(
    column(12,
           textOutput("status")
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Load data at startup
  data_loaded <- reactiveVal(FALSE)
  analysis_data <- reactiveValues()
  species_choices <- reactiveVal(NULL)
  
  # Load climate variable names (used to define var1 etc.)
  clim.names <- read_csv(paste0(DATA.LOC, "ClimateNames.csv"), show_col_types = FALSE) %>%
    filter(filename == CLIM.VAR.USE)
  
  var1 <<- clim.names$values[grep("pre", clim.names$values)]
  var.delt <<- clim.names$values[grep("d", clim.names$values)]
  
  var.label <<- clim.names$label[clim.names$values == var1]
  var.delt.label <<- clim.names$label[clim.names$values == var.delt]
  
  
  
  # Load all data on startup
  observe({
    withProgress(message = 'Loading data...', value = 0, {
      
      
      
      # Load species names and create choices
      incProgress(0.2, detail = "Loading species information...")
      spp_names_data <- readxl::read_xlsx(paste0(DATA.LOC, "FullSppNames.xlsx"))
      
      # Create species choices in the format "Genus species - Common name"
      choices_df <- spp_names_data %>%
        filter(paste0("X", SPCD) %in% SEL.SPP) %>%
        mutate(
          label = paste0(GENUS, " ", SPECIES, " - ", COMMON_NAME),
          value = paste0("X", SPCD)
        ) %>%
        arrange(GENUS, SPECIES)
      
      choices_list <- setNames(choices_df$value, choices_df$label)
      species_choices(choices_list)
      
      # Load analysis arrays
      incProgress(0.4, detail = "Loading analysis arrays...")
      if (file.exists(paste0(DATA.LOC, "analysis.arrays.zip"))) {
        arrays_data <- readRDS(unzip(paste0(DATA.LOC, "analysis.arrays.zip"), 
                                     "analysis.arrays.RDS"))
        file.remove("analysis.arrays.RDS")  # Clean up unzipped file
        
        analysis_data$arrays_mort <- arrays_data$arrays.mort
        analysis_data$arrays_grow <- arrays_data$arrays.grow
      } else {
        showNotification("Analysis arrays file not found", type = "error")
        return()
      }
      
      # Load domain analysis outputs
      incProgress(0.6, detail = "Loading growth analysis...")
      growth_file <- paste0(save.loc.fcn(1), "Domain_Analysis_Output.RDS")
      if (file.exists(growth_file)) {
        analysis_data$growth_output <- readRDS(growth_file)
      } else {
        showNotification("Growth analysis file not found", type = "error")
        return()
      }
      
      incProgress(0.8, detail = "Loading mortality analysis...")
      mort_file <- paste0(save.loc.fcn(2), "Domain_Analysis_Output.RDS")
      if (file.exists(mort_file)) {
        analysis_data$mort_output <- readRDS(mort_file)
      } else {
        showNotification("Mortality analysis file not found", type = "error")
        return()
      }
      
      incProgress(1.0, detail = "Data loading complete!")
      data_loaded(TRUE)
    })
  })
  
  # Update species choices when data is loaded
  observe({
    if (data_loaded() && !is.null(species_choices())) {
      updateSelectInput(session, "species", 
                        choices = species_choices(),
                        selected = species_choices()[[1]])
    }
  })
  
  # Generate plot
  output$main_plot <- renderPlot({
    if (!data_loaded() || is.null(input$species) || input$species == "") {
      return(NULL)
    }

    tryCatch({
      # Get analysis type
      k <- as.numeric(input$analysis_type)
      
      # Get the appropriate data
      if (k == 1) {
        # Growth
        mort_grow_dat <- analysis_data$arrays_grow
        plt_dat <- analysis_data$growth_output
      } else {
        # Mortality
        mort_grow_dat <- analysis_data$arrays_mort
        plt_dat <- analysis_data$mort_output
      }
      

      plt_dat2 <- cm2.fcn(k, plt_dat$domain.summaries) # Transforming growth data to cm2
      domain.matrix <- mort_grow_dat$domain.matrix
      domain.n <- mort_grow_dat$domain.n
      quant.lims <- mort_grow_dat$quant.lims
      quant.lims.delta <<- mort_grow_dat$quant.lims.delt
      
      virid.use <<- viridis_pal(option = "H", begin = 0.1, end = 0.9)(n_domain)  # Get colors for plotting
      
      
#     browser()
      # Generate the plot using pair.plts.fcn with SHINYAPP.IN.USE = TRUE
      plot_result <- pair.plts.fcn(
        sppnum.to.plot = input$species, 
        use.dat = plt_dat2, 
        domain.matrix = domain.matrix,
        quant.lims = quant.lims, 
        domain.n = domain.n, 
        k = k, 
        SHINYAPP.IN.USE = TRUE,
        SHINY.FONTSIZE = SHINY.FONTSIZE
      )
      
      return(plot_result)
      
    }, error = function(e) {
      # Create a simple error plot
      plot(1, 1, type = "n", xlab = "", ylab = "", main = "Error occurred", 
           axes = FALSE, frame.plot = FALSE)
      text(1, 1, "There was an error in the code or there was insufficient data for selected species and analysis type", 
           cex = 1.2, col = "red")
    })
  })
  
  # Status output
  output$status <- renderText({
    if (!data_loaded()) {
      "Loading data..."
    } else if (is.null(input$species) || input$species == "") {
      "Please select a species"
    } else {
      analysis_type_name <- ifelse(input$analysis_type == 1, "Growth", "Mortality")
      species_name <- names(species_choices())[species_choices() == input$species]
      paste("Displaying", analysis_type_name, "analysis for", species_name)
    }
  })
}

# Run the app
shinyApp(ui = ui, server = server)
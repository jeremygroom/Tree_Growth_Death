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
  
  titlePanel("Tree Growth and Mortality by Climatic Water Deficit Domain"),
  
  
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
  ),
  fluidRow(column(width = 10, offset = 0,
                  div(htmlOutput("box9"),  style = 'color:black; border-radius: .5em; font-size:20px; padding: 20px; margin: 20px;' , # background-color:#b8d9b9;
                      h2("Dashboard Synopsis"),
                      h3("Overview"),
                      p("This page offers a data visualization tool for the manuscript 'Tree growth, mortality, and climatic water deficit in west-coast states, USA' by JD Groom (", 
                        a("Groom Analytics LLC", href = "https://www.groomanalytics.com/"), ") and B Frank (USDA Forest Service; 2026).
                         We examined the average growth and mortality rates for 23 tree species across California, Oregon, and Washington. The tree data were obtained from study plots measured by the ", a("USDA Forest 
                        Service Forest Inventory and Analysis Program", href = "https://research.fs.usda.gov/programs/fia"), ", or FIA. For these figures we look at tree growth and mortality according to the water stress, measured as Climatic 
                        Water Deficit (CWD), that each FIA sample plot experienced. Go ahead and select a species at the top along with an analysis type (Growth or Mortality)."),
                      h3("Figure: FIA plot distributions by CWD"), 
                      p("Since we can only directly measure tree growth and death rates by visiting individual trees over time, it is possible that the water stress a tree experiences can change over time as well. 
                       Some places get drier while others receive more rain or snow. Our data comes from FIA plots that were visited twice, 10 years apart. The x-axis of the middle figure, \"FIA plot distribution by CWD\", 
                       captures CWD value at the initial study plot visit.  The figure's y-axis depicts the change in CWD between visits. 
                         For all species we define a plot as having an Increasing (I) CWD value if the y-axis value is above 5 mm, and a Stable (S) CWD value if it is below 5 mm. The horizontal blue dashed line is
                         at 5 mm of CWD change between the two plot visits.  For the x-axis, the initial CWD values are broken into species-specific quartiles, defined by the vertical orange lines.  The 
                         25% of data on the left are Low (L) initial CWD values, the middle 50% of points are 
                        Medium (M), and the upper 25% are High (H). We thus have six domains, IL (Increasing, Low), IM (Increasing, Medium), IH (Increasing, High), SL (Stable, Low), 
                        SM (Stable, Medium), and SH (Stable, High). Domain SL (yellow) includes plots with little increase in CWD at the lower range of initial CWD values. These plots arguably experience the least water stress. 
                        IH (light green) is the opposite; these FIA plots have the greatest water stress. This figure shows, for each species, the relationship and distribution of plot values across the six domains."),
                      h3("Figure: CWD domain estimates"),
                        p("The leftmost figure gives us the estimated 10-year growth or mortality rates (depending on which you selected) for the trees in each of the domains.  That is, the dark blue column for IL 
                          is the mean growth or mortality estimate for the trees in the FIA plots represented by the dark-blue points in the middle figure \"FIA plot distrubution by CWD\". We use the same 
                          domain abbreviations (IL, IM, etc.) as before.  The error bars are the 95% bootstrapped confidence intervals"),
                      h3("Figure: FIA fuzzed plot locations"),
                      p("The right-most figure provides a map of the approximate geographic location of each of the FIA plots in which your selected species was found.  The color scheme of the points remains 
                        consistent with the domain labels at the bottom of the middle figure.  We can see, for instance, that ", em("Abies amablis"), " SL plots are mostly found in northern Washington in the Cascades range
                         and in the upper Olympic Peninsula.")
                      
                  )
  ))
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
      
      
    # browser()
      # Generate the plot using pair.plts.fcn with SHINYAPP.IN.USE = TRUE
      plot_result <- pair.plts.fcn(
        sppnum.to.plot = input$species, 
        use.dat = plt_dat2, 
        domain.matrix = domain.matrix,
        quant.lims = quant.lims, 
        domain.n = domain.n, 
        k = k, 
        SHINYAPP.IN.USE = TRUE,
        SHINY_FONTSIZE = SHINY.FONTSIZE
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
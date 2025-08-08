################################################################################

# PWS Herring Mid-year Management Tool

# Author: CL Roberts

################################################################################


#------------------------------ front matter ----------------------------------#

# attach packages ---
library(shiny)
library(bslib)
library(gt)
library(zip)

# source utility scripts ----
source("helper.r")
source("setup.r")

#------------------------------ user interface --------------------------------#

ui  <- page_sidebar(

  tags$head(
    tags$style(
      HTML("
        .centered-card {
          display: flex; align-items: center; 
        }
        .shiny-output-error-warning {
          display: flex; align-items: center; color: red;
        }
        body {
          background-image: url('SoA-background.jpg');
          /* background-image: linear-gradient(to right, rgba(255,255,255, 0.7) 0 100%), url('SoA-background.jpg'); */
          background-size: cover;
          background-repeat: no-repeat;
          background-attachment: fixed;
        }
      ")
    ),
    tags$script(
      HTML('
        $(document).ready(function() {
          $(".navbar .container-fluid")
            .append("<a href=https://github.com/cl-roberts/pws-herring-basa target=_blank> <img id=flag src=AK-flag.svg align=right height=45px> </a>");
        });
      ')
    )
  ),


  # app theme ----
  theme = bs_theme(
      bootswatch = "sandstone", 
      font_scale = .8, 
      "bslib-page-sidebar-title-color" = theme_offwhite,
      "bslib-page-sidebar-title-bg" = theme_navy,
      "bslib-sidebar-bg" = theme_navy,
      "bslib-border-top-color" = theme_gold,
      "bslib-sidebar-button-collapse-toggle" = theme_gold
      ) |>
      bs_add_rules(
        "
        .navbar {
          border-top: 3px solid #edbd03; 
        }
        "),

  # app title ----
  title = h2(paste("PWS Herring", curr_year,"Mid-year Management Tool")),

  # sidebar panel for inputs ----
  sidebar = sidebar(title = NULL, width = 350,

    # required fields ----
    card(card_header(h4("Required Fields")),

      # NAA input ----
      h6(paste("Spring", curr_year, "catch-at-age (millions)")),

      fluidRow( 

        splitLayout(
          numericInput("caaAge3", "Age 3:", value = NA, min = 0),
          numericInput("caaAge4", "Age 4:", value = NA, min = 0),
          numericInput("caaAge5", "Age 5:", value = NA, min = 0),
          numericInput("caaAge6", "Age 6:", value = NA, min = 0)
        ),

        splitLayout(
          numericInput("caaAge7", "Age 7:", value = NA, min = 0),
          numericInput("caaAge8", "Age 8:", value = NA, min = 0),
          numericInput("caaAge9", "Age 9+:", value = NA, min = 0)
        )    

      ),

      # WAA input ----
      h6(paste("Spring", curr_year, "weight-at-age (g)")),

      fluidRow( 

        splitLayout(
          numericInput("waaAge3", "Age 3:", value = NA, min = 0),
          numericInput("waaAge4", "Age 4:", value = NA, min = 0),
          numericInput("waaAge5", "Age 5:", value = NA, min = 0),
          numericInput("waaAge6", "Age 6:", value = NA, min = 0)
        ),

        splitLayout(
          numericInput("waaAge7", "Age 7:", value = NA, min = 0),
          numericInput("waaAge8", "Age 8:", value = NA, min = 0),
          numericInput("waaAge9", "Age 9+:", value = NA, min = 0)
        )    

      ),

      # female proportion and MDM input ----
      fluidRow( 

        splitLayout(
          h6(HTML("Female proportion")),
          h6(HTML("Observed MDM")),
        ),

        splitLayout(
          numericInput("propFemale", NULL, value = 0.5, min = 0.001, max = 0.999),
          numericInput("mdmObserved", NULL, value = NA, min = 0)
        )    

      )

    ),

    card(card_header(h4("Optional Fields")),

      h6(paste("Expected fall", curr_year+1, "food/bait catch (millions)")),

      # expected fall catch-at-age
      fluidRow(

        splitLayout(
          numericInput("fallCaaAge2", "Age 2:", value = NA, min = 0),
          numericInput("fallCaaAge3", "Age 3:", value = NA, min = 0),
          numericInput("fallCaaAge4", "Age 4:", value = NA, min = 0),
          numericInput("fallCaaAge5", "Age 5:", value = NA, min = 0)
        ),

        splitLayout(
          numericInput("fallCaaAge6", "Age 6:", value = NA, min = 0),
          numericInput("fallCaaAge7", "Age 7:", value = NA, min = 0),
          numericInput("fallCaaAge8", "Age 8:", value = NA, min = 0),
          numericInput("fallCaaAge9", "Age 9+:", value = NA, min = 0)
        )
        
      ),

      # Observed current year NAA agecomp ----
      h6(paste("Spring", curr_year, "numbers-at-age composition (proportion)")),

      fluidRow( 

        splitLayout(
          numericInput("naaAgecomp3", "Age 3:", value = NA, min = 0, max = 1),
          numericInput("naaAgecomp4", "Age 4:", value = NA, min = 0, max = 1),
          numericInput("naaAgecomp5", "Age 5:", value = NA, min = 0, max = 1),
          numericInput("naaAgecomp6", "Age 6:", value = NA, min = 0, max = 1)
        ),

        splitLayout(
          numericInput("naaAgecomp7", "Age 7:", value = NA, min = 0, max = 1),
          numericInput("naaAgecomp8", "Age 8:", value = NA, min = 0, max = 1),
          numericInput("naaAgecomp9", "Age 9+:", value = NA, min = 0, max = 1)
        )    

      ),

    ),

    # action and download buttons ----
    actionButton("doCalculations", "Calculate", 
                  icon = icon("circle-right")),
    downloadButton("downloadTbl", "Download Estimates"),
    downloadButton("downloadFigs", "Download Figures")
  
  ),

  navset_card_tab(

    # one-year biomass forecast tab ----
    nav_panel(paste(curr_year, "Biomass Forecast"),

      h2("Welcome"),
      p("This application is intended to assist in the mid-year management of Prince William Sound herring."),

      fluidRow(
        column(width = 4,

          div(
            card(
              class = "centered-card",
              card_header(h4(paste("Spring", curr_year, "Forecast"))),
              h3(Btilde_forecast_med)
            )      
          ),

          div(
            card(
              class = "centered-card",
              card_header(h4("95% Credible Interval")), 
              h3(Btilde_forecast_ci)
            )
          ),

          div(
            card(
              class = "centered-card",
              card_header(h4("Probability below threshold")), 
              h3(prob_below_threshold)
            )
          )

        ),

        column(width = 8, plotOutput("biomassForecastPlot"))

      )
    ),

    # agecomp forecast tab ----
    nav_panel(paste(curr_year, "Agecomp Forecast"), 

        h2(paste(curr_year, "age composition forecast")),

        layout_columns(

          gt(filter(agecomp_tbl, Type == "Numbers")) |>
            fmt_number(decimals = 2) |>
            tab_header(
              title = md("*Numbers-at-age Composition*")
            ) |>
            cols_hide(Type),

          gt(filter(agecomp_tbl, Type == "Biomass")) |>
            fmt_number(decimals = 2) |>
            tab_header(
              title = md("*Biomass-at-age Composition*")
            ) |>
            cols_hide(Type),

          plotOutput("naaAgecompPlot"),
          plotOutput("baaAgecompPlot"),

          col_widths = c(6, 6, 6, 6),
          row_heights = c(1, 1)

        )

    ),

    # MDM forecast tab ----
    nav_panel(paste(curr_year, "MDM Forecast"), 
    
      h2(paste(curr_year, "mile-days milt forecast")),
      p("A low probability indicates that the observed value mile-days milt was substantially less than what the bayesian age-structure assessment (BASA) model expected, and suggests that the spring 2025 mature biomass estimate may have been overforcasted. This could be due to a smaller than expected recruitment, for example."),

      fluidRow(

        column(width = 4,

          div(
            card(
              class = "centered-card",
              card_header(h4(paste("Spring", curr_year, "Forecast"))),
              h5(textOutput("mdmEstimate"))
            )      
          ),

          div(
            card(
              class = "centered-card",
              card_header(h4("95% Credible Interval")), 
              h5(textOutput("mdmCI"))
            )
          ),

          div(
            card(
              class = "centered-card",
              card_header(h4("Probability below observed")), 
              h5(textOutput("mdmProb"))
            )
          )

        ), 

        column(width = 8,
          plotOutput("mdmPlot")
        )

      )
    
    ),

    # two-year biomass forecast tab ----
    nav_panel(paste(curr_year+1, "Biomass Forecast"), 

      h2(paste(curr_year+1, "biomass forecast")),
      p("This is a two-year forecast. Caution is advised."),

      fluidRow(
        column(width = 4,

          div(
            card(
              class = "centered-card",
              card_header(h4(paste("Spring", curr_year+1, "Forecast"))),
              h3(textOutput("biomassTwoyearEstimate"))
            )      
          ),

          div(
            card(
              class = "centered-card",
              card_header(h4("95% Credible Interval")), 
              h3(textOutput("biomassTwoyearCI"))
            )
          ),

          div(
            card(
              class = "centered-card",
              card_header(h4("Probability below threshold")), 
              h3(textOutput("biomassTwoyearProb"))
            )
          )

        ),

        column(width = 8, plotOutput("biomassPlotTwoyear"))

      )    

    ),

    # debugging panel ----
    # nav_panel("Debug", 
    #   card(
    #     card_header(h4("Debugging pane")), 
    #     textOutput('debug')
    #   )
    # )

  )

)

#--------------------------------- server -------------------------------------#

server <- function(input, output) {

  # set up reactive table for downloading
  # boolean values keep track of whether plots are rendered or not
  rv <- reactiveValues(
    tbl = forecastTbl,
    naaAgecompPlotShown = TRUE,
    baaAgecompPlotShown = TRUE,
    mdmPlotShown = FALSE,
    biomassPlotTwoyearShown = FALSE
  )

  # static one-year biomass forecast plot ----
  output$biomassForecastPlot <- renderPlot({
      return(biomass_forecast_plot)
  })

  # reactive current year WAA ----
  currYearWAA <- eventReactive(input$doCalculations, {

    currYearWAA <- c(
      input$waaAge3, input$waaAge4, input$waaAge5, input$waaAge6,
      input$waaAge7, input$waaAge8, input$waaAge9
    )

    return(currYearWAA)

  })

  # reactive current year spring catch ----
  springCatch <- eventReactive(input$doCalculations, {

    springCatch <- cbind(
        input$caaAge3, input$caaAge4, input$caaAge5, input$caaAge6,
        input$caaAge7, input$caaAge8, input$caaAge9
    ) 
    springCatch <- ifelse(is.na(springCatch), 0, springCatch)

    return(springCatch)

  })

  # reactive current year spring catch ----
  expectedFallCatch <- eventReactive(input$doCalculations, {

    expectedFallCatch <- cbind(
        input$fallCaaAge2, input$fallCaaAge3, input$fallCaaAge4, input$fallCaaAge5, 
        input$fallCaaAge6, input$fallCaaAge7, input$fallCaaAge8, input$fallCaaAge9
    ) 
    expectedFallCatch <- ifelse(is.na(expectedFallCatch), 0, expectedFallCatch)

    return(expectedFallCatch)

  })

  # collect observed NAA age comps ----
 naaAgecompInput <- eventReactive(input$doCalculations, {

    naaAgecompInput <- c(
        input$naaAgecomp3, input$naaAgecomp4, input$naaAgecomp5, 
        input$naaAgecomp6, input$naaAgecomp7, input$naaAgecomp8, 
        input$naaAgecomp9
    )

    return(naaAgecompInput)

  })

  # make observed NAA age comps df ----
  naaAgecompObserved <- reactive({

    naaAgecompObserved <- data.frame(
        Type = "Numbers", 
        Age = paste("Age", c(3:8, "9+")),
        Proportion = naaAgecompInput()
    )

    return(naaAgecompObserved)

  })

  # calculate biomass age comps ----
  baaAgecompObserved <- reactive({

    agecompsCurrYearWAA <- naaAgecompObserved()$Proportion*currYearWAA()
    baaAgecompObserved <- data.frame(
        Type = "Biomass", 
        Age = paste("Age", c(3:8, "9+")),
        Proportion = agecompsCurrYearWAA / sum(agecompsCurrYearWAA)
    )

    return(baaAgecompObserved)

  })

  # make numbers-at-age comp plot ----
  naaAgecompPlot <- reactive({

    naaAgecompPlot <- naa_agecomp_fig

    if (input$doCalculations > 0) {  
      if (!all(is.na(naaAgecompObserved()$Proportion))) {

        naaAgecompPlot <- naaAgecompPlot +
          geom_point(data = naaAgecompObserved(),
            aes(x = Age, y = Proportion, shape = "Observed"), 
                color = theme_gold, size = 4, na.rm = TRUE) +
          scale_shape_manual(values = c("Observed" = 3)) +
          labs(shape = NULL) +
          theme(legend.position.inside = c(.8, .8), 
                legend.position = "inside",
                legend.background = element_rect(fill = "transparent"))
        
      }
    }

    return(naaAgecompPlot)

  })

  # make biomass age comp plot ----
  baaAgecompPlot <- reactive({

    baaAgecompPlot <- baa_agecomp_fig

    if (input$doCalculations > 0) {  
      if (!all(is.na(baaAgecompObserved()$Proportion))) {

        baaAgecompPlot <- baaAgecompPlot +
          geom_point(data = baaAgecompObserved(),
            aes(x = Age, y = Proportion, shape = "Observed"), 
                color = theme_gold, size = 4, na.rm = TRUE) +
          scale_shape_manual(values = c("Observed" = 3)) +
          labs(shape = NULL) +
          theme(legend.position.inside = c(.8, .8), 
                legend.position = "inside",
                legend.background = element_rect(fill = "transparent"))

      }
    }

    return(baaAgecompPlot)

  })


  # add observed numbers-at-age comp to agecomp plots ----
  output$naaAgecompPlot <- renderPlot({

    # check that agecomps sum to 1
    if (input$doCalculations > 0) {  
      naaAgecompSum <- sum(na.omit(naaAgecompObserved()$Proportion))
      validate(errorClass = "warning",
        need((abs(naaAgecompSum-1) < .01) | all(is.na(naaAgecompInput())), 
            message = paste0("Numbers-at-age composition must sum to 1, not ", naaAgecompSum, "!"))
      )
    }

    naaAgecompPlot()
  
  })

  # add observed biomass age comp to agecomp plots ----
  output$baaAgecompPlot <- renderPlot({

    # check that observed age classes with nonzero proportions have observed weights
    if (input$doCalculations > 0) {  

      naaAgeCompZeroToNA <- ifelse(naaAgecompInput() == 0, NA, naaAgecompInput())
      currYearWAAZeroToNA <- ifelse(currYearWAA() == 0, NA, currYearWAA())

      if (!all(is.na(naaAgeCompZeroToNA)) & !all(is.na(currYearWAAZeroToNA))) {
   
        naaAgeCompHaveObservations <- which(is.na(naaAgeCompZeroToNA))
        currYearWAAHaveObservations <- which(is.na(currYearWAAZeroToNA))
   
        validate(errorClass = "warning",
          need(all(naaAgeCompHaveObservations %in% currYearWAAHaveObservations), 
              message = paste0("Some age classes with observed weights don't have nonzero agecomp proportions!"))
        )

        validate(errorClass = "warning",
          need(all(currYearWAAHaveObservations %in% naaAgeCompHaveObservations), 
              message = paste0("Some observed age classes don't have weights!"))
        )

      }
    }

    baaAgecompPlot()
  })

  # calculate expected post-fishery MDM ----
  mdmPosterior <- eventReactive(input$doCalculations, {

    if (any(is.na(currYearWAA()))) {
      springYield <- apply(maturity, MARGIN = 1, FUN = \(x) sum(x*springCatch()*waa_forecast))
    } else {
      springYield <- apply(maturity, MARGIN = 1, FUN = \(x) sum(x*springCatch()*currYearWAA()))
    }
    
    postfisheryBiomass <- Btilde_forecast - springYield
    mdmPosterior <- postfisheryBiomass * (1-input$propFemale) / exp(logmdm_c)
    
    return(mdmPosterior)
  })

  # render MDM summary ----
  output$mdmEstimate <- renderText({ 
    round(median(mdmPosterior()), 1) 
  })

  output$mdmCI <- renderText({ 
    paste(round(quantile(mdmPosterior(), c(.025, .975)), 1), collapse = ", ") 
  })

  output$mdmProb <- renderText({ 
    paste0(round(100*sum(mdmPosterior() < input$mdmObserved) / n_iters, 1), "%") 
  })

  # make MDM histogram ----
  mdmPlot <- reactive({

    mdmPlotData <- data.frame(mdm = mdmPosterior()) 
    
    mdmPlot <- ggplot(mdmPlotData) +
      geom_histogram(aes(x = mdm, y = after_stat(count / sum(count))), 
                     fill = "NA", color = "black", bins = 50) +
      xlab("Mile-days Milt") +
      ylab("Frequency") +
      labs(title = "MDM Posterior Distribution", color = NULL) 

    if (!is.na(input$mdmObserved)) {
      mdmPlot <- mdmPlot +
        geom_vline(aes(xintercept = input$mdmObserved, color = "Observed MDM")) +
        scale_color_manual(values = theme_gold) +
        theme(legend.position.inside = c(.8, .8), legend.position = "inside")
    }

    return(mdmPlot)

  })

  # render MDM plot ----
  output$mdmPlot <- renderPlot({
    mdmPlot()
  })

  # calculate two-year biomass forecast ----
  biomassForecastTwoyear <- eventReactive(input$doCalculations, {

    postfisheryNAA <- apply(N_a_forecast, MARGIN = 1, FUN = \(x) x - springCatch()) |>
      t()
    fallNAA <- postfisheryNAA*survival_forecast

    twoyearNAA <- cbind(
      (exp(mean_log_rec)/survival_forecast_age2 - expectedFallCatch()[1])*survival_forecast_age2, 
      (fallNAA[,1:5] - expectedFallCatch()[2:6]) * survival_forecast[,1:5],
      (fallNAA[,6]-expectedFallCatch()[7])*survival_forecast[,6] + (fallNAA[,7]-expectedFallCatch()[8])*survival_forecast[,7]
    )

    if (any(is.na(currYearWAA()))) {
      twoyearWAAForecast <- waa_forecast
    } else {
      twoyearWAAForecast <- window_average(rbind(waa, currYearWAA()), waa_average_years)
    }

    biomassForecastTwoyear <- apply(twoyearNAA*maturity, MARGIN = 1, FUN = \(x) x %*% twoyearWAAForecast)  

    return(biomassForecastTwoyear)

  })

  # convert two-year biomass forecast to tons ----
  biomassForecastTwoyearTons <- eventReactive(input$doCalculations, {

    biomassForecastTwoyearTons <- biomassForecastTwoyear()*1.10231

    return(biomassForecastTwoyearTons)

  })

  # render two-year forecast summary ----
  output$biomassTwoyearEstimate <- renderText({ 
    median(biomassForecastTwoyearTons()) |>   
      round() |>
      prettyNum(big.mark = ",")
  })

  output$biomassTwoyearCI <- renderText({ 
    quantile(biomassForecastTwoyearTons(), c(.025, .975)) |>
      round() |>
      prettyNum(big.mark = ",") |>
      paste(collapse = ", ")
  })

  output$biomassTwoyearProb <- renderText({ 
    paste0(
      round(100*sum(biomassForecastTwoyearTons() < threshold) / n_iters, 1), 
      "%"
    ) 
  })

  # make two-year biomass forecast histogram ----
  biomassPlotTwoyear <- reactive({
    
    biomassPlotTwoyearData <- data.frame(biomass = biomassForecastTwoyearTons()) 
    
    if (input$doCalculations > 0) {
      biomassPlotTwoyear <- ggplot(biomassPlotTwoyearData) +
        geom_histogram(aes(x = biomass, y = after_stat(count / sum(count))), 
                      fill = "NA", color = "black", bins = 50) +
        xlab("Mature Biomass (tons)") +
        ylab("Frequency") +
        labs(title = paste(curr_year+1, "Biomass Forecast Posterior Distribution"), 
            color = NULL) +
        geom_vline(aes(xintercept = threshold, color = "22,000-ton Threshold")) +
        scale_color_manual(values = theme_gold) +
        theme(legend.position.inside = c(.8, .8), legend.position = "inside")
    }

    return(biomassPlotTwoyear)

  })

  # render two-year biomass plot ----
  output$biomassPlotTwoyear <- renderPlot({ 
    biomassPlotTwoyear()
  })

  # add rows to download table upon calculations ----
  observeEvent(input$doCalculations, {

    # remove previously added rows
    if (input$doCalculations) {
        rv$tbl <- forecastTbl
    }

    # add row for MDM forecast
    rv$tbl <- rv$tbl |> 
      add_row(
        Quantity = "Mile-days Milt Forecast",
        Year = curr_year, 
        Age = "3+", 
        Units = "Mile-days milt", 
        `50%` = median(mdmPosterior()), 
        `2.5%` = quantile(mdmPosterior(), .025),
        `97.5%` = quantile(mdmPosterior(), .975)
      )

    # add row for 2026 biomass forecast
    rv$tbl <- rv$tbl |> 
      add_row(
        Quantity = "Mature Biomass Forecast",
        Year = curr_year+1, 
        Age = "3+", 
        Units = "Tons", 
        `50%` = median(biomassForecastTwoyear()), 
        `2.5%` = quantile(biomassForecastTwoyear(), .025),
        `97.5%` = quantile(biomassForecastTwoyear(), .975)
      )

  })

  # tells download handler which plots to download ---- 
  observeEvent(input$doCalculations, {

    # plot NAA agecomp? ----
    naaAgecompSum <- sum(na.omit(naaAgecompObserved()$Proportion))
    if((abs(naaAgecompSum-1) < .01) | all(is.na(naaAgecompInput()))) {
      rv$naaAgecompPlotShown <- TRUE
    } else {
      rv$naaAgecompPlotShown <- FALSE
    }

    # plot BAA agecomp? ----
    naaAgeCompZeroToNA <- ifelse(naaAgecompInput() == 0, NA, naaAgecompInput())
    currYearWAAZeroToNA <- ifelse(currYearWAA() == 0, NA, currYearWAA())
 
    if (!all(is.na(naaAgeCompZeroToNA)) & !all(is.na(currYearWAAZeroToNA))) {
   
      naaAgeCompHaveObservations <- which(is.na(naaAgeCompZeroToNA))
      currYearWAAHaveObservations <- which(is.na(currYearWAAZeroToNA))

      if(all(naaAgeCompHaveObservations %in% currYearWAAHaveObservations)) {
        rv$baaAgecompPlotShown <- TRUE
      } else {
        rv$baaAgecompPlotShown <- FALSE
      }
      
      if(all(currYearWAAHaveObservations %in% naaAgeCompHaveObservations)) {
        rv$baaAgecompPlotShown <- TRUE
      } else {
        rv$baaAgecompPlotShown <- FALSE
      }

    } 

    # plot MDM or 2-year biomass forcast histograms? ----
    if (input$doCalculations > 0) {
      rv$mdmPlotShown <- TRUE
      rv$biomassPlotTwoyearShown <- TRUE
    }

  })

  # download table ----
  output$downloadTbl <- downloadHandler(

    filename = function() {
      paste0("PWS-herring-mid-year-management-estimates-", Sys.Date(), ".csv")
    },

    content = function(file) {
      write.csv(rbind(rv$tbl, agecompTbl), file, row.names = FALSE)
    }
  
  )

  # download table ----
  output$downloadFigs <- downloadHandler(

    filename = function() {
      paste0("PWS-herring-mid-year-management-figures-", Sys.Date(), ".zip")
    },

    content = function(file) {
      
      temp_dir <- paste0(tempdir(), "/plots")
      unlink(temp_dir, recursive = TRUE)
      dir.create(temp_dir)

      ggsave(paste0(curr_year, "-mature-biomass-forecast.png"), plot = biomass_forecast_plot, path = temp_dir, device = "png")
      if(rv$naaAgecompPlotShown) {
        ggsave("naa-agecomp.png", plot = naaAgecompPlot(), path = temp_dir, device = "png")
      }
      if(rv$baaAgecompPlotShown) {
        ggsave("baa-agecomp.png", plot = baaAgecompPlot(), path = temp_dir, device = "png")
      }
      if(rv$mdmPlotShown) {
        ggsave("mile-days-milt-forecast.png", plot = mdmPlot(), path = temp_dir, device = "png")
      }
      if(rv$biomassPlotTwoyearShown) {
        ggsave(paste0(curr_year+1, "-mature-biomass-forecast.png"), plot = biomassPlotTwoyear(), path = temp_dir, device = "png")
      }

      zipr(file, files = paste(temp_dir, list.files(temp_dir), sep = "/"))
    
      unlink(temp_dir, recursive = TRUE)

    }
  )

  # prints output to UI for debugging purposes ----
  # output$debug <- renderText({
  #   rv$naaAgecompPlotShown
  # })

}

# open shiny app ----
shinyApp(ui = ui, server = server)
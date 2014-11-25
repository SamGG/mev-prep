# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
# http://rstudio.github.io/shiny/tutorial/#uploads

library(RColorBrewer)

source("tdms.R")

library(shiny)

shinyServer(function(input, output, session) {

  # Read the input file
  # TODO: check that file is TDMS; otherwise output an error
  rec_data <- reactive({
    tdms <- NULL
    if (is.null(input$file$datapath)) return(NULL)
    tdms <- read.tdms(input$file$datapath)
    tdms
  })
  
  # Table of sample's parameters
  # TODO: populate a menu of paramters within certain constraints
  pData <- reactive({
    tdms <- rec_data()
    if (is.null(tdms)) return(NULL)
    tdms$pData
  })

  # Process, ie transform, center, scale
  proc_data <- reactive({
    tdms <- rec_data()
    if (is.null(tdms)) return(NULL)
    # Transform
    transf_param = as.numeric(input$transf_param)
    if (transf_param < 0 | is.na(transf_param)) transf_param = 0
    tdms = switch(input$transf,
                  "none" = tdms,
                  "log2" = tdms.transform(tdms, "log2", transf_param),
                  "asinh" = tdms.transform(tdms, "asinh", transf_param),
                  ... = NULL
    )
    # Center
    tdms = switch(input$center,
                  "mean" = tdms.center(tdms, method = "mean"),
                  "median" = tdms.center(tdms, method = "median"),
                  "minmax" = tdms.center(tdms, method = "minmax"))
    # Scale
    sc.perc = as.numeric(input$scale_perc)
    if (sc.perc < 0 | is.na(sc.perc) | !is.numeric(sc.perc)) sc.perc = 5
    sc.perc = sc.perc / 100
    sc.group = NA
    sc.filter = "none"
    sc.min.size = 3
    # Return
    tdms = tdms.scale(tdms, method = input$scale, percent = sc.perc,
                      group = sc.group, filter = sc.filter, min.size = sc.minsize)
    tdms
  })
  
  # TODO: add logging
  # TODO: manage suffix for saving
  # TODO: re-add download
  # TODO: use tabset
  
  # Post-processing
  post_data <- reactive({
    tdms = proc_data()
    if (is.null(tdms)) return(NULL)
    if (input$summarize != "0") {
      probs = switch(input$summarize,
                     "med" = c(0.50),
                     "5"   = c(0.05, 0.50, 0.95),
                     "20"  = c(0.20, 0.50, 0.80))
      tdms = tdms.summarize(tdms, group = tdms$pData$Strain, probs = probs)
      # suffix = "-summ"
    }
    if (input$sort_within == "1") {
      # tdms = tdms.sort(tdms, group = tdms$pData$Strain, probs = probs)
      # suffix = "-sort"
    }
    
    if (input$transpose != "0") {
      tdms = tdms.t(tdms)
      # suffix = "-trps"
    }
    tdms
  })

  # Check consistency
  observe({
    # If minimax, no grouping
    if (input$scale == "minmax" & input$scale_group != "none")
      updateSelectInput(session, "scale_group", selected = "none")
    # File input
  })

  # Various outputs
  output$textOutp <- renderText({ input$scale })

  output$summary <- renderText({ input$action; input$file$datapath; str(rec_data()) })
  
  output$summ <- renderTable({ pData() })

  # ???
  output$contents <- renderText({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    inFile$datapath
  })

  # Download
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$dataset, '.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput(), file)
    }
  )

  # Heatmap view
  # TODO: add a toggle between raw/proc/post?
  output$heatmap <- renderPlot({
    #tdms <- rec_data()
    #tdms <- proc_data()
    tdms <- post_data()
    if (is.null(tdms)) return(NULL)
    col = rev(brewer.pal(11, "RdYlBu"))
    #heatmap(tdms$exprs, Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)
    # TODO: limit rows/columns
    .exprs = tdms$exprs
    maxr = seq(min(nrow(.exprs), 20))
    maxc = seq(min(ncol(.exprs), 20))
    heatmap(.exprs[maxr,maxc], Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)
  })

})


# TODO: order all columns by factor
# TODO: final rescaling, add an histogram and a colorscale below
# TODO: ensure that zero is at the middle of the scale

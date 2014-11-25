# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("CytoPreparePercent"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      fileInput("file", label = strong("Input file"), multiple = FALSE, accept = NULL),

      # Shortcuts
      selectInput("shortcut", label = strong("Shortcuts"),
                  choices = list("Manual settings" = 1),
                  selected = 1),

      # Transform
      selectInput("transf", label = strong("Transform"),
                  choices = list("none",
                                 "log2",
                                 "asinh"),
                  selected = "none"),
      textInput("transf_param", label = "soft floor",
                value = "0.1"),

      # Center
      selectInput("center", label = strong("Center"),
                  choices = list("mean",
                                 "median",
                                 "minimax (middle)" = "minmax"),
                  selected = "minmax"),

      # Scale
      selectInput("scale", label = strong("Scale"),
                  choices = list("sd (std dev)" = "sd",
                                 "mad (med abs diff)" = "mad",
                                 "minimax (half scale)" = "minmax"),
                  selected = "minmax"),
      selectInput("scale_group", label = "Use groups:",
                  choices = list("None (discovery)" = "none",
                                 "All" = "all",
                                 "Exclude Ref" = "excl.ref",
                                 "Only Ref" = "only.ref"),
                  selected = "none"),
      sliderInput("scale_perc", label = "Percentile (minimax):",
                  min = 0, max = 25, value = 5),

      hr(),
      #downloadButton('downloadData', 'Download'),
      
      # Post-processing
      radioButtons("summarize", label = strong("Summarize"),
                   choices = list("No" = 0, "Median" = "med",
                                  "5-50-95%" = 5, "20-50-80" = 20), selected = 0),

      radioButtons("sort_within", label = strong("Sort within groups"),
                   choices = list("No" = 0, "Yes" = 1), selected = 0, inline = TRUE),      

      radioButtons("transpose", label = strong("Transpose finally"),
                   choices = list("No" = 0, "Yes" = 1), selected = 0, inline = TRUE)
    ),

    # Show a plot of the generated distribution
    mainPanel(
      htmlOutput("summary"),
      tableOutput("summ"),
      verbatimTextOutput("textOutp"),
      htmlOutput('contents'),
      plotOutput('heatmap')
    )
  )
))

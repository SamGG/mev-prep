# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

# Conditional content in SidebarPanel
# http://stackoverflow.com/questions/16225260/are-tabs-with-their-own-sidebars-and-mainpanels-possible-in-shiny

library(shiny)

# small textInput boxes side-by-side
# http://stackoverflow.com/questions/20637248/shiny-4-small-textinput-boxes-side-by-side
textInputRow <- function (inputId, label, value = "") {
  div(style="display:inline-block",
      tags$label(label, `for` = inputId),
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}

shinyUI(fluidPage(
  titlePanel("mev-prep - Prepare multivariate data for MeV"),
  sidebarLayout(
  sidebarPanel( witdh=6,

conditionalPanel(condition="input.conditionedPanels==1",
  selectInput("demo_file", label = "Upload a file or select a demo",
             choices = c("Upload", names(tdms.examples()))),
  conditionalPanel(
    condition = "input.demo_file == 'Upload'",
    fileInput("file", label = NA, multiple = FALSE, accept = NULL)),
  br(),
  #tags$label(strong("Grouping factor")),
  uiOutput("factorSelectInput")
),

conditionalPanel(condition="input.conditionedPanels==3",
  # Shortcuts
  #                      selectInput("shortcut", label = strong("Shortcuts"),
  #                                  choices = list("Manual settings" = 1),
  #                                  selected = 1),
  # Transform
  selectInput("transf", label = strong("Transform"),
             choices = list("none",
                            "log2",
                            "asinh"),
             selected = "none"),
  textInput("transf_param", label = "soft floor",
           value = "0.1"),
  # Center
  selectInput("center", label = tagList(strong("Center"), em("line"), "using global"),
             choices = list("none",
                            "mean",
                            "median",
                            "minimax (middle)" = "minmax"),
             selected = "median"),
  # Scale
  selectInput("scale", label = tagList(strong("Scale"), em("line"), "using "),
             choices = list("none"
                            , "sd (std dev)" = "sd"
                            , "mad (med abs diff)" = "mad"
                            #, "minimax (half scale)" = "minmax"
                            ),
             selected = "sd"),
  conditionalPanel(
   condition = "input.scale != 'minmax' && input.scale != 'none'",
   selectInput("scale_group", label = "onto groups",
               choices = list("None (discovery)" = "none"
                              , "Every" = "all"
#                             , "Pool" = "pool"
#                             , "Exclude Ref" = "excl.ref",
#                             , "Ref. only" = "only.ref"
                              ),
               selected = "none")
  ),
  conditionalPanel(
   condition = "input.scale == 'minmax'",
   sliderInput("scale_perc", label = "Percentile (minimax):",
               min = 0, max = 25, value = 5)
  ),
# Set Zero
uiOutput("color.zero.dialog")
#,textInput("color.zero.value", label = "Set zero at",          value = "0")
),

conditionalPanel(condition="input.conditionedPanels==4",
#  checkboxInput("adjust.scale.auto", label = "Automatically adjust scale", value = TRUE),
  selectInput("scale_adjust", label = "Color scale limits",
              choices = list("Out-of-scale data percentage" = "percent"
                             , "Manuel setting" = "manual"
              ),
              selected = "percent"),
  conditionalPanel(condition="input.scale_adjust=='percent'",
  sliderInput("adjust.scale.perc", label = NA,
              min = 0, max = 15, value = 3)),
  textInputRow("optZmin", label = "Data Min.", value = "-3"),
  textInputRow("optZmax", label = "Data Max.", value = "3"),
#  tags$label(strong("Color scale zero")),
  checkboxInput("zero_center_palette", "Bind zero to center of color scale", value = TRUE),
  # Scale Color
  selectInput("color.transf", label = strong("Color function"),
             choices = list("square","linear", "sqrt"),
             selected = "linear")
),

conditionalPanel(condition="input.conditionedPanels==5",
#                      radioButtons("sort_within", label = strong("Sort within groups"),
#                                   choices = list("No" = 0, "Yes" = 1), selected = 0, inline = TRUE),
  selectInput("summarize", label = strong("Summarize"),
             choices = list("No" = 0, "Median" = "med",
                            "5-50-95%" = 5, "20-50-80%" = 20, "5-20-50-80-95%" = 205), selected = 0),
  checkboxInput("summarize.interleave", label = "Interleave",
               value = FALSE),
  radioButtons("transpose", label = strong("Transpose"),
              choices = list("No" = 0, "Yes" = 1), selected = 0, inline = TRUE),
  br(),tags$label(strong("Download")),
  downloadButton("downloadData", label = "Download")
),

# conditionalPanel(condition="input.conditionedPanels==6",
#                  selectInput("cl.distance", label = strong("Distance"),
#                              choices = list("Eucl" = "eucl",
#                                             "Eucl(sqrt)" = "eucl_sqrt",
#                                             "Eucl(sqrt(sqrt))" = "eucl_sqrt_sqrt"
#                                             ), selected = "eucl_sqrt"),
#                  br(),tags$label(strong("Download")),
#                  downloadButton("cl.downloadData", label = "Download")
# ),

conditionalPanel(condition="input.conditionedPanels==8",
  selectInput("palette", label = strong("Color scheme"), selected = "gne BuWhRd",
             choices = list("gne BuWhRd", "BuRd", "BuYlRd", "GnYlRd", "BrOrYl", "Oranges", "BGBr", "mev BuBkYl", "mev GrBkRd")),
  tags$label(strong("Data columns odering")),
  checkboxInput("sort_byfactor", "Auto-sort columns according grouping factor", value = TRUE),
  tags$label(strong("Data selection")),br(),
  textInputRow("data.sel.rows", label = "Rows",
           value = "40"),
  textInputRow("data.sel.cols", label = "Cols",
           value = "20")
),

conditionalPanel(condition="input.conditionedPanels==9",
  helpText(strong("Help"))
)

  ),

  mainPanel(
    #titlePanel("cytoPrepare - Normalize percentages"),
    tabsetPanel(id = "conditionedPanels", type = "pills",

      tabPanel("Input", value=1,
               htmlOutput("input_info"), br(),
               strong("Table of factors"),
               tableOutput("input_pData"),
               strong("Table of frequencies"),
               tableOutput("input_factorFreq")
      ),

      tabPanel("Adjust", value=3,
               textOutput("normalize_info"),
               strong("Heatmap of adjusted data"),
               plotOutput("heatmapNorm"),
               strong("Notice"),
               htmlOutput("noticeNorm")
      ),

      tabPanel("Colorize", value=4,
               strong("Heatmap of colorized data"),
               plotOutput("heatmapCalib"),
               plotOutput("colorbarCalib", height = "50px"),
               strong("Notice"),
               htmlOutput("noticeCalib")
      ),

      tabPanel("Finalize", value=5,
               strong("Heatmap of finalized data"),
               plotOutput("heatmapPostProc"),
               strong("Notice"),
               htmlOutput("noticePostProc")
      ),

#       tabPanel("Clustering", value=5,
#                strong("Heatmap of clusterized data"),
#                plotOutput("heatmapCluster")),

      tabPanel("Options", value=8,
               strong("Notice"),
               htmlOutput("noticeOptions")
      ),

      tabPanel("Help", value=9,
               strong("Notice"),
               htmlOutput("noticeHelp")
      )

  )

))))

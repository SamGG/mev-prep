# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
# http://rstudio.github.io/shiny/tutorial/#uploads
# http://shiny.rstudio.com/reference/shiny/latest/fileInput.html

library(RColorBrewer)

source("tdms.R")

# TODO: check 'use group' with 'grouping factor'



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
# TODO: populate a menu of paramaters within certain constraints
# ie at least 5 points to compute sd or merge all sd
# TODO: define a reference group for zero and sd
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
  if (transf_param < 0 || is.na(transf_param)) transf_param = 0
  tdms = switch(input$transf,
                "none" = tdms,
                "log2" = tdms.transform(tdms, "log2", transf_param),
                "asinh" = tdms.transform(tdms, "asinh", transf_param),
                ... = NULL
  )
  # Center
  tdms = switch(input$center,
                "none"   = tdms.center(tdms, method = "none"),
                "mean"   = tdms.center(tdms, method = "mean"),
                "median" = tdms.center(tdms, method = "median"),
                "minmax" = tdms.center(tdms, method = "minmax"))
  # Scale
  sc.perc = as.numeric(input$scale_perc)
  if (sc.perc < 0 | is.na(sc.perc) | !is.numeric(sc.perc)) sc.perc = 5
  sc.perc = sc.perc / 100
  sc.group = input$factorSelected
  if (input$scale_group == "none" || input$scale == "none" || input$scale == "minmax") {
    sc.group = NA
  }
  sc.filter = input$scale_group # "none" # TODO!
  print(sc.filter)
  sc.min.size = 3 # TODO!
  # Return
  tdms = tdms.scale(tdms, method = input$scale, percent = sc.perc,
                    group = sc.group, filter = sc.filter, min.size = sc.minsize)
  if (input$sort_byfactor) 
    tdms = tdms.sort.byfactor(tdms, input$factorSelected)
  tdms
})

calib_data <- reactive({
  tdms = proc_data()
  if (is.null(tdms)) return(NULL)
  # TODO: get the mean of the group
  print(input$color.zero)
  if (input$color.zero != "none")
    tdms = tdms.center.ref(tdms, method = "median", grouping = input$factorSelected, group = input$color.zero)
  print(round(tdms$exprs[14:17,], 2))
  
  # Apply shift
  color.zero.value = as.numeric(input$color.zero.value)
  .exprs = tdms$exprs - color.zero.value
  # TODO: heatmap with a fixed zero
  # Color function
  .sign  = sign(.exprs)
  .exprs = .sign * .exprs
  .exprs = switch(input$color.transf,
         linear = .exprs,
         square = .exprs^2,
         sqrt   = sqrt( .exprs ),
         sqrt_sqrt = sqrt( sqrt( .exprs ) ),
         ...    = .exprs )
  .exprs = .sign * .exprs
  tdms$exprs = .exprs
  # Done
  tdms
})


# Post-processing
post_data <- reactive({
  tdms = proc_data()
  if (is.null(tdms)) return(NULL)
  if (input$summarize != "0") {
    probs = switch(input$summarize,
                   "med" = c(0.50),
                   "5"   = c(0.05, 0.50, 0.95),
                   "20"  = c(0.20, 0.50, 0.80),
                   "205" = c(0.05, 0.20, 0.50, 0.80, 0.95) )
    tdms = tdms.summarize(tdms, group = input$factorSelected, probs = probs, inter.leave = input$summarize.interleave)
    # suffix = "-summ"
    # if (input$summarize.interleave) suffix = "-summ-inter"
  }
#     if (input$sort_within == "1") {
#       # tdms = tdms.sort(tdms, group = tdms$pData$Strain, probs = probs)
#       #suffix = "-sort"
#     }
  
  if (input$transpose != "0") {
    tdms = tdms.t(tdms)
    # suffix = "-trps"
  }
  tdms
})

# Generic heatmap
.heatmap <- function(tdms, Rowv = NA) {
  if (is.null(tdms)) return(NULL)
  .exprs = tdms$exprs
  .pct = as.numeric(input$adjust.scale.perc) / 100
  print(quantile(.exprs, probs = c(.pct, 1-.pct), na.rm = T))
  # Adjust color scale
  if (input$adjust.scale.auto && !input$zero_center_palette) {
    .pct = as.numeric(input$adjust.scale.perc) / 100
    .rng = quantile(.exprs, probs = c(.pct, 1-.pct), na.rm = T)
    .exprs = ((.exprs - .rng[1]) / (.rng[2] - .rng[1]) - 0.5) * 6  # scale to -3, +3
  }
  print(quantile(.exprs, probs = c(.pct, 1-.pct), na.rm = T))
  #
  if (input$zero_center_palette) {
    .pct = as.numeric(input$adjust.scale.perc) / 100
    .rng = quantile(.exprs, probs = c(.pct, 1-.pct), na.rm = T)
    .abs = max(abs(.rng))
    .rng = c( -.abs, .abs )
    .exprs = ((.exprs - .rng[1]) / (.rng[2] - .rng[1]) - 0.5) * 6  # scale to -3, +3
  }
  print(quantile(.exprs, probs = c(.pct, 1-.pct), na.rm = T))
  # Reduce the size of frame
  maxr = seq(min(nrow(.exprs), optMaxRow()))
  maxc = seq(min(ncol(.exprs), optMaxCol()))
  .exprs = .exprs[maxr,maxc]
  # Remove out of scale values as they appear as missing
  .exprs[.exprs > optZmax()] = optZmax()
  .exprs[.exprs < optZmin()] = optZmin()
  # Update the color scale
  zlim = c(optZmin(), optZmax())
  # Reset the character size
  cex = 0.2 + 1/log10(max(maxr, maxc))
  # Plot
  heatmap(.exprs, Rowv = Rowv, Colv = NA, col = colPal(), revC = TRUE, 
          cexCol = cex, cexRow = cex, scale = "none", 
          zlim = c(optZmin(), optZmax()))
}

# Check consistency
observe({
  # If minimax, no grouping
  # No factor, no summarization possible
  if ((is.null(input$factorSelected) || input$factorSelected == "none") && 
      (is.null(input$summarize) || input$summarize != "0")) {
    updateSelectInput(session, "summarize", selected = "0")
    updateCheckboxInput(session, "summarize.interleave", value = FALSE) # TODO: not acting
  }
  # File input
  # Transform parameter can't be empty
  #if (input$transf_param == "") input$transf_param = 0.1
})

# Various outputs
output$textOutp <- renderText({ input$scale })

output$summary <- renderText({ input$action; input$file$datapath; str(rec_data()) })

##### INPUT TAB ######

# Parameter table
output$input_pData <- renderTable({ pData() })

# Frequency table
output$input_factorFreq <- renderTable({
  .pData = pData()
  if (is.null(.pData)) NULL
  z = lapply(.pData[-1], summary, maxsum=7L)
  # Build a data frame of fraquency; ugly but avoid dependency to plyr
  do.call("rbind",
          lapply(names(z), 
                 function(i) { 
                   do.call("rbind", 
                           lapply(names(z[[i]]), 
                                  function(j) 
                                    data.frame(factor=i, value=j,
                                               freq=z[[i]][j]) )) } ))
  }) #, caption = "Table of frequeny of parameters")

# Information
output$input_info <- renderText({
  if (is.null(rec_data()))
    return("Please select a TDMS data file<br>")
  c(sprintf("File is \"%s\"<br>", input$file$name))
})

##### NORMALIZING TAB ######

# Heatmap view
output$heatmapNorm <- renderPlot({
  # Get the data
  tdms <- proc_data()
  .heatmap(tdms)
})

# Updating the factor selection popup
output$factorSelectInput <- renderUI({
  .pData <- pData()
  factors <- list("None available" = "none")
  if (!is.null(.pData) && ncol(tdms$pData > 1))
    factors = c("None" = "none", colnames(.pData)[-1])
  selectInput("factorSelected", label = strong("Grouping factor"), choices = factors)
})

# Notice
output$noticeNorm <- reactive({
  .text = "This panel allows transforming, then centering, then scaling."
  .text = c(.text, switch(input$transf,
      none = "Transform: no transformation is applied.",
      log2 = "Transform: log2( raw + soft_floor ). The soft_floor parameter avoid the over-spreading of low levels.",
      asinh = "Transform: asinh( raw / soft_floor ). The soft_floor parameter avoid the over-spreading of low levels.") )
  .text = c(.text, switch(input$center,
      none = "Center: no transformation is applied.",
      ...  = "Centering") )
  if (input$scale!="none") {
    .text = c(.text, "Scaling")
    if (input$scale == 'minmax') {
      .text = c(.text, "Minmax")
    } else {
      .text = c(.text, "Not Minmax")
    }
  }
  .text = c(.text, "Grouping factor: if a factor is selected, groups are deduced and the data ")
  if (input$sort_byfactor) {
    .text = c(.text, "The columns are ordered by the selected factor.")
  }
  .text = c(.text, "Heatmap colors are adjusted in the Options tab.")
  if (is.null(.text)) return("")
  .text = paste(.text, collapse = "</li><li>")
  sprintf("<ul><li>%s</li></ul>", .text)
})

# Information
output$normalize_info <- renderText({
  ""
})

# TODO: cluster rows
# TODO: color scale coefficient, slider?
# TODO: check color scale is finite

# TODO: Adjust color: constrast slider or auto-adjust -3%, +3%
# TODO: tab color adjust?


##### CALIBRATION TAB #####

# Heatmap view
output$heatmapCalib <- renderPlot({
  # Get the data
  print(round(tdms$exprs[14:16,], 2))
  tdms <- calib_data()
  print(round(tdms$exprs[14:16,], 2))
  .heatmap(tdms)
})

# Notice
output$noticeCalib <- reactive({
  .text = "This panel allows transforming, then centering, then scaling."
  .text = c(.text, switch(input$transf,
                          none = "Transform: no transformation is applied.",
                          log2 = "Transform: log2( raw + parameter ). The parameter avoid the over-spreading of low levels.",
                          asinh = "Transform: asinh( raw / parameter ). The parameter avoid the over-spreading of low levels.") )
  .text = c(.text, switch(input$center,
                          none = "Center: no transformation is applied.",
                          ...  = "Centering") )
  if (input$scale!="none") {
    .text = c(.text, "Scaling")
    if (input$scale == 'minmax') {
      .text = c(.text, "Minmax")
    } else {
      .text = c(.text, "Not Minmax")
    }
  }
  .text = c(.text, "Grouping factor: if a factor is selected, groups are deduced and the data ")
  if (input$sort_byfactor) {
    .text = c(.text, "The columns are ordered by the selected factor.")
  }
  .text = c(.text, "Heatmap colors are adjusted in the Options tab.")
  if (is.null(.text)) return("")
  .text = paste(.text, collapse = "</li><li>")
  sprintf("<ul><li>%s</li></ul>", .text)
})

# Updating the factor selection popup
output$color.zero.dialog <- renderUI({
  .pData = pData()
  # if (is.null(.pData)) NULL
  z = NULL
  if (!is.null(input$factorSelected) && input$factorSelected != "none")
    z = unique(.pData[ ,input$factorSelected])
  print(z)
  z = as.character(z)
  # TODO: filter to get enough sample per group
  groups <- list("None available" = "none")
  if (!is.null(.pData) && length(z))
    groups = c("None" = "none", z)
  print(groups)
  selectInput("color.zero", label = tagList(strong("Zero"), em("line for")),
              choices = groups)
})


output$colorbarCalib <- renderPlot({
  old.mar = par("mar")
  par(mar = c(2, 20, 0.1, 20))
  image(x=seq(optZmin(), optZmax(), length.out = 254), z = matrix(1:254, nc = 1),
        col = colPal(), xaxp=c(-3, 3, 6), axes = FALSE, xlab = "")
  axis(1, at = seq(optZmin(), optZmax(), 1), cex.axis = 0.8)
  par(mar = old.mar)
})

##### CLUSTERING TAB ######

# TODO: finalize
# Heatmap view of clustering
output$heatmapCluster <- renderPlot({
  # Get the data
  tdms <- proc_data()
  .exprs = tdms$exprs
  .exprs = sign(.exprs) * sqrt(abs(.exprs))
  .dist = dist(.exprs, method = "euclidean")
  .rowv = as.dendrogram( hclust(.dist, method = "average") )
  .heatmap(tdms, .rowv)
})


##### POST-PROCESSING TAB ######

# Heatmap view of post-processing
output$heatmapPostProc <- renderPlot({
  # Get the data
  tdms <- post_data()
  .heatmap(tdms)
})

# Notice
output$noticePostProc <- reactive({
  .text = NULL
  if (input$summarize!="0") {
    .text = c(.text, "Summarization allows comparing the distribution of values between groups: quantiles are computed and displayed with color. If a shift is observed, it might be meaningful. Summarization is a color representation of boxplots.")
  } else {
    .text = c(.text, "To apply Summarization, a grouping factor must be selected (Normalize panel).")
  }
  if (input$summarize.interleave) {
    .text = c(.text, "Interleave option organizes summarization by quantiles instead of groups. A pattern alternating high and low values means that quantiles have very different levels in each group. It might be meaningful.")
  }
  if (is.null(.text)) return("")
  .text = paste(.text, collapse = "</li><li>")
  sprintf("<ul><li>%s</li></ul>", .text)
})

# Download
output$downloadData <- downloadHandler(
  filename = function() {
    if (is.null(post_data())) return(NULL)  # TODO: does not return if no file loaded
    # TODO: manage suffix for saving
    .filename = input$file$name
    .filename = gsub(.filename, '/\\..{3,4}+$/', '', perl = TRUE)
    paste(.filename, '-', Sys.Date(), '.txt', sep='')
  },
  content = function(file) {
    tdms <- post_data()
    if (is.null(tdms)) return(NULL)
    write.tdms(tdms, file)
  }
)

##### CONTROL-OPTIONS TAB ######

# Notice
output$noticeOptions <- reactive({
  .text = "This panel allows adjusting various options."
  .text = c(.text, "Various Brewer diverging palettes are available.")
  if (input$adjust.scale.auto) {
    .text = c(.text, "The color range is adjusted to the data automatically. A user defined quantile is allowed at both end, leading to color saturation.")
  }
  .text = c(.text, "To reach a better interactivty, only a small part of the data is displayed. This data selection is defined as rows and columns from the top left sheet.")
  if (is.null(.text)) return("")
  .text = paste(.text, collapse = "</li><li>")
  sprintf("<ul><li>%s</li></ul>", .text)
})

# Color palette
colPal <- reactive({
  # Palettes are obtained from the library(RColorBrewer) using the following code
  # cat("c(\"", paste(brewer.pal(9, "Oranges"), collapse = '\",\"'), "\")", sep="")
  .cp = switch(input$palette,
               "BuRd" = rev(c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")),
               "BuYlRd" = rev(c("#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4")),
               "GnYlRd" = rev(c("#D73027", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#D9EF8B", "#A6D96A", "#66BD63", "#1A9850")),
               "BrOrYl" = rev(c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")),
               "Oranges" = rev(c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704")),
               "BGBr" = c("#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E")
               )
  colorRampPalette(.cp, space = "Lab")(64)
})

optZmin <- function() -3
optZmax <- function() +3

# TODO: add color scale to heatmaps

# TODO: add logging
# TODO: limit rows/columns
optMaxRow <- function() 40
optMaxCol <- function() 20

# TODO: windsoring

# TODO: final rescaling, add a tab?
#       add an histogram and a colorscale below
#       ensure that zero is at the middle of the scale

})

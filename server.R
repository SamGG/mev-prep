# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
# http://rstudio.github.io/shiny/tutorial/#uploads
# http://shiny.rstudio.com/reference/shiny/latest/fileInput.html

library(RColorBrewer)

library(shiny)

shinyServer(function(input, output, session) {

# Read the input file
# TODO: check that file is TDMS; otherwise output an error
rec_data <- reactive({
  tdms <- NULL
  if (input$demo_file == "Upload" && !is.null(input$file$datapath))
    tdms <- read.tdms(input$file$datapath)
  if (input$demo_file %in% names(tdms.examples())) # !is.null(input$demo_file) &&
    tdms <- read.tdms(tdms.examples()[[input$demo_file]]$file)
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
  sc.min.size = 3 # TODO: put in options?
  tdms = tdms.scale(tdms, method = input$scale, percent = sc.perc,
                    group = sc.group, filter = sc.filter, min.size = sc.minsize)
  if (input$sort_byfactor)
    tdms = tdms.sort.byfactor(tdms, input$factorSelected)

  # TODO: get the mean of the group
  print(input$color.zero)
  if (input$color.zero != "none")
    tdms = tdms.center.ref(tdms, method = "median", grouping = input$factorSelected, group = input$color.zero)

  tdms
})

calib_data <- reactive({
  tdms = proc_data()
  if (is.null(tdms)) return(NULL)

  # Apply shift
#   color.zero.value = as.numeric(input$color.zero.value)
#   .exprs = tdms$exprs - color.zero.value
  # Color function
  .exprs = tdms$exprs
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
  tdms = calib_data()
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

#optZmin <- -3
#optZmax <- +3

tdms = reactive({
  print(c("Panel = ", input$conditionedPanels))
  if (input$conditionedPanels == 3)
    proc_data()
  else if (input$conditionedPanels == 4)
    calib_data()
  else if (input$conditionedPanels == 5)
    post_data()
  else
    NULL
})

optRange = reactive({
  if (input$scale_adjust == 'percent') {
    .pct = as.numeric(input$adjust.scale.perc) / 100
    print(.pct)
    .exprs = tdms()$exprs
    .rng = quantile(.exprs, probs = c(.pct, 1-.pct), na.rm = T)
    if (input$zero_center_palette) {
      .abs = max(abs(.rng))
      .rng = c( -.abs, .abs )
    }
    .rng = round(as.numeric(.rng), 2)
    updateTextInput(session, "optZmin", value = .rng[1])
    updateTextInput(session, "optZmax", value = .rng[2])
  }
  .rng
})

optZmin <- reactive({
   if (input$scale_adjust == 'percent')
    .zmin = optRange()[1]
   else {
     .zmin = as.numeric(input$optZmin)
     .zmax = as.numeric(input$optZmax)
     if (input$zero_center_palette) {
       .zmin = -abs(as.numeric(input$optZmin))
       updateTextInput(session, "optZmax", value = abs(.zmin))
     }
     else if (.zmin >= .zmax) {
       updateTextInput(session, "optZmin", value = .zmax * 0.99)
     }
   }
   .zmin
})

optZmax <- reactive({
   if (input$scale_adjust == 'percent')
    .zmax = optRange()[2]
   else {
     .zmin = as.numeric(input$optZmin)
     .zmax = as.numeric(input$optZmax)
     if (input$zero_center_palette) {
       .zmax = abs(as.numeric(input$optZmin))
       updateTextInput(session, "optZmax", value = .zmax)
     }
     else if (.zmax <= .zmin) {
       updateTextInput(session, "optZmax", value = .zmin * 1.01)
     }
   }
   .zmax
})

# TODO: add color scale to heatmaps


# Generic heatmap
.heatmap <- function(Rowv = NA) {
  if (is.null(tdms())) return(NULL)
  .exprs = tdms()$exprs
  print(is.numeric(.exprs))
  optZmin = optZmin()
  optZmax = optZmax()
    print(c(optZmin, optZmax))
  # Reduce the size of frame
  if (optMaxRow() == 0)
    maxr = seq(nrow(.exprs))
  else if (optMaxRow()>0)
    maxr = seq(1, min(nrow(.exprs), optMaxRow()))
  else
    maxr = nrow(.exprs) + 1 - seq(min(nrow(.exprs), -optMaxRow()), 1)
  if (optMaxCol() == 0)
    maxc = seq(ncol(.exprs))
  else if (optMaxCol()>0)
    maxc = seq(1, min(ncol(.exprs), optMaxCol()))
  else
    maxc = ncol(.exprs) + 1 - seq(min(ncol(.exprs), -optMaxCol()), 1)
  .exprs = .exprs[maxr, maxc]
  # Remove out of scale values as they appear as missing
  .exprs[.exprs > optZmax] = optZmax
  .exprs[.exprs < optZmin] = optZmin
  # Update the color scale
  zlim = c(optZmin, optZmax)
  # Reset the character size
  cex = 0.2 + 1/log10(max(maxr, maxc))
  # Colorize columns
  # TODO: grouping = input$factorSelected, group = input$color.zero
  # Plot
  heatmap(.exprs, Rowv = Rowv, Colv = NA, col = colPal(), revC = TRUE,
          cexCol = cex, cexRow = cex, scale = "none",
          #          ColSideColors = NULL,
          zlim = zlim)
}


# Generic heatmap
.heatmap2 <- function(tdms, Rowv = NA) {
  if (is.null(tdms)) return(NULL)
  .exprs = tdms$exprs
  .pct = as.numeric(input$adjust.scale.perc) / 100
  print(quantile(.exprs, probs = c(.pct, 1-.pct), na.rm = T))
  # Adjust color scale
  if (input$scale_adjust == "percent") {
    if (input$zero_center_palette) {
      .pct = as.numeric(input$adjust.scale.perc) / 100
      .rng = quantile(.exprs, probs = c(.pct, 1-.pct), na.rm = T)
      .abs = max(abs(.rng))
      .rng = c( -.abs, .abs )
      #.exprs = ((.exprs - .rng[1]) / (.rng[2] - .rng[1]) - 0.5) * 6  # scale to -3, +3
    } else {
      .pct = as.numeric(input$adjust.scale.perc) / 100
      .rng = quantile(.exprs, probs = c(.pct, 1-.pct), na.rm = T)
      #.exprs = ((.exprs - .rng[1]) / (.rng[2] - .rng[1]) - 0.5) * 6  # scale to -3, +3
    }
    optZmin(.rng[1])
    optZmax(.rng[2])
    print(quantile(.exprs, probs = c(.pct, 1-.pct), na.rm = T))
  }
  # Reduce the size of frame
  if (optMaxRow() == 0)
    maxr = seq(nrow(.exprs))
  else if (optMaxRow()>0)
    maxr = seq(1, min(nrow(.exprs), optMaxRow()))
  else
    maxr = nrow(.exprs) + 1 - seq(min(nrow(.exprs), -optMaxRow()), 1)
  if (optMaxCol() == 0)
    maxc = seq(ncol(.exprs))
  else if (optMaxCol()>0)
    maxc = seq(1, min(ncol(.exprs), optMaxCol()))
  else
    maxc = ncol(.exprs) + 1 - seq(min(ncol(.exprs), -optMaxCol()), 1)
  .exprs = .exprs[maxr, maxc]
  # Remove out of scale values as they appear as missing
  .exprs[.exprs > optZmax()] = optZmax()
  .exprs[.exprs < optZmin()] = optZmin()
  # Update the color scale
  zlim = c(optZmin(), optZmax())
  # Reset the character size
  cex = 0.2 + 1/log10(max(maxr, maxc))
  # Colorize columns
  # TODO: grouping = input$factorSelected, group = input$color.zero
  # Plot
  heatmap(.exprs, Rowv = Rowv, Colv = NA, col = colPal(), revC = TRUE,
          cexCol = cex, cexRow = cex, scale = "none",
#          ColSideColors = NULL,
          zlim = zlim)
}

# Check consistency
observe({
  # If minimax, no grouping
  # No factor, no summarization possible
  if (is.null(input$factorSelected) || input$factorSelected == "none") {
    if (is.null(input$summarize) || input$summarize != "0") {
      updateSelectInput(session, "summarize", selected = "0")
      updateCheckboxInput(session, "summarize.interleave", value = FALSE) # TODO: not acting
    }
    if (input$scale_group != "none")
      updateSelectInput(session, "scale_group", selected = "none")
    if (input$scale_group != "color.zero")
      updateSelectInput(session, "color.zero", selected = "none")
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
  if (!is.null(input$demo_file) && input$demo_file %in% names(tdms.examples()))
    c(input$demo_file)
  else
    c(sprintf("File is \"%s\"<br>", input$file$name))
})

##### ADJUST TAB ######

# Heatmap view
output$heatmapNorm <- renderPlot({
  # Get the data
  #tdms <- proc_data()
  .heatmap()
})

# Updating the factor selection popup
output$factorSelectInput <- renderUI({
  .pData <- pData()
  factors <- list("None available" = "none")
  if (!is.null(.pData) && ncol(.pData > 1))
    factors = c("None" = "none", colnames(.pData)[-1])
  selectInput("factorSelected", label = strong("Group columns by factor"), choices = factors)
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
  selectInput("color.zero", label = tagList(strong("Zero"), em("line for group")),
              choices = groups)
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
      mean = "Centering to the average.",
      median = "Centering to the median.",
      minmax = "Centering to the minmax.",
      ...  = "Centering undefined") )
  print(length(input$scale))
  .text = c(.text, switch(input$scale,
      none = "Scaling: none, ignore dispersion.",
      sd = "Scaling using standard deviation function.",
      mad = "Scaling using median absolute deviation.",
      minmax = "Scaling to the overall minmax.",
      ...  = "Scaling undefined") )
  .text = c(.text, "Grouping factor: if a factor is selected, columns are grouped and the data are scaled relatively the within group dispersion.")
  .text = c(.text, switch(input$scale_group,
         none = "The global dispersion is used, ignoring groups.",
         all = "The average within group dispersion is used.",
         only.ref = "The dispersion is computed on the reference group alone."))
  # TODO:        excl.ref = "")
  if (input$sort_byfactor) {
    .text = c(.text, "The columns are ordered by the selected factor.")
  }
  .text = c(.text, "Heatmap colors are adjusted according to the Colorize tab.")
  if (is.null(.text)) return("")
  .text = paste(.text, collapse = "</li><li>")
  sprintf("<ul><li>%s</li></ul>", .text)
})

# Information
output$normalize_info <- renderText({
  ""
})


##### COLORIZE TAB #####

# TODO: color scale coefficient, slider?
# TODO: check color scale is finite

# TODO: final rescaling, add a tab?
#       add an histogram and a colorscale below
#       ensure that zero is at the middle of the scale

# TODO: Adjust color: constrast slider or auto-adjust -3%, +3%

# Heatmap view
output$heatmapCalib <- renderPlot({
  # Get the data
  #print(round(tdms$exprs[14:16,], 2))
  #tdms <- calib_data()
  #print(round(tdms$exprs[14:16,], 2))
  .heatmap()
})

# Notice
output$noticeCalib <- reactive({
  .text = "This panel allows colorizing the data."
  .text = c(.text, ifelse(input$zero_center_palette,
                         "The zero value is bind to the center of the color scale.",
                         "The zero value is not bind to the center of the color scale. The color scale is mapped to the robust minimum and maximum.") )
  .text = c(.text, "A percentage of the lowest or highest data are considered out of scale and are removed before computing minimum and maximum of the data to be mapped onto the color scale." )
  .text = c(.text, "A function that changes the dynamic of the data allows to reinforce lowest data (sqrt - squarre root) or the highest data (squarre).")
  if (is.null(.text)) return("")
  .text = paste(.text, collapse = "</li><li>")
  sprintf("<ul><li>%s</li></ul>", .text)
})

output$colorbarCalib <- renderPlot({
  tdms <- calib_data()
  old.mar = par("mar")
  par(mar = c(2, 20, 0.1, 20))
  image(x=seq(round(optZmin(), 2), round(optZmax(), 2), length.out = 254), z = matrix(1:254, nc = 1),
        col = colPal(), xaxp=c(round(optZmin(), 2), round(optZmax(), 2), 4), axes = FALSE, xlab = "")
  axis(1, at = seq(optZmin(), optZmax(), 1), cex.axis = 0.8)
  par(mar = old.mar)
})


##### CLUSTERING TAB ######

# TODO: to do or not to do?
# TODO: cluster rows

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
  #tdms <- post_data()
  .heatmap()
})

# Notice
output$noticePostProc <- reactive({
  .text = NULL
  if (input$summarize!="0") {
    .text = c(.text, "Summarization allows comparing the distribution of values between groups: quantiles are computed and displayed with color. If a shift is observed, it might be meaningful. Summarization is a color representation of boxplots.")
  } else {
    .text = c(.text, "To apply Summarization, a grouping factor must be selected (Input panel).")
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
    if (!is.null(input$demo_file) && input$demo_file %in% names(tdms.examples()))
      .filename = input$demo_file
    else
      .filename = input$file$name
    .filename = gsub('\\..{3,4}$', '', .filename)
    .filename = gsub(' ', '_', .filename)
    paste(.filename, '-', Sys.Date(), '.txt', sep='')
  },
  content = function(file) {
    tdms <- post_data()
    # TODO: round via options parameter
    if (is.null(tdms)) return(NULL)
    write.tdms(tdms, file)
  }
)

##### CONTROL-OPTIONS TAB ######

# Notice
output$noticeOptions <- reactive({
  .text = "This panel allows adjusting various options."
  .text = c(.text, "Various Brewer diverging palettes are available.")
  if (input$scale_adjust == "percent") {
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
               "BGBr" = c("#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E"),
               "gne BuWhRd" = c("#0133EE", "#EEEEEE", "#FE0101"),
               "mev BuBkYl" = c("#0101FE", "#111111", "#FEFE01"),
               "mev GrBkRd" = c("#01FE01", "#111111", "#FE0101")
               )
  colorRampPalette(.cp, space = "Lab")(64)
})

# TODO: add logging

# Limit rows/columns
optMaxRow <- reactive({
  sel = as.numeric(input$data.sel.rows)
  if (is.null(sel) || is.na(sel)) return(40)
  if (sel != 0 && abs(sel)<10) return(40)
  sel
})
optMaxCol <- reactive({
  sel = as.numeric(input$data.sel.cols)
  if (is.null(sel) || is.na(sel) || abs(sel)<10) return(20)
  if (sel != 0 && abs(sel)<10) return(40)
  sel
})


##### HELP TAB ######

# Notice
output$noticeHelp <- reactive({
  .text = c("The main aim is to offer supervised visualization of the data.",
            "The first aim is to rocess the data in order to visualize them adaquately in MeV and in relation to the statistical analyses. A lot of data analyses aim to show difference between groups. Whereas statistical tests transform data in order to relate difference to the dispersion, visualisations do not. Visualizations usually consider data in a discovery way. Once centered, data are scaled to the overall dispersion which encompass the within group dispersion and the between group dispersion. Such a process over-estimates the dispersion when real differences are observed. The scaling is thus too strong and may mask interestering differences.",
            "The second aim is to propose a heatmap alternative to the boxplots. When features under interest become too numerous, boxplots may become difficult to compare. A heatmap gives a more synthetic view. The boxplot limits are replaced by quantiles. Those quantiles replace the indivual data of each group. Quantiles being in the same amount whether a group is small or big, it is easier to compare groups when they are of various size or of large size. Moreover quantiles allow a rapid comparison of groups. Finally, quantiles could be clustered using euclidean distance in order to show some links.",
            "Here, we exploit the TDMS file format and the MeV software. The interface allows to process the data interactively. The user can tune various parameters in order to understand their influence ont the visualisation.",
            "Data to color mapping is also a major point that helps catching the real sense of the data. In particular, binding zero to the center of the diverging color scale makes visualisation readly comprehensive.")
  .text = paste(.text, collapse = "</li><li>")
  sprintf("<ul><li>%s</li></ul>", .text)
})

##### TODO

# TODO: windsoring

})

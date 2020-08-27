library(shiny)
library(shinyjs)

jsfile <- "getdata.js"
cssfile <- "style.css"
js_scroll_file <-'scroll.js'

shinyUI(
    fluidPage(
      titlePanel("MicAff Application!"),
      sidebarLayout(
        sidebarPanel(
          fileInput("read.affymetrix.files", "Choose files for analysis",
                    multiple = TRUE,
                    accept = c(".CEL", ".cel")),
          radioButtons("normalization.algorithm", "Normalization algorithm:",
                       c("MAS-5" = "mas5",
                         "RMA" = "rma")),
          checkboxInput("performing QC-Stats", "Perform QC-Stats?", FALSE),
          checkboxInput("performing.clustering", "Perform clustering?", FALSE),
          checkboxInput("performing.NUSE.RLE", "Perform NUSE and RLE charts?", FALSE),
          downloadButton("downloading.preprocessed.data", "Download preprocessed data")
          ),
        mainPanel(
          plotOutput("distPlot")
        )
      )
    )
)
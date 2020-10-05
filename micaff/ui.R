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
          radioButtons("normalization.algorithm", "Normalization algorithm:",
                       c("MAS-5" = "mas5",
                         "RMA" = "rma")),
          fileInput("read.affymetrix.files", "Choose files for analysis",
                    multiple = TRUE,
                    accept = c(".CEL", ".cel")),
          #checkboxInput("performing QC-Stats", "Perform QC-Stats?", FALSE),
          #checkboxInput("performing.clustering", "Perform clustering?", FALSE),
          #checkboxInput("performing.NUSE.RLE", "Perform NUSE and RLE charts?", FALSE),
          actionButton("calculate.stats", "Calculate statistics and plot charts!"),
          downloadButton("downloading.preprocessed.data", "Download preprocessed data")
          ),
        mainPanel(
          plotOutput("boxplot"),
          plotOutput("qc.stats.plot", width = "100%", height = "100%"),
          #plotOutput("clustering.plot"),
          plotOutput("nuse.plot"),
          plotOutput("rle.plot")
        )
      )
    )
)
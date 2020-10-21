library(shiny)
library(shinyjs)

jsfile <- "getdata.js"
cssfile <- "style.css"
js_scroll_file <-'scroll.js'

shinyUI(
    fluidPage(
      titlePanel("MicAff Application!"),
      
      fluidRow(
        column(3,
               radioButtons("normalization.algorithm", "Normalization algorithm:",
                       c("MAS-5" = "mas5",
                         "RMA" = "rma"))
               ),
        column(3,
               fileInput("read.affymetrix.files", "Choose files for analysis",
                    multiple = TRUE,
                    accept = c(".CEL", ".cel"))
               ),
        column(3,
               selectInput(
            inputId = "input.control.labels", label = "Click white box below to select the controls.",
            choices =  c("first you have to load data!", "--none--"),
            multiple = TRUE, selected = c("first you have to load data!", "--none--"))
               ),
        column(3,
               numericInput("num.genes", "Number of relevant genes", 50, min = 2, max = 150, step = 1, width = NULL)
               ),
        actionButton("calculate.stats", "Calculate statistics and plot charts!"),
        downloadButton("downloading.normalized.data", "Download normalized data"),
        downloadButton("downloading.site", "Download the sites state.")
      ),
      hr(),
 
      tabsetPanel(
        tabPanel("Quality control charts", 
                 plotOutput("boxplot"),
                 plotOutput("qc.stats.plot", width = "100%", height = "100%"),
                 plotOutput("nuse.plot"),
                 plotOutput("rle.plot")
                 ), 
        tabPanel("Checking normalization results", 
                 plotOutput("boxplot.norm"),
                 plotOutput("ma.plot"),
                 #plotOutput("volcano"),
                 #plotOutput("dendrogram"),
                 plotOutput("volcano.moderated"),
                 plotOutput("dentrogram.moderated")
                 )
      )
      # titlePanel("Quality control charts:"),
      # titlePanel("Checking normalization results:"),
    )
)
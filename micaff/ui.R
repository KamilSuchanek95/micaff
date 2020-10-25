library(shiny)
library(shinythemes)


#library(shinyjs)

#jsfile <- "getdata.js"
#cssfile <- "style.css"
#js_scroll_file <-'scroll.js'

shinyUI(
    fluidPage(theme = shinytheme("united"),
              br(),
      titlePanel("MicAff Application!"),
      br(),
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
               )
      ),
      br(),
      fluidRow(
        actionButton("calculate.stats", "Calculate statistics and plot charts!"),
        downloadButton("downloading.normalized.data", "Download normalized data"),
        downloadButton("downloading.site", "Download the sites state.")
      ),
      br(), hr(), br(),
      tabsetPanel(
        tabPanel("Quality control charts", 
                 br(),
                 plotOutput("boxplot"),
                 br(),
                 plotOutput("qc.stats.plot", width = "100%", height = "100%"),
                 br(),
                 plotOutput("nuse.plot"),
                 br(),
                 plotOutput("rle.plot")
                 ), 
        tabPanel("Checking normalization results", 
                 br(),
                 plotOutput("boxplot.norm"),
                 br(),
                 plotOutput("ma.plot", height = "800px")
                 ),
        tabPanel("Statistics summaries", 
                 #plotOutput("volcano"),
                 #plotOutput("dendrogram"),
                 br(),
                 plotOutput("volcano.moderated", height = "800px"),
                 br(),
                 plotOutput("dentrogram.moderated", height = "800px")
                 )
      )
    )
)
library(shiny)
library(shinythemes)

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
        column(4,
               numericInput("num.genes", "Number of relevant genes", 50, min = 2, max = NA, step = 1, width = NULL),
               actionButton("update.statistics.plots", "Update relevant genes"),
               hr(),
               numericInput("p.val.threshold", "Threshold for p value", 0.5, min = 0, max = 1, step = 0.0001, width = NULL),
               numericInput("f.c.threshold", "Threshold for fold change", 2, min = NA, max = NA, step = 0.1, width = NULL),
               actionButton("update.statistics.for.thresholds", "Update relevant genes for thresholds")
               )
      ),
      br(),
      fluidRow(
        actionButton("calculate.stats", "Calculate statistics and plot charts"),
        downloadButton("downloading.normalized.data", "Download normalized data"),
        downloadButton("downloading.pvals", "Download statistics")
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
                 br(),
                 plotOutput("volcano.moderated", height = "800px"),
                 br(),
                 plotOutput("dentrogram.moderated", height = "800px")
                 ),
        tabPanel("Table of relevant genes",
                 br(),
                 tableOutput('table.relevant')
                 )
      )
    )
)
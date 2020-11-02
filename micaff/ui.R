library(shiny)
library(shinythemes)
library(DT)

shinyUI(
    fluidPage(theme = shinytheme("united"),
              tags$style(HTML("
              .tabbable > .nav > li > a[data-value='Statistics summaries'] {background-color: orange; color:black}
              .tabbable > .nav > li > a[data-value='Table of relevant genes'] {background-color: orange; color:black}
                              ")),
              br(),
      titlePanel("MicAff Application!"),
      br(),
      fluidRow(
        column(3,
               radioButtons("normalization.algorithm", "Normalization algorithm:",
                       c("MAS-5" = "mas5",
                         "RMA" = "rma")),
               br(),
               fileInput("read.affymetrix.files", "Choose CEL files for analysis",
                         multiple = TRUE,
                         accept = c(".CEL", ".cel")),
               br(),
               checkboxInput("do.qc", "Perform quality control", value = TRUE, width = NULL),
               br(),
               checkboxInput("do.check", "Check the normalization results", value = TRUE, width = NULL)
               ),
        column(3,
               selectInput(
            inputId = "input.control.labels", label = "Click white box below to select the controls.",
            choices =  c("first you have to load data!", "--none--"),
            multiple = TRUE, selected = c("first you have to load data!", "--none--"))
               ),
        column(4,
               numericInput("num.genes", "Number of relevant genes", 50, min = 2, max = NA, step = 1, width = NULL),
               fluidRow(
                 column(7,actionButton("update.statistics.plots", "Update relevant genes")),
                 column(1,checkboxInput("FDR.number", "FDR", value = FALSE, width = NULL))
               ),
               hr(),
               fluidRow(
                 column(10, numericInput("p.val.threshold", "Threshold for p value", 0.05, min = 0, max = 1, step = 0.0001, width = NULL)),
                 column(2, checkboxInput("FDR", "FDR", value = FALSE, width = NULL))
               ),
               numericInput("f.c.threshold", "Threshold for fold change", 2, min = NA, max = NA, step = 0.1, width = NULL),
               actionButton("update.statistics.for.thresholds", "Update relevant genes for thresholds")
               )
      ),
      br(),
      fluidRow(
        column(4,actionButton("calculate.stats", "Calculate statistics and plot charts")),
        column(4,downloadButton("downloading.normalized.data", "Download normalized data")),
        column(3,downloadButton("downloading.pvals", "Download statistics"))
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
        tabPanel("Tabs for a specific number of genes",
                 br(),
                 uiOutput('title.specify'),
                 br(),
                 tabsetPanel(
                   tabPanel("Statistics summaries",
                            br(),
                            plotOutput("volcano.moderated", height = "800px"),
                            br(),
                            plotOutput("dendrogram.moderated", height = "800px")
                   ),
                   tabPanel("Table of relevant genes",
                            br(),
                            DT::dataTableOutput('table.relevant')
                   )
                 )
                 ),
        tabPanel("Tabs for genes within thresholds",
                 br(),
                 uiOutput('title.threshold'),
                 br(),
                 tabsetPanel(
                   tabPanel("Statistics summaries",
                            br(),
                            plotOutput("volcano.moderated.threshold", height = "800px"),
                            br(),
                            plotOutput("dendrogram.moderated.threshold", height = "800px")
                   ),
                   tabPanel("Table of relevant genes",
                            br(),
                            DT::dataTableOutput('table.threshold.p.val')
                   )
                 )
        ),
        tabPanel("Table of all genes",
                 br(),
                 DT::dataTableOutput('table.all')
        )
      )
    )
)
# BiocManager::install("affy")
library(affy)
# BiocManager::install("simpleaffy")
library(simpleaffy)

options(shiny.maxRequestSize=30000*1024^10)
data <- 0
exprs.mas5 <- 0
exprs.rma <- 0

shinyServer(function(input, output) {
  
  observeEvent(input$read.affymetrix.files, {
    name <- input$read.affymetrix.files$name
    datapath <- input$read.affymetrix.files$datapath
    
    data <- affy::read.affybatch(datapath)
  
    norm.alg <- switch(input$normalization.algorithm,
                        mas5 = "mas5",
                        rma = "rma")
    print(norm.alg)
    })
})

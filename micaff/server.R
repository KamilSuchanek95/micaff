# BiocManager::install("affy")
library(affy)
# BiocManager::install("simpleaffy")
library(simpleaffy)

options(shiny.maxRequestSize=30000*1024^10)
data = 0
normali

shinyServer(function(input, output) {
  input$
  observeEvent(input$read.affymetrix.files, {
    name <- input$read.affymetrix.files$name
    datapath <- input$read.affymetrix.files$datapath
    
    data <<- affy::read.affybatch(datapath)
  
    normalization.algorithm <<- switch(input$normalization.algorithm,
                                       mas5 = "mas5",
                                       rma = "rma")
    print(normalization.algorithm)
    })
  
})

library(affy)

options(shiny.maxRequestSize=30000*1024^2)


shinyServer(function(input, output) {
    
  observeEvent(input$file1, {
    name <- input$file1$name
    datapath <- input$file1$datapath
    
    data <- affy::read.affybatch(datapath)
    
    d_mas5 = affy::mas5(data)
    #  tworzenie macierzy ekpresji
    d_exprs = affy::exprs(d_mas5)
    d_log2 = log(d_exprs, 2)
    
    e <- 4
    
    })
})

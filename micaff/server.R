# BiocManager::install("affy")
library(affy)
# BiocManager::install("simpleaffy")
library(simpleaffy)
library(latticeExtra)
# BiocManager::install("affyPLM")
library(affyPLM)
# BiocManager::install("limma")
library(limma)

options(shiny.maxRequestSize=30000*1024^10)
data <- 0
data.mas5 <- 0
data.rma <- 0

display.report <- function(norm.alg, data, output){
  if(norm.alg == "mas5"){
    data.norm = simpleaffy::call.exprs(data, "mas5")
  } else{
    data.norm = simpleaffy::call.exprs(data, "rma")
  }
  ###
  qc = simpleaffy::qc(data, data.norm)
  output$qc.stats.plot <- renderPlot({
    simpleaffy::plot(qc)})
  ### 
  dataPLM = fitPLM(data)
  output$nuse.plot <- renderPlot({ 
    boxplot(dataPLM, main = "NUSE", outline = FALSE, col="lightblue", 
            las=3, whisklty=0, staplelty=0)})
  output$rle.plot <- renderPlot({
    Mbox(dataPLM, main = "RLE", outline = FALSE, col="mistyrose", 
         las=3, whisklty=0, staplelty=0)})
  ### 
  return(data.norm)
}

shinyServer(function(input, output) {
  
  observeEvent(input$read.affymetrix.files, {
    
    name <- input$read.affymetrix.files$name
    datapath <- input$read.affymetrix.files$datapath
    
    data <- affy::read.affybatch(datapath)
    
    norm.alg <- switch(input$normalization.algorithm,
                       mas5 = "mas5",
                       rma = "rma")
    if(norm.alg == "mas5"){
      data.mas5 = simpleaffy::call.exprs(data, "mas5")
      exprs.mas5 = exprs(data.mas5)
      ###
      qc = simpleaffy::qc(data, data.mas5)
      output$qc.stats.plot <- renderPlot({
        simpleaffy::plot(qc)})
      ##
      dd = dist2(log2(exprs(data)))
      dd.row <- as.dendrogram(hclust(as.dist(dd)))
      row.ord <- order.dendrogram(dd.row)
      legend = list(top=list(fun=latticeExtra::dendrogramGrob,
                             args=list(x=dd.row, side="top")))
      output$clustering.plot <- renderPlot({
        levelplot(dd[row.ord, row.ord],
                  scales=list(x=list(rot=90)), xlab="",
                  ylab="", legend=legend)})
      ##
      dataPLM = fitPLM(data)
      output$nuse.plot <- renderPlot({ 
        graphics::boxplot(dataPLM, main="NUSE", outline = FALSE, col="lightblue", 
                las=3, whisklty=0, staplelty=0)})
      output$rle.plot <- renderPlot({
        Mbox(dataPLM, main="RLE", outline = FALSE, col="mistyrose", 
             las=3, whisklty=0, staplelty=0)})
    } else{
      data.rma = simpleaffy::call.exprs(data, "rma")
      exprs.rma = exprs(data.rma)
      ###
    }
  })
})

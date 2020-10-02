# BiocManager::install("affy")
library(affy)
library("genefilter")
library("RColorBrewer")
# BiocManager::install("simpleaffy")
library(simpleaffy)
library(latticeExtra)
# BiocManager::install("affyPLM")
library(affyPLM)
# BiocManager::install("limma")
library(limma)

My_NUSE_data <- function(x,type=c("plot","values","stats","density"),ylim=c(0.9,1.2),...){
  
  compute.nuse <- function(which){
    nuse <- apply(x@weights[which,],2,sum)
    1/sqrt(nuse)
  }
  
  type <- match.arg(type)
  model <- x@model.description$modelsettings$model
  ## if (type == "values" || type == "stats" || type == "density"){
  
  if (x@model.description$R.model$which.parameter.types[3] == 1 & x@model.description$R.model$which.parameter.types[1] == 0 ){
    grp.rma.se1.median <- apply(se(x), 1,median,na.rm=TRUE)
    grp.rma.rel.se1.mtx <- sweep(se(x),1,grp.rma.se1.median,FUN='/')
  } else {
    # not the default model try constructing them using weights.
    which <-indexProbesProcessed(x)
    ses <- matrix(0,length(which) ,4)
    
    for (i in 1:length(which))
      ses[i,] <- compute.nuse(which[[i]])
    
    
    grp.rma.se1.median <- apply(ses, 1,median)
    grp.rma.rel.se1.mtx <- sweep(ses,1,grp.rma.se1.median,FUN='/')
  }
  if (type == "values"){
    return(grp.rma.rel.se1.mtx)
  } else if (type == "density"){
    plotDensity(grp.rma.rel.se1.mtx,xlim=ylim,...)
  } else if (type=="stats"){
    Medians <- apply(grp.rma.rel.se1.mtx,2,median)
    Quantiles <- apply(grp.rma.rel.se1.mtx,2,quantile,prob=c(0.25,0.75))
    nuse.stats <- rbind(Medians,Quantiles[2,] - Quantiles[1,], Quantiles[1,], Quantiles[2,])
    rownames(nuse.stats) <- c("median","IQR","25%", "75%")
    return(nuse.stats)
  }
  if (type == "plot"){	
    boxplot(data.frame(grp.rma.rel.se1.mtx),ylim=ylim,range=0,...)
  }
}
My_RLE_Plot <- function(dataPLM) {
  par(mar = c(7,5,2,2))
  brewer.cols <- brewer.pal(num.probes, "Set1")
  medianchip <- apply(coefs(dataPLM), 1, median)
  M <- sweep(coefs(dataPLM),1,medianchip,FUN='-')
  my.quantile = apply(M,2,quantile)
  my.mini = min(my.quantile["25%",])
  my.maxi = max(my.quantile["75%",])
  y.lim = c(my.mini, my.maxi)
  graphics::boxplot(M, main="RLE Plot", col = brewer.cols,
                    outline = FALSE, ylim = y.lim, 
                    las=3, whisklty=0, staplelty=0)
  grid()
  abline(0,0)
}
My_NUSE_Plot <- function(dataPLM) {
  par(mar = c(7,5,2,2))
  brewer.cols <- brewer.pal(num.probes, "Set1")
  nuse.data = My_NUSE_data(dataPLM, type = "stats")
  my.mini = min(nuse.data["25%",])
  my.maxi = max(nuse.data["75%",])
  y.lim = c(my.mini, my.maxi)
  NUSE(dataPLM, main="NUSE", ylim = y.lim,
       outline = FALSE, las=3, whisklty=3, staplelty=0,
       col = brewer.cols)
  grid()
}
My_Box_Plot <- function() {
  par(mar = c(7,5,2,2))
  brewer.cols <- brewer.pal(num.probes, "Set1")
  BiocGenerics::boxplot(data, col = brewer.cols, las = 3,
                        ylab = "Unprocessed log (base 2)scale Probe Intensities")
  grid()
}
My_check_and_normalise_data <- function(norm.alg, data){
  if(norm.alg == "mas5"){
    data.norm = simpleaffy::call.exprs(data, "mas5")
  } else{
    data.norm = simpleaffy::call.exprs(data, "rma")
  }
  return(data.norm)
}
options(shiny.maxRequestSize=30000*1024^10)
data <- 0
data.mas5 <- 0
data.rma <- 0
num.probes <- 0

display.report <- function(norm.alg, data, output){
  
  data.norm <- My_check_and_normalise_data(norm.alg, data)
  ###
  output$boxplot <- renderPlot({
    My_Box_Plot()
  })
  ###
  qc = simpleaffy::qc(data)
  output$qc.stats.plot <- renderPlot({
    simpleaffy::plot.qc.stats(qc)})
  ###
  rrr = {
    #dd = dist2(log2(exprs(data)))
    #dd.row <- as.dendrogram(hclust(as.dist(dd)))
    #row.ord <- order.dendrogram(dd.row)
    #legend = list(top=list(fun=latticeExtra::dendrogramGrob,
    #                       args=list(x=dd.row, side="top")))
    #output$clustering.plot <- renderPlot({
    #  levelplot(dd[row.ord, row.ord],
    #            scales=list(x=list(rot=90)), xlab="",
    #            ylab="", legend=legend)})
  }
  ### 
  dataPLM = fitPLM(data)
  output$nuse.plot <- renderPlot({
    My_NUSE_Plot(dataPLM)
  })
  output$rle.plot <- renderPlot({
    My_RLE_Plot(dataPLM)})
  ### 
  return(data.norm)
}

shinyServer(function(input, output) {
  
  observeEvent(input$read.affymetrix.files, {
    
    name <- input$read.affymetrix.files$name
    datapath <- input$read.affymetrix.files$datapath
    
    data <- affy::read.affybatch(datapath)
    sampleNames(data) = name
    sampleNames(data) = sub("\\.CEL$", "", sampleNames(data))
    num.probes = length(sampleNames(data))
    
    norm.alg <- switch(input$normalization.algorithm,
                       mas5 = "mas5",
                       rma = "rma")
    
    data.norm = display.report(data = data, norm.alg = norm.alg, output = output)
    
    output$downloading.preprocessed.data <- downloadHandler(
      filename =  function() {paste("data_", norm.alg,".txt", sep = "")},
      content = function(file) {write.table(exprs(data.norm), file = file)}
    )
    
    
    
  })
})

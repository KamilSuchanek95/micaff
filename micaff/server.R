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
library(shiny)
library("pheatmap")


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
  par(mar = c(10,5,2,2))
  brewer.cols <- brewer.pal(num.probes, "Set1")
  medianchip <- apply(coefs(dataPLM), 1, median)
  M <- sweep(coefs(dataPLM),1,medianchip,FUN='-')
  my.quantile = apply(M,2,quantile)
  my.mini = min(my.quantile["25%",])
  my.maxi = max(my.quantile["75%",])
  y.lim = c(my.mini, my.maxi)
  graphics::boxplot(M, main="RLE Plot", col = brewer.cols,
                    outline = FALSE, ylim = y.lim, 
                    las=3, whisklty=3, staplelty=0)
  grid()
  abline(0,0)
}
My_NUSE_Plot <- function(dataPLM) {
  par(mar = c(10,5,2,2))
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
My_Box_Plot <- function(data, num.probes, ylab) {
  par(mar = c(10,5,2,2))
  brewer.cols <- brewer.pal(num.probes, "Set1")
  BiocGenerics::boxplot(data, col = brewer.cols, las = 3,
                        ylab = ylab)
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
My_MA_Plot <- function(data.norm, control){
  exprs.eset <- exprs(data.norm)
  
  control.array = match(control, sampleNames(data.norm))
  Difference <- rowMeans(exprs.eset[,-control.array]) - rowMeans(exprs.eset[,control.array])
  Average <- rowMeans(exprs.eset)
  plot(Average, Difference)
  grid()
  lines(lowess(Average, Difference),
        col = 'red', lwd = 4)
  abline(h = -2)
  abline(h = 2)
}
My_volcano_moderated <- function(data.norm, control, number.of.relevant.genes){
  exprs.eset <- exprs(data.norm)
  control.array = match(control, sampleNames(data.norm))
  Difference <- rowMeans(exprs.eset[,-control.array]) - rowMeans(exprs.eset[,control.array])
  control.idx = which(sampleNames(data.norm) %in% control)
  fac = c()
  fac[1:num.probes] = "Test"
  fac[control.idx] = "Control"
  population.groups <- factor(fac)
  design <- model.matrix(~population.groups)
  fit <- lmFit (data.norm, design)
  fit.eBayes <- eBayes (fit)
  
  lodd <- -log10(fit.eBayes$p.value[,2])
  o1 <- order(abs(Difference), decreasing = TRUE)[1:number.of.relevant.genes]
  oo2 <- order(abs (fit.eBayes$t[,2]), decreasing = TRUE)[1:number.of.relevant.genes]
  oo <- union (o1, oo2)
  ii <- intersect (o1, oo2)
  plot (Difference[-oo], lodd[-oo],
        cex = .25, xlim = c (-3,3),
        ylim = range(lodd), xlab = 'Average (log) Fold-change',
        ylab = 'LOD score – Negative log10 of P-value')
  points(Difference [oo2], lodd [oo2],
         pch = 5, col ='red', cex = 1,
         lwd = 1)
  points (Difference[o1], lodd[o1],
          pch = 18, col = 'blue', cex = 1, lwd = 1)
  
  title("Volcano Plot with moderated t-statistics")
  grid()
  abline(v = c(-1,1), lwd = 0.5)
  abline(h= c(-log10(0.05),-log10(0.01),-log10(0.001)), 
         col=c('red','blue','black'), lwd = 0.5)
  #legend("topleft", , pch = c(5, 18), col = c('red', 'blue'), cex = 0.75, bg="transparent",
  #       legend = c("First 50 with the lowest p value","First 50 with the highest Fold-change value"))
  My_list <- list("fit.eBayes" = fit.eBayes, "ii" = ii)
  return(My_list)
}
My_heatmap <- function(data.norm, num.probes, control, fit.eBayes_ii){
  par(mar = c(10,5,2,2))
  ii = fit.eBayes_ii$ii
  exprs.eset = exprs(data.norm)
  ii.mat <- exprs.eset[ii,]
  # ii.df <- data.frame(ii.mat)
  # brewer.cols <- brewer.pal(num.probes, "Set1")
  # hmcol <- colorRampPalette(brewer.pal(num.probes, 'Greys'))(256)
  # spcol = c()
  # spcol[1:num.probes] = 'grey10'
  # spcol[which(dimnames(ii.mat)[[2]] %in% control)] = 'grey80'
  # heatmap (ii.mat, col = hmcol, ColSideColors = spcol)#,margins = c (10,15))
  
  where_control = which(sampleNames(data.norm) %in% control)
  col_groups = c()
  col_groups[1:num.probes] = "Test"
  col_groups[where_control] = "Control"
  mat_col = data.frame(group = col_groups)
  # mat_colors = list(group = brewer.pal(2,"Set1"))
  # mat_colors$group = mat_colors$group[1:2]
  rownames(mat_col) = colnames(data.norm)
  pheatmap(mat = abs(ii.mat), annotation_col = mat_col)# , clustering_method = "complete")
}
display.report <- function(data, data.norm, input, output, num.probes){

  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Creating report", value = 0)
  n <- 8
  progress$inc(1/n, detail = "boxplot of unnormalized data")
  
  output$boxplot <- renderPlot({
    My_Box_Plot(data = data, num.probes = num.probes, ylab = "Unprocessed log (base 2)scale Probe Intensities")
  })
  ###
  progress$inc(1/n, detail = "calculating and plot qc report")
  par()
  qc = simpleaffy::qc(data)
  output$qc.stats.plot <- renderImage({
    outfile <- tempfile(fileext = '.png')
    width2 <- 100
    if(num.probes < 20){width2 <- 200}
    png(outfile, width = 480, height = ceiling(num.probes/4) * width2)
    simpleaffy::plot.qc.stats(qc)
    dev.off()
    list(src = outfile,
         alt = "This is alternate text")
    }, deleteFile = TRUE)
  ### 
  progress$inc(2/n, detail = "calculating PLM and plot NUSE and RLE")
  dataPLM = fitPLM(data)
  output$nuse.plot <- renderPlot({
    My_NUSE_Plot(dataPLM)
  })
  output$rle.plot <- renderPlot({
    My_RLE_Plot(dataPLM)})
  ### 
  progress$inc(1/n, detail = "boxplot of normalized data")
  output$boxplot.norm <- renderPlot({
    My_Box_Plot(data = data.norm, num.probes = num.probes, ylab = "Normalized log (base 2)scale Probe Intensities")
  })
  ###
  progress$inc(1/n, detail = "MA plot")
  control = input$input.control.labels
  output$ma.plot <- renderPlot({
    My_MA_Plot(data.norm = data.norm, control = control)
  })
  ###
  progress$inc(1/n, detail = "Volcano-plot with moderated statistics")
  number.of.relevant.genes = input$num.genes
  output$volcano.moderated <- renderPlot({
    fit.eBayes_ii <<- My_volcano_moderated(data.norm, control = control, number.of.relevant.genes = number.of.relevant.genes)
  })
  ###
  progress$inc(1/n, detail = "heatmap plot")
  output$dentrogram.moderated <- renderPlot({
    My_heatmap(data.norm = data.norm, num.probes = num.probes, control = control, fit.eBayes_ii = fit.eBayes_ii)
  })
  #write.table(fit.eBayes$t,"My_fit_eBayes") 
}

options(shiny.maxRequestSize=30000*1024^10)
.GlobalEnv$data <- 0
.GlobalEnv$data.norm <- 0
.GlobalEnv$num.probes <- 0
.GlobalEnv$fit.eBayes_ii <- 0

shinyServer(function(input, output, session) {
  
  observeEvent(input$calculate.stats, {
    display.report(data = data, data.norm = data.norm, input = input, output = output,
                   num.probes = num.probes)
  })
  
  observeEvent(input$read.affymetrix.files, {
    
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Reading and normalizing data", value = 0)
    n <- 4
    progress$inc(1/n, detail = "reading data")
    name <- input$read.affymetrix.files$name
    datapath <- input$read.affymetrix.files$datapath
    data <<- affy::read.affybatch(datapath)
    sampleNames(data) <<- name
    sampleNames(data) <<- sub("\\.CEL$", "", sampleNames(data))
    num.probes <<- length(sampleNames(data))
    
    progress$inc(1/n, detail = "updating sample names")
    updateSelectInput(session = session, inputId = "input.control.labels",
                      choices = c(sampleNames(data)))
    
    norm.alg <- switch(input$normalization.algorithm,
                       mas5 = "mas5",
                       rma = "rma")
    progress$inc(1/n, detail = "normalizing data")
    data.norm <<- My_check_and_normalise_data(norm.alg, data)
    
    progress$inc(1/n, detail = "preparing to share normalized data")
    output$downloading.normalized.data <- downloadHandler(
      filename =  function() {paste("data_", norm.alg,".txt", sep = "")},
      content = function(file) {write.table(exprs(data.norm), file = file)}
    )
    
    #output$downloading.site <- downloadHandler(
    #)
  })
})

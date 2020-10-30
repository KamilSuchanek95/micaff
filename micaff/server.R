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

options(shiny.maxRequestSize=30000*1024^10)
.GlobalEnv$data <- 0
.GlobalEnv$data.norm <- 0
.GlobalEnv$num.probes <- 0
.GlobalEnv$fit.eBayes_ii <- list()
.GlobalEnv$norm.alg <- 0
.GlobalEnv$tab <- 0

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
My_Box_Plot <- function(data, num.probes, ylab, main) {
  par(mar = c(10,5,2,2))
  brewer.cols <- brewer.pal(num.probes, "Set1")
  BiocGenerics::boxplot(data, col = brewer.cols, las = 3,
                        ylab = ylab, main = main)
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
  plot(Average, Difference, ylab = "M - log2 fold change", xlab = "A - log2 mean expression", main = "MA plot")
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
  text(Difference[ii], lodd[ii], row.names(data.norm)[ii], pos = 3)
  title("Volcano Plot with moderated t-statistics")
  grid()
  abline(v = c(-2,-1,1,2), lwd = 0.5)
  abline(h= c(-log10(0.05),-log10(0.01),-log10(0.001)), 
         col=c('red','blue','black'), lwd = 0.5)
  #legend("topleft", , pch = c(5, 18), col = c('red', 'blue'), cex = 0.75, bg="transparent",
  #       legend = c("First 50 with the lowest p value","First 50 with the highest Fold-change value"))
  My_list <- list("fit.eBayes" = fit.eBayes, "ii" = ii)
  return(My_list)
}
My_volcano_moderated_threshold <- function(data.norm, control, my_index){
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

  plot (Difference[-my_index], lodd[-my_index],
        cex = .25, xlim = c (-3,3),
        ylim = range(lodd), xlab = 'Average (log) Fold-change',
        ylab = 'LOD score – Negative log10 of P-value')
  points(Difference [my_index], lodd [my_index],
         pch = 5, col ='red', cex = 1,
         lwd = 1)
  points (Difference[my_index], lodd[my_index],
          pch = 18, col = 'blue', cex = 1, lwd = 1)
  text(Difference[my_index], lodd[my_index], row.names(data.norm)[my_index], pos = 3)
  title("Volcano Plot with moderated t-statistics")
  grid()
  abline_h <- input$p.val.threshold
  abline_v <- input$f.c.threshold
  abline(v = c(-abline_v,abline_v), lwd = 0.5, col = c('red'))
  abline(h = -log10(abline_h), 
         col=c('red'), lwd = 0.5)
}
My_heatmap <- function(data.norm, num.probes, control, my_intersect){
  par(mar = c(10,5,2,2))
  ii = my_intersect
  exprs.eset = exprs(data.norm)
  ii.mat <- exprs.eset[ii,]
  
  where_control = which(sampleNames(data.norm) %in% control)
  col_groups = c()
  col_groups[1:num.probes] = "Test"
  col_groups[where_control] = "Control"
  mat_col = data.frame(group = col_groups)
  
  rownames(mat_col) = colnames(data.norm)
  pheatmap(mat = abs(ii.mat), annotation_col = mat_col, main = "Heatmap")# , clustering_method = "complete")
  }
My_FDR_select <- function(input){
  if(input$FDR)
  {
    my_column = "adj.P.Val"
  }
  else
  {
    my_column = "P.Value"
  }
  return(my_column)
}
display.report <- function(data, data.norm, input, output, num.probes, progress, n, norm.alg){
  progress$inc(1/n, detail = "boxplot of unnormalized data")
  output$boxplot <- renderPlot({
    My_Box_Plot(data = data, num.probes = num.probes, ylab = "Unprocessed log (base 2)scale Probe Intensities", main = "Boxplot of unnormalized data")
  })
  ###
  progress$inc(1/n, detail = "calculating and plot qc report")
  par()
  if(FALSE){qc = simpleaffy::qc(data)
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
    My_RLE_Plot(dataPLM)})}
  ### 
  progress$inc(1/n, detail = "boxplot of normalized data")
  output$boxplot.norm <- renderPlot({
    My_Box_Plot(data = data.norm, num.probes = num.probes, ylab = "Normalized log (base 2)scale Probe Intensities", main = "Boxplot of normalized data")
  })
  ###
  progress$inc(1/n, detail = "MA plot")
  control = input$input.control.labels
  output$ma.plot <- renderPlot({
    My_MA_Plot(data.norm = data.norm, control = control)
  })
  ###
  output$title.specify <- renderUI({
    h3("Statistics and charts for specyfic number of genes")
  })
  
  progress$inc(1/n, detail = "Volcano-plot with moderated statistics")
  number.of.relevant.genes = input$num.genes
  fit.eBayes_ii <<- My_volcano_moderated(data.norm, control = control, number.of.relevant.genes = number.of.relevant.genes)
  output$volcano.moderated <- renderPlot({
    My_volcano_moderated(data.norm, control = control, number.of.relevant.genes = number.of.relevant.genes)
  })
  ###
  progress$inc(1/n, detail = "heatmap plot")
  output$dendrogram.moderated <- renderPlot({
    My_heatmap(data.norm = data.norm, num.probes = num.probes, control = control, my_intersect = fit.eBayes_ii$ii)
  })
  ###
  tab <<- topTable(fit.eBayes_ii$fit.eBayes, coef = 2, adjust.method = "BH", 
                   number = length(exprs(data.norm)), sort.by = "none")
  
  progress$inc(1/n, detail = "preparing to share statistics")
  output$downloading.pvals <- downloadHandler(
    filename =  function() {paste("p-vals_", norm.alg,".txt", sep = "")},
    content = function(file) {write.table(tab, file = file)}
  )
  ###
  output$table.relevant <- DT::renderDataTable({tab[fit.eBayes_ii$ii,]})
  ###
  output$table.all <- DT::renderDataTable({tab})
  ###
  # tutaj trzeba zrobic volcano i heatmap warunki oraz wybieranie po FDR lub nie
  my_column = My_FDR_select(input = input)
  my_index_p = base::which(tab[,c(my_column)] < input$p.val.threshold)
  my_index_fc = base::which(abs(tab[,"logFC"]) > input$f.c.threshold)
  my_index = intersect(my_index_p, my_index_fc)
  my_size = length(my_index)
  if(my_size < 1){
    output$title.threshold <- renderUI({
      h3("None of the genes are below the given threshold!")
    })
  }else if(my_size == 1){
    output$title.threshold <- renderUI({
      h3(paste("Only one of the genes are below the given threshold: ",row.names(tab[my_index,])))
    })
    output$table.threshold.p.val <- DT::renderDataTable({tab[my_index,]})
  }else{
    output$title.threshold <- renderUI({
      h3("Statistics and charts for genes within thresholds")
      })
    output$volcano.moderated.threshold <- renderPlot({
        My_volcano_moderated_threshold(data.norm, control = control, my_index = my_index)
      })
    output$dendrogram.moderated.threshold <- renderPlot({
      My_heatmap(data.norm = data.norm, num.probes = num.probes, control = control, my_intersect = my_index)
    })
    output$table.threshold.p.val <- DT::renderDataTable({tab[my_index,]})
  }
}

shinyServer(function(input, output, session) {
  
  observeEvent(input$calculate.stats, {
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Creating report", value = 0)
    n <- 10
    display.report(data = data, data.norm = data.norm, input = input, output = output,
                   num.probes = num.probes, progress = progress, n = n, norm.alg = norm.alg)
    progress$inc(1/n, detail = "rendering charts")
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
    
    norm.alg <<- switch(input$normalization.algorithm,
                       mas5 = "mas5",
                       rma = "rma")
    progress$inc(1/n, detail = "normalizing data")
    data.norm <<- My_check_and_normalise_data(norm.alg, data)
    
    progress$inc(1/n, detail = "preparing to share normalized data")
    output$downloading.normalized.data <- downloadHandler(
      filename =  function() {paste("data_", norm.alg,".txt", sep = "")},
      content = function(file) {write.table(exprs(data.norm), file = file)}
    )
    
  })
  
  observeEvent(input$update.statistics.plots, {
    number.of.relevant.genes = input$num.genes
    control = input$input.control.labels
    output$volcano.moderated <- renderPlot({
      fit.eBayes_ii <<- My_volcano_moderated(data.norm, control = control, number.of.relevant.genes = number.of.relevant.genes)
    })
    ###
    output$dendrogram.moderated <- renderPlot({
      My_heatmap(data.norm = data.norm, num.probes = num.probes, control = control, my_intersect = fit.eBayes_ii$ii)
    })
    output$table.relevant <- DT::renderDataTable({tab[fit.eBayes_ii$ii,]})
  })
  
  observeEvent(input$update.statistics.for.thresholds,{
    control = input$input.control.labels
    my_column <<- ""

    my_column = My_FDR_select(input = input)
    
    my_index_p = which(tab[,c(my_column)] < input$p.val.threshold)
    my_index_fc = which(abs(tab[,"logFC"]) > input$f.c.threshold)
    my_index = intersect(my_index_p, my_index_fc)
    my_size = length(my_index)
    if(my_size < 1){
      output$title.threshold <- renderUI({
        h3("None of the genes are below the given threshold!")
      })
    }else if(my_size == 1){
      output$title.threshold <- renderUI({
        h3(paste("Only one of the genes are below the given threshold: ",row.names(tab[my_index,])))
      })
      output$table.threshold.p.val <- DT::renderDataTable({tab[my_index,]})
    }else{
      output$title.threshold <- renderUI({
        h3("Statistics and charts for genes within thresholds")
      })
      output$volcano.moderated.threshold <- renderPlot({
        My_volcano_moderated_threshold(data.norm, control = control, my_index = my_index)
      })
      output$dendrogram.moderated.threshold <- renderPlot({
        My_heatmap(data.norm = data.norm, num.probes = num.probes, control = control, my_intersect = my_index)
      })
      output$table.threshold.p.val <- DT::renderDataTable({tab[my_index,]})
    }
  })
  
})

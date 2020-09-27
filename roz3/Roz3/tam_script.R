library("genefilter")
library("RColorBrewer")
library("affy")
library("simpleaffy")
library("affyPLM")

# functions 
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
# Wczytanie danych
cel = affy::ReadAffy(celfile.path = "/home/kamil/Pulpit/magisterka/dane tam/cel files")

# Dostosowanie nazw prób
sampleNames(cel) = sub("\\.CEL$", "", sampleNames(cel))
# liczba prób
num.probes = length(sampleNames(cel))

# Normalizacja danych
eset.mas5 = simpleaffy::call.exprs(cel, "mas5")
eset.rma = simpleaffy::call.exprs(cel, "rma")

# Wykres pudełkowy
par(mar = c(7,5,2,2))
brewer.cols <- brewer.pal(num.probes, "Set1")
BiocGenerics::boxplot(cel, col = brewer.cols, las = 3,
  ylab = "Unprocessed log (base 2)scale Probe Intensities")
grid()

# Histogramy gęstości
hist(cel, col = brewer.cols, lty = 1,
  xlab = "Log (base 2) Intensities", lwd = 3)
samp.leg.names <- sampleNames(cel)
legend (10,.4, legend = samp.leg.names,
  lty = 1,
  col = brewer.cols, lwd = 3)

# Wykres metryk Affymetrix
qc = simpleaffy::qc(cel, eset.mas5)
journalpng(file = "brain_Stats_QC_3_5.png", height = ceiling(num.probes/4))
simpleaffy::plot.qc.stats(qc)
dev.off()
journalpng(file = "brain_Stats_QC_3_M.png", height = ceiling(num.probes/4))
simpleaffy::plot.qc.stats(qc, usemid = T)
dev.off()

# Degradacja RNA
RNAdeg <- AffyRNAdeg(cel)
plotAffyRNAdeg(RNAdeg)

# RLE i NUSE
dataPLM = fitPLM(cel)#, background.method = "MAS")
My_RLE_Plot(dataPLM = dataPLM)
My_NUSE_Plot(dataPLM)

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

###

exprs.eset <- exprs(eset.mas5)

exprs.eset.df <- data.frame(exprs.eset)
boxplot(exprs.eset.df,col = brewer.cols,las=3)

controls = c("GSM439779", "GSM439780", "GSM439782", "GSM439783", "GSM439787",
             "GSM439788", "GSM439789", "GSM439792", "GSM439793", "GSM439797",
             "GSM439799", "GSM439802", "GSM439804", "GSM439808", "GSM439810",
             "GSM439812", "GSM439814", "GSM439816", "GSM439818", "GSM439819",
             "GSM439822", "GSM439825", "GSM439826")

library ('limma')

control.array = match(controls, sampleNames(data))
Difference <- rowMeans(exprs.eset[,-control.array]) - rowMeans(exprs.eset[,control.array])
Average <- rowMeans(exprs.eset)
plot(Average, Difference)
lines(lowess(Average, Difference),
      col = 'red', lwd = 4)
abline(h = -2)
abline(h = 2)

control.array = match(controls, sampleNames(data))
Difference <- rowMeans(exprs.eset[,-control.array]) - rowMeans(exprs.eset[,control.array])
control.idx = which(sampleNames(cel) %in% controls)
fac = c()
fac[1:num.probes] = "Test"
fac[control.idx] = "Control"
population.groups <- factor(fac)
design <- model.matrix(~population.groups)
fit <- lmFit (eset.mas5, design)
fit.eBayes <- eBayes (fit)

lodd <- -log10(fit.eBayes$p.value[,2])
o1 <- order(abs(Difference), decreasing = TRUE)[1:50]
oo2 <- order(abs (fit.eBayes$t[,2]), decreasing = TRUE)[1:50]
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



library("genefilter")
library("RColorBrewer")
library("affy")
library("simpleaffy")
library("affyPLM")
library("pheatmap")
library(limma)

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
data = affy::ReadAffy(celfile.path = "/home/kamil/Pulpit/CEL_files_rename")
# Dostosowanie nazw prób
sampleNames(data) = sub("\\.CEL$", "", sampleNames(data))
# liczba prób
num.probes = length(sampleNames(data))

# Normalizacja danych
eset.mas5 = simpleaffy::call.exprs(data, "mas5")
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
control = c( "C002_Control_F_87", "C005_Control_M_91", "C006_Control_F_71", "C007_Control_M_54",
             "C014_Control_M91",  "C015_Control_M_91", "C013_Control_M_61", "C011_Control_F_89",
             "C017_Control_M_72", "C016_Control_M_60", "C018_Control_F_67", "C010_Control_F_90",
             "C019_Control_F_46", "C020_Control_M_25", "C022_Control_F_68", "C024_Control_F_78",
             "C025_Control_F_88", "C023_Control_M_38", "C026_Control_M_58", "C021_Control_M_25",
             "C027_Control_M_90", "C008_Control_F_94", "C028_Control_F_54")
data = cel
control.array = match(control, sampleNames(data))
Difference <- rowMeans(exprs.eset[,-control.array]) - rowMeans(exprs.eset[,control.array])
Average <- rowMeans(exprs.eset)
plot(Average, Difference)
lines(lowess(Average, Difference),
      col = 'red', lwd = 4)
abline(h = -2)
abline(h = 2)

control.array = match(control, sampleNames(data))
Difference <- rowMeans(exprs.eset[,-control.array]) - rowMeans(exprs.eset[,control.array])
control.idx = which(sampleNames(data) %in% control)
fac = c()
fac[1:num.probes] = "Test"
fac[control.idx] = "Control"
population.groups <- factor(fac)
design <- model.matrix(~population.groups)
fit <- lmFit (eset.mas5, design)
fit.eBayes <- eBayes (fit)

lodd <- -log10(fit.eBayes$p.value[,2])
o1 <- order(abs(Difference), decreasing = TRUE)[1:30]
oo2 <- order(abs (fit.eBayes$t[,2]), decreasing = TRUE)[1:30]
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
text(Difference[ii], lodd[ii], row.names(data)[ii], pos = 3)
title("Volcano Plot with moderated t-statistics")
grid()
abline(v = c(-1,1), lwd = 0.5)
abline(h= c(-log10(0.05),-log10(0.01),-log10(0.001)), 
       col=c('red','blue','black'), lwd = 0.5)
#legend("topleft", , pch = c(5, 18), col = c('red', 'blue'), cex = 0.75, bg="transparent",
#       legend = c("First 50 with the lowest p value","First 50 with the highest Fold-change value"))

par(mar = c(10,10,10,10))
ii = fit.eBayes_ii$ii
exprs.eset = exprs(data.norm)
ii.mat <- exprs.eset[ii,]
ii.df <- data.frame(ii.mat)
brewer.cols <- brewer.pal(num.probes, "Set1")
hmcol <- colorRampPalette(brewer.pal(num.probes, 'Greys'))(256)
spcol = c()
spcol[1:num.probes] = 'grey10'
spcol[which(dimnames(ii.mat)[[2]] %in% control)] = 'grey80'
heatmap (ii.mat, col = hmcol, ColSideColors = spcol)#,margins = c (10,15))
heatmap (ii.mat, ColSideColors = spcol)

where_control = which(sampleNames(data.norm) %in% control)
col_groups = c()
col_groups[1:num.probes] = "Test"
col_groups[where_control] = "Control"

mat_col = data.frame(group = col_groups)

mat_colors = list(group = brewer.pal(2,"Set1"))
mat_colors$group = mat_colors$group[1:2]

rownames(mat_col) = colnames(data.norm)

pheatmap(mat = abs(ii.mat), annotation_col = mat_col, clustering_method = "complete")


tab = topTable(fit.eBayes, coef = 2, adjust.method = "BH", number = 50,#length(row.names(data)), 
               sort.by = "p")
tab

topgenes = tab[tab[, "adj.P.Val"] < 0.05]
dim(topgenes)
rrr=tab[,"P.Value"] < 0.05
topgenes = tab[rrr]
dim(topgenes)

library(hgu133plus2.db)
library(AnnotationDbi)
my_gen <- mapIds(hgu133plus2.db, "1558631_at", column = c("SYMBOL"),keytype="PROBEID")
                     
# kobirty 1558631_at

im = order(tab$P.Value)
imtab = tab[im,]
View(imtab)

# w publikacji podają ponas 1600 instotnych genów >0.05 i intensywnosc >30(mediana)
ppp = which(tab[,"AveExpr"] > 4.9)
pv = which(tab[,"P.Value"] < 0.05)
iii = intersect(ppp, pv)
length(iii)

mediany.norm = apply(exprs(data.norm),1,median)
ppp2 = which(mediany.norm > log2(30))
iii2 = intersect(ppp2, pv)

 annotatedTopTable <- function(topTab, anotPackage)
   {
     topTab <- cbind(PROBEID=rownames(topTab), topTab)
     myProbes <- rownames(topTab)
     thePackage <- eval(parse(text = anotPackage))
     geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID",
                                                  "GENENAME"))
     annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="PROBEID",
                              by.y="PROBEID", sort = TRUE)
     return(annotatedTopTab)
     }
g = annotatedTopTable(tab, paste(data@annotation, ".db", sep = ""))#"hgu133plus2.db")
g
order_pval = order(g[,"P.Value"])
g[order_pval,]

My_volcano_moderated_threshold <- function(data.norm, control, my_index, tab){
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
  
  #lodd <- -log10(fit.eBayes$p.value[,2])
  lodd <- -log10(tab[,"adj.P.Val"])
  
  plot (Difference[-my_index], lodd[-my_index],
        cex = .25, xlim = c (-3,3),
        ylim = range(lodd), xlab = 'Average (log) Fold-change',
        ylab = 'LOD score – Negative log10 of P-value')
  points(Difference [my_index], lodd [my_index],
         pch = 5, col ='red', cex = 1,
         lwd = 1)
  points (Difference[my_index], lodd[my_index],
          pch = 18, col = 'blue', cex = 1, lwd = 1)
  #text(Difference[my_index], lodd[my_index], row.names(data.norm)[my_index], pos = 3)
  title("Volcano Plot with moderated t-statistics")
  grid()

  abline(v = c(-2,2), lwd = 0.5, col = c('red'))
  abline(h = -log10(0.05), 
         col=c('red'), lwd = 0.5)
}

subb = substr(colnames(data), start = 1, stop = 4)
idx = which(subb == "norm")
control = colnames(data)[idx]
thresholded_p.val = which(tab[,"adj.P.Val"] < 0.05)
thresholded_f.c = which(abs(tab[,"logFC"]) > 2)
iii = intersect(thresholded_f.c, thresholded_p.val)
My_volcano_moderated_threshold(data.norm, control, iii, tab)
new_tab = tab[iii,]
g = annotatedTopTable(new_tab, "hgu95av2.db")
my_order = order(g[,"adj.P.Val"])
g[my_order,]

library("genefilter")
library("RColorBrewer")
library("affy")
library("simpleaffy")
library("affyPLM")


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

# PLM i NUSE
dataPLM = fitPLM(cel, normalize.method = "mas5")
affyPLM::boxplot(dataPLM, main="NUSE", #ylim = c(0.95, 1.22),
        # outline = FALSE, col="lightblue", 
        las=3, whisklty=0, staplelty=0)

Mbox(dataPLM, main="RLE", #ylim = c(-0.4, 0.4),
     # outline = FALSE, col="mistyrose", 
     las=3, whisklty=0, staplelty=0)

###
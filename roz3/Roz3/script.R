if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")


BiocManager::install("CLL")
BiocManager::install("genefilter")

library("genefilter")
library("affy")
library("CLL")
data("CLLbatch")
data("disease")

sampleNames(CLLbatch)
sampleNames(CLLbatch) = sub("\\.CEL$", "", sampleNames(CLLbatch))

head(disease)

rownames(disease) = disease$SampleID

mt = match(rownames(disease), sampleNames(CLLbatch))


vmd = data.frame(labelDescription = c("Sample ID", "Disease status: progressive or stable disease"))

phenoData(CLLbatch) = new("AnnotatedDataFrame", data = disease[mt, ], varMetadata = vmd)

CLLbatch = CLLbatch[, !is.na(CLLbatch$Disease)]

# BiocManager::install("affyQCReport")
library("affyQCReport")
library("simpleaffy")
# qc stats
eset.mas5 = call.exprs(CLLbatch, "mas5")
# tak szybciej, wynik ten sam
my_qc = simpleaffy::qc(CLLbatch, eset.mas5)
journalpng(file = "Stats_QC.png")
qc.stat.plot = simpleaffy::plot(my_qc, type=)
dev.off()

# CEL = affy::ReadAffy(celfile.path = "/home/kamil/Pulpit/magisterka/dane tam/cel files")
# CEL.mas5 = call.exprs(CEL, "mas5")
# CEL.qc = simpleaffy::qc(CEL, CEL.mas5)
# journalpng(file = "brain_Stats_QC_3_5.png", height = 10)
# simpleaffy::plot(CEL.qc)
# dev.off()
# journalpng(file = "brain_Stats_QC_3_M.png", height = 10)
# simpleaffy::plot(CEL.qc, usemid = T)
# dev.off()


dd = dist2(log2(exprs(CLLbatch)))
diag(dd)
dd.row <- as.dendrogram(hclust(as.dist(dd)))
row.ord <- order.dendrogram(dd.row)
library("latticeExtra")
legend = list(top=list(fun=latticeExtra::dendrogramGrob,
                         args=list(x=dd.row, side="top")))
lp = levelplot(dd[row.ord, row.ord],
                 scales=list(x=list(rot=90)), xlab="",
                 ylab="", legend=legend)
journalpng(file = "heatmap.png", height = 5, width = 5)
lp
dev.off()

library("affyPLM")
dataPLM = fitPLM(CLLbatch)
boxplot(dataPLM, main="NUSE", ylim = c(0.95, 1.22),
        outline = FALSE, col="lightblue", las=3,
        whisklty=0, staplelty=0)
Mbox(dataPLM, main="RLE", ylim = c(-0.4, 0.4),
     outline = FALSE, col="mistyrose", las=3,
     whisklty=0, staplelty=0)
# usunięcie tego, czego nie chcemy.
badArray = match("CLL1", sampleNames(CLLbatch))
CLLB = CLLbatch[, -badArray]
dataPLM = fitPLM(CLLB)
boxplot(dataPLM, main="NUSE", ylim = c(0.95, 1.22),
        outline = FALSE, col="lightblue", las=3,
        whisklty=0, staplelty=0)
Mbox(dataPLM, main="RLE", ylim = c(-0.4, 0.4),
     outline = FALSE, col="mistyrose", las=3,
     whisklty=0, staplelty=0)

###

CLLrma = rma(CLLB)
e = exprs(CLLrma)
dim(CLLrma)

where_progress = which(CLLrma$Disease %in% "progres.")
where_stable = which(CLLrma$Disease %in% "stable")
SampleID_progress = CLLrma$SampleID[where_progress]
SampleID_stable = CLLrma$SampleID[where_stable]
exprs.progres = exprs(CLLrma[,SampleID_progress])
exprs.stable = exprs(CLLrma[,SampleID_stable])

# tablica intensywności
average.intensity = data.frame(stable = rowMeans(exprs.stable), progres = rowMeans(exprs.progres))
# tablica logarytmów intensywności
logs.intensity = e
average.logs = data.frame(stable = rowMeans(exprs.stable), progres = rowMeans(exprs.progres))
# testowanie istotności różnicy logarytmów stable vs progres
log.ttests = genefilter::rowttests(as.matrix(logs.intensity),CLLrma$Disease)
# indeksy, gdzie p jest istotne oraz tablica tylko z nimi
idx = which(log.ttests$p.value < 0.05, c(TRUE))

all.logs = average.logs$progres - average.logs$stable
plot(rowMeans(average.intensity),all.logs, pch = ".", main = "Scatterplot of log-ratio vs mean intensity")
abline(h=0)

plot(rank(average.intensity$progres),all.logs, pch = ".",main = "Scatterplot of log-ratio vs rank of mean intensity")
abline(h=0)

lod1 = -log10(log.ttests$p.value)
plot(all.logs, lod1, pch=".", xlab="log-ratio",ylab=expression(-log[10]~p),main='Volcano plot (t-tests)')
abline(h=-log10(0.05),col='red')
abline(h=-log10(0.01),col='blue')

library("limma")
design = model.matrix(~CLLrma$Disease)
CLLlim = lmFitC(CLLrma, design)
CLLeb = limma::ebayes(CLLlim)
## kopia dla bayes
#idx = which(intensity.ttests$p.value < 0.05, c(TRUE))
lod2 = -log10(CLLeb$p.value[,2])
plot(all.logs, lod2, pch=".", xlab="log-ratio",ylab=expression(-log[10]~p),main='Volcano plot (Bayesian)')
abline(h=-log10(0.05),col='red')
abline(h=-log10(0.01),col='blue')












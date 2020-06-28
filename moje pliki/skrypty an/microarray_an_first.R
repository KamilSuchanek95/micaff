library(affy)

#### first ####

# wczytaj pliki .CEL
affy.data = affy::ReadAffy(celfile.path = 'moje pliki/dane/Su_CELs')
# normalizacja dla każdej sondy
eset.mas5 = affy::mas5(affy.data)
#  tworzenie macierzy ekpresji
exprSet.nologs = affy::exprs(eset.mas5)
# lista chipów
colnames(exprSet.nologs)
# zmiana nazw kolumn (bez .CEL)
colnames(exprSet.nologs) = c("brain.1", "brain.2", 
                             "fetal.brain.1", "fetal.brain.2",
                             "fetal.liver.1", "fetal.liver.2", 
                             "liver.1", "liver.2")
# zmiana rozkładu na normalny oraz ulogicznienie danych 
# by były ujemną/dotatnią ekpresją za pomocą log2
exprSet = log(exprSet.nologs, 2)
# zapisz do pliku
write.table(exprSet, file="Su_mas5_matrix.csv", quote=F, sep="\t")
# ewentualne obecności/nieobecności dla zestwów sond
data.mas5calls = affy::mas5calls(affy.data)
data.mas5calls.calls = affy::exprs(data.mas5calls)
# zapis do pliku
write.table(data.mas5calls.calls, file="Su_mas5calls.txt", quote=F, sep="\t")

# średnie dla 2 chipów
brain.mean = apply(exprSet[, c("brain.1", "brain.2")], 1, mean)
fetal.brain.mean = apply(exprSet[, c("fetal.brain.1", "fetal.brain.2")], 1, mean)
liver.mean = apply(exprSet[, c("liver.1", "liver.2")], 1, mean)
fetal.liver.mean = apply(exprSet[, c("fetal.liver.1", "fetal.liver.2")], 1, mean)
# stosunek fetal/adult to logarytmy więc odejmowanie odpowiada dzieleniu
brain.fetal.to.adult = fetal.brain.mean - brain.mean 
liver.fetal.to.adult = fetal.liver.mean - liver.mean
# zapis danych do jednej struktury
all.data = cbind (exprSet, brain.mean, fetal.brain.mean, liver.mean, fetal.liver.mean, 
                  brain.fetal.to.adult, liver.fetal.to.adult) 
# podgląd
colnames (all.data)
# zapis do następnej części
write.table (all.data, file = "Microarray_Analysis_data_1_SOLUTION.txt", quote = F, sep = "\t")

# example writing: write.table(expression.plus.pvals, "Su_mas5_DE_analysis.txt", sep="\t", quote=F)

#### second ####
# odczytanie pliku z ekpresjami log2
exprSet = read.delim("Su_mas5_matrix.txt")

# dataset.1 = exprSet[1, c("brain.1", "brain.2")]
# dataset.2 = exprSet[1, c("fetal.brain.1", "fetal.brain.2")]
# t.test.gene.1 = t.test(dataset.1, dataset.2, "two.sided")
# t.test.gene.1$p.value
# testy t Welcha dla troche małych grup (tylko po 2 zestawy...)
brain.p.value.all.genes = apply(exprSet, 1, function(x) { t.test(x[1:2], x[3:4]) $p.value } )
liver.p.value.all.genes = apply(exprSet, 1, function(x) { t.test(x[5:6], x[7:8]) $p.value } )
# odczytanie obecnych/nieobecnych z pliku
data.mas5calls.calls = read.delim("Su_mas5calls.txt")

# dla jednego genu
# AP.gene.1 = paste(data.mas5calls.calls[1,], collapse="")
# stworzenie ramki z konkatenacją AAAA APAAA itd.
AP = apply(data.mas5calls.calls, 1, paste, collapse="")

genes.present = names(AP[AP != "AAAAAAAA"])
# jak wiele zestawów sond jest oznaczonych jako obecne?
length(genes.present)
# weźmiemy pod uwagę gdzie gen jest obecny w choć 1 sondzie
exprSet.present = exprSet[genes.present,]
# surowe wartości p dla obecnych genów
brain.raw.pvals.present = brain.p.value.all.genes[genes.present]
liver.raw.pvals.present = liver.p.value.all.genes[genes.present]
# wartości p po korekcji...
brain.fdr.pvals.present = p.adjust(brain.raw.pvals.present, method="fdr")
liver.fdr.pvals.present = p.adjust(liver.raw.pvals.present, method="fdr")
# sortowanie po fdr
brain.fdr.pvals.present.sorted = 
  brain.fdr.pvals.present[order(brain.fdr.pvals.present)]
liver.fdr.pvals.present.sorted = 
  liver.fdr.pvals.present[order(liver.fdr.pvals.present)]
# złączenie ramki obecnych ekspresji log2 i surowych wartości p ich testowania
expression.plus.pvals = cbind(exprSet.present, brain.raw.pvals.present, 
                              brain.fdr.pvals.present, liver.raw.pvals.present, liver.fdr.pvals.present)
write.table(expression.plus.pvals, "Su_mas5_DE_analysis.txt", sep="\t", quote=F)
# wartości FDR są okropne, ale nie mamy więcej danych więc weźmiemy surowe p dla powiedzmy 0.01
# jednak musimy liczyc się z pewnymi wynikami FP
brain.DE.probesets = names(brain.raw.pvals.present[brain.raw.pvals.present < 0.01])
liver.DE.probesets = names(liver.raw.pvals.present[liver.raw.pvals.present < 0.01])

# odczytanie podsumowania poprzedniej części
all.data = read.delim("Microarray_Analysis_data_1_SOLUTION.txt")
# zabieramy obecne wiersze dla p>0.01, z dodatkową kolumną... 
brain.DE.log2.ratios = all.data[brain.DE.probesets, 
                                c("brain.fetal.to.adult", "liver.fetal.to.adult")]
liver.DE.log2.ratios = all.data[liver.DE.probesets, 
                                c("brain.fetal.to.adult", "liver.fetal.to.adult")]

write.table(brain.DE.log2.ratios, "brain.DE.log2.ratios.txt", sep="\t", quote=F)
write.table(liver.DE.log2.ratios, "liver.DE.log2.ratios.txt", sep="\t", quote=F)
# unia tych wektorów dla wszystkich wartości, jeśli tylko w jednym z 2 porównań uznaliśmy różnicę
# w ekspresji za istotną, czyli jeśli w brain jest istotna to ten sam gen również przepisano dla liver
union.DE.probesets = sort(union(brain.DE.probesets, liver.DE.probesets))
union.DE.log2.ratios = all.data[union.DE.probesets, 
                                c("brain.fetal.to.adult", "liver.fetal.to.adult")]
write.table(union.DE.log2.ratios, "union.DE.log2.ratios.txt", sep="\t", quote=F)
# unia, ale uwzględniająca istotne tylko dla mózgu i wątroby (4 zamiast 2 kolumn wyżej)
union.DE.log2.ratios.uni = transform(merge(brain.DE.log2.ratios[c("brain.fetal.to.adult")],
                                           liver.DE.log2.ratios[c("liver.fetal.to.adult")],
                                           by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
write.table(union.DE.log2.ratios.uni, "union.DE.log2.ratios.uni.txt", sep="\t", quote=F)

# klastrowanie


a = affy::ReadAffy('C:/Users/kamis/Desktop/micaff/moje pliki/dane/Su_CELs/Brain_1.CEL','C:/Users/kamis/Desktop/micaff/moje pliki/dane/Su_CELs/Brain_2.CEL')
a2 = affy::mas5(a)
#  tworzenie macierzy ekpresji
a3 = affy::exprs(a2)
a4 = log(a3, 2)
plot1 = affy::cdfName(a)

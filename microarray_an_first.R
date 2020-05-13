library(affy)
# wczytaj pliki .CEL
affy.data = affy::ReadAffy(celfile.path = "Su_CELs")
# normalizacja dla każdej sondy
eset.mas5 = affy::mas5(affy.data)
#  tworzenie macierzy ekpresji
exprSet.nologs = exprs(eset.mas5)
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
data.mas5calls = mas5calls(affy.data)
data.mas5calls.calls = exprs(data.mas5calls)
# zapis do pliku
write.table(data.mas5calls.calls, file="Su_mas5calls.csv", quote=F, sep="\t")

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
write.table (all.data, file = "Microarray_Analysis_data_1_SOLUTION.csv", quote = F, sep = "\ t")




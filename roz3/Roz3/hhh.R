dataPLM = fitPLM(cel, background.method = "MAS")

Mbox(dataPLM, main="RLE Plot", col = brewer.cols,
     outline = FALSE,
     las=3, whisklty=0, staplelty=0)
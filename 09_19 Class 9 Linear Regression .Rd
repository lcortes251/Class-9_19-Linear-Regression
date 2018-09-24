Github Class Exercise 1; Exercise 2R
Lines 73-76, 80-86, 88-89, 112, 114-122 
library(affy)
> setwd("/Users/Laura/Desktop/Estrogen")
> targetsFile <- "estrogen.txt"
> pd <- read.AnnotatedDataFrame(targetsFile,header=TRUE,sep="",row.names=1)
> ER <- pData(pd)$estrogen
> Time <- factor(pData(pd)$time.h)
> design <- model.matrix(~ER+Time)
> design

design2 <- model.matrix(~ER*Time)
> design2

#matrix was made in order to make data into binary numbers. Design matrix function has a formula ER + Time. Design 2 object has ER * time.
For additive matrix, there is no interaction. Look up later. I think the p values are added? And if one is significant  and the other isn't, it'll still be significant. 

For multiplication, if there is an interaction matrix then both have to be significant? for it to be significant??  

raw <-ReadAffy(celfile.path = "/Users/Laura/Desktop/Estrogen", filenames=rownames(pData(pd)),phenoData = pd)
> raw
eset <- rma(raw)
library(limma)
fit1 <- lmFit(eset, design)
> fit1 <- eBayes(fit1)
> topTable(fit1, coef=2)

fit2 <- lmFit(eset, design2)
> 
> fit2 <- eBayes(fit2)
> topTable(fit2, coef=2)

#we get raw data and put it in RMA normalization. 
#eset is expression data. 
#lmfit is linear fit. 
#ER is y axis and time is x axis
#eBayes is a bayesian function/algoritm to calculate p values and fold change. 

#annotation of probe ids 

#download gene filter


source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")

library(genefilter)

library(GEOquery)
library(limma)
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33126/matrix/GSE33126_series_matrix.txt.gz"

filenm <- "GSE33126_series_matrix.txt.gz"
if(!file.exists(filenm)) download.file(url, destfile=filenm)

gse <- getGEO(filename=filenm)

gse.expfilter <- varFilter(gse)

#make volcano plot on fit2 

volcanoplot(fit2, coef = 1L, style = "p-value", highlight = 0, names = fit$genes$ID, xlab = "Log2 Fold Change", ylab=NULL, pch=16, cex=0.35)




















## ----setup, include=FALSE------------------------------------------------
# smaller font size for chunks
opts_chunk$set(size = 'footnotesize')

## ------------------------------------------------------------------------
library("snpStats")
load("datos/snpsSNPstats.RData")
snpsSNPstats

## ------------------------------------------------------------------------
sum <- col.summary(snpsSNPstats)
dim(sum)
head(sum)

## ------------------------------------------------------------------------
Callrate <- sum$Call.rate
selectCallRate <- Callrate > 0.8
length(selectCallRate)
head(selectCallRate)

## ------------------------------------------------------------------------
MAF <- sum$MAF
selectMAF <- MAF > 0.01
length(selectMAF)
head(selectMAF)

## ------------------------------------------------------------------------
selectMAFCAllrete <- selectMAF & selectCallRate
head(selectMAFCAllrete)
table(selectMAFCAllrete)

## ------------------------------------------------------------------------
snpnames <- colnames(snpsSNPstats)
length(snpnames)
head(snpnames)
selsnpnames <- snpnames[selectMAFCAllrete]
length(selsnpnames)
head(selsnpnames)

## ------------------------------------------------------------------------
NewsnpsSNPstats<-snpsSNPstats[,selsnpnames]
NewsnpsSNPstats

## ------------------------------------------------------------------------
selhw<-abs(sum$z.HWE) < 6
NewsnpsSNPstats<-snpsSNPstats[, selhw & selectMAF & selectCallRate]
save(NewsnpsSNPstats, file="NewsnpsSNPstats.RData")


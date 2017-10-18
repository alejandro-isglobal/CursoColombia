## ----setup, include=FALSE------------------------------------------------
# smaller font size for chunks
opts_chunk$set(size = 'footnotesize')

## ----message=FALSE-------------------------------------------------------
library(snpStats)

## ----plot1,fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
load("datos/NewsnpsSNPstats.RData")
xxmat <- xxt(NewsnpsSNPstats)
evv <- eigen(xxmat)
pcs <- evv$vectors[,1:2]
plot(pcs)

## ----plot2, fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
ids<-read.table("datos/20130606_g1k.ped", sep="\t", header=TRUE)

rownames(ids)<-ids$Individual.ID
pops<-ids[rownames(NewsnpsSNPstats),]$Population

plot(pcs, col=as.numeric(pops), pch=16)
legend("topright", legend=levels(pops), pch=16, col=1:length(levels(pops)), cex=0.7)


## ----setup, include=FALSE------------------------------------------------
# smaller font size for chunks
opts_chunk$set(size = 'footnotesize')

## ----message=FALSE-------------------------------------------------------
library(snpStats)

## ------------------------------------------------------------------------
load("datos/NewsnpsSNPstats.RData")
ls()

## ------------------------------------------------------------------------
snps <- NewsnpsSNPstats[,c(90,91)]
ld(snps, stats=c("D.prime", "R.squared"),depth=1)

## ----plotLD, fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
LD <- ld(NewsnpsSNPstats, stats=c("D.prime", "R.squared"),depth=379)
image(LD$R.squared, lwd=0)

## ------------------------------------------------------------------------
ids <- read.table("datos/20130606_g1k.ped", sep="\t", header=TRUE)

rownames(ids) <- ids$Individual.ID
pops <- ids[rownames(NewsnpsSNPstats),]$Population
head(pops)
GBR <- pops=="GBR"
head(GBR)

## ----plotLDGRB, fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----

snpsGBR <- NewsnpsSNPstats[GBR,]
LDGBR <- ld(snpsGBR, stats=c("D.prime", "R.squared"),depth=379)
image(LDGBR$R.squared, lwd=0)

## ----plotLDYRI, fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold', echo=FALSE----

YRI <- pops=="YRI"
snpsYRI <- NewsnpsSNPstats[YRI,]
LDYRI <- ld(snpsYRI, stats=c("D.prime", "R.squared"),depth=379)
image(LDGBR$R.squared, lwd=0)
image(LDYRI$R.squared, lwd=0)

## ----eval=FALSE----------------------------------------------------------
## YRI <- pops=="YRI"
## snpsYRI <- NewsnpsSNPstats[YRI,]
## LDYRI <- ld(snpsYRI, stats=c("D.prime", "R.squared"),depth=379)
## image(LDYRI$R.squared, lwd=0)

## ------------------------------------------------------------------------
sumSnps <- col.summary(NewsnpsSNPstats)
head(sumSnps)

## ----fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
p<-sumSnps$MAF
hist(1-p^2-(1-p)^2)

## ------------------------------------------------------------------------
lv <- levels(pops)
x <- lv[1]
x
whichpop <- pops==x
POPsnps <- NewsnpsSNPstats[whichpop,]
sumSnps <- col.summary(POPsnps)
p <- sumSnps$MAF

## ----fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
hist(1-p^2-(1-p)^2)

## ----echo=FALSE, fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
lv <- levels(pops)
x <- lv[26]
x
whichpop <- pops==x
POPsnps <- NewsnpsSNPstats[whichpop,]
sumSnps <- col.summary(POPsnps)
p <- sumSnps$MAF
hist(1-p^2-(1-p)^2)

## ----eval=FALSE, fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
## lv <- levels(pops)
## x <- lv[26]
## x
## whichpop <- pops==x
## POPsnps <- NewsnpsSNPstats[whichpop,]
## sumSnps <- col.summary(POPsnps)
## p <- sumSnps$MAF
## hist(1-p^2-(1-p)^2)

## ------------------------------------------------------------------------
hetpop<-sapply(levels(pops), function(x)
 {
    whichpop <- pops==x
    POPsnps <- NewsnpsSNPstats[whichpop,]
    sumSnps <- col.summary(POPsnps)
    p <- sumSnps$MAF
    1-p^2-(1-p)^2
  })
hetpop[1:5,1:5]


## ----fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
boxplot(hetpop)

## ----plotfst, fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
FST <- Fst(NewsnpsSNPstats,pops)
hist(FST$Fst)

## ----echo=FALSE,fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
CEUandYRI <- pops%in%c("CEU","YRI")
snpsCEUandYRI<-NewsnpsSNPstats[CEUandYRI,]
popsCEUandYRI <- pops[CEUandYRI]
FST <- Fst(snpsCEUandYRI,popsCEUandYRI)
hist(FST$Fst)
mean(FST$Fst, na.rm=TRUE)

## ----eval=FALSE----------------------------------------------------------
## CEUandYRI <- pops%in%c("CEU","YRI")
## snpsCEUandYRI<-NewsnpsSNPstats[CEUandYRI,]
## popsCEUandYRI <- pops[CEUandYRI]
## FST <- Fst(snpsCEUandYRI,popsCEUandYRI)
## hist(FST$Fst)
## mean(FST$Fst, na.rm=TRUE)


## ----setup, include=FALSE------------------------------------------------
# smaller font size for chunks
opts_chunk$set(size = 'footnotesize')

## ----library-------------------------------------------------------------
library("snpStats")

load("datos/snp.RData")
snp

## ------------------------------------------------------------------------
phenos<-read.table("datos/phenosCont.txt",header=TRUE)
head(phenos)

## ------------------------------------------------------------------------
results <- snp.rhs.tests(phenos$caco~1, snp.data=snp)
results[1:10]

## ----plot1,fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
qq.chisq(chi.squared(results), 1)

## ----plot2,fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
resultsAd <- snp.rhs.tests(phenos$caco~phenos$pop, snp.data=snp)
qq.chisq(chi.squared(resultsAd), 1)

## ----plot3,fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
snpnum<-matrix(as.numeric(snp),ncol=ncol(snp))
d<-dist(snpnum, method="manhattan")
mds<- cmdscale(d,eig=TRUE, k=2)$points
plot(mds,col=phenos$pop)

## ----plot4,fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
resultsAdPCA <- snp.rhs.tests(phenos$caco~mds[,1]+mds[,2], snp.data=snp)

qq.chisq(chi.squared(resultsAdPCA), 1)

## ------------------------------------------------------------------------
mod<-glm(phenos$caco~snpnum[,1]+mds[,1]+mds[,2], family="binomial")
summary(mod)

## ------------------------------------------------------------------------
summod<-summary(mod)
summod$coeff
summod$coeff[2,c(1,4)]

## ------------------------------------------------------------------------
pvals<-sapply(1:ncol(snpnum), function(j)
{  
  mod<-glm(phenos$caco~snpnum[,j]+mds[,1]+mds[,2], family="binomial")
  summod<-summary(mod)
  summod$coeff
  summod$coeff[2,c(1,4)]
})
head(t(pvals))

## ----plot5,fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
plot(p.value(resultsAdPCA), pvals[2,])

## ------------------------------------------------------------------------
pvalsPop<-sapply(1:ncol(snpnum), function(j)
{  
  mod<-glm(phenos$caco~snpnum[,j]+phenos$pop, family="binomial")
  summod<-summary(mod)
  summod$coeff
  summod$coeff[2,c(1,4)]
})


## ----plot6,fig.width=5, fig.height=5, out.width='.45\\linewidth', fig.show='hold'----
plot(pvals[2,],pvalsPop[2,])


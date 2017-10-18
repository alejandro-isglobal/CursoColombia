## ----setup, include=FALSE------------------------------------------------
# smaller font size for chunks
opts_chunk$set(size = 'footnotesize')

## ----eval=FALSE----------------------------------------------------------
## source("https://bioconductor.org/biocLite.R")
## biocLite("snpStats")

## ----library-------------------------------------------------------------
library("snpStats")

## ------------------------------------------------------------------------
snp<-read.plink("datos/mydata")

## ------------------------------------------------------------------------
bim<-read.table("datos/mydata.bim",as.is=TRUE)

rs<-bim[,2]
tb<-table(rs)
dup<-tb[tb>1]
selrs<-rs[!rs%in%names(dup)]

head(rs)

## ------------------------------------------------------------------------
snp<-read.plink("datos/mydata", select.snps = selrs)
names(snp)

## ------------------------------------------------------------------------
snp$genotypes

## ------------------------------------------------------------------------
head(snp$fam)

## ------------------------------------------------------------------------
head(snp$map)

## ----eval=FALSE----------------------------------------------------------
## save(snp, file="datos/mydata.RData")

## ----loadData------------------------------------------------------------
load("datos/snp.RData")

## ----snpStats------------------------------------------------------------
snp

## ----eval=FALSE----------------------------------------------------------
## source("http://bioconductor.org/biocLite.R")
## biocLite("VariantAnnotation")

## ----message=FALSE-------------------------------------------------------
library(VariantAnnotation)
fl<-"datos/17.43921017-43972966.ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"
vcf <- readVcf(fl, "hg19")
vcf

## ------------------------------------------------------------------------
genos<-geno(vcf)
names(genos)
dim(genos$GT)
genos$GT[1:5,1:5]

## ------------------------------------------------------------------------
snps<-genos$GT
snps[snps=="0|0"]<-0
snps[snps=="1|1"]<-2
snps[snps!=0 & snps !=2]<-1
snps[1:5,1:5]
save(snps, file="snpsMAPT.RData")

## ------------------------------------------------------------------------
library(snpStats)
snpsnew<-t(snps)
snpsnew[snps=="0"] <- 1
snpsnew[snps=="1"] <- 2
snpsnew[snps=="2"] <- 3

snpsSNPstats <- new("SnpMatrix", snpsnew)
print(as(snpsSNPstats[1:5,1:5], 'character'))
save(snpsSNPstats, file="snpsSNPstats.RData")


############################################################################
# Barn Swallow Z-linked and autosomal population differentiation
############################################################################

### Goal: compare distributions of population differentiation among Barn
### swallow subspecies and hybrid zones.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('../fst_results/')

library(data.table)
library(dplyr)
library(hexbin)
library(RColorBrewer)

### Read in Z-linked and autosomal scaffold lists---------------------------

list.z <- read.table('../processing_files/hirundo_rustica_scaffold_list.Z.txt',header=F)
list.a <- read.table('../processing_files/hirundo_rustica_scaffold_list.auto.txt',header=F)

### ------------------------------------------------------------------------
### Relative population differentiation (Fst)
### ------------------------------------------------------------------------

### Read in data------------------------------------------------------------

er.gu <- read.table('erythrogaster_gutturalis.100kb.windowed.weir.fst',header=T)
er.ru <- read.table('erythrogaster_rustica.100kb.windowed.weir.fst',header=T)
er.sa <- read.table('erythrogaster_savignii.100kb.windowed.weir.fst',header=T)
er.tr <- read.table('erythrogaster_transitiva.100kb.windowed.weir.fst',header=T)
er.ty <- read.table('erythrogaster_tytleri.100kb.windowed.weir.fst',header=T)
gu.sa <- read.table('gutturalis_savignii.100kb.windowed.weir.fst',header=T)
gu.tr <- read.table('gutturalis_transitiva.100kb.windowed.weir.fst',header=T)
gu.ty <- read.table('gutturalis_tytleri.100kb.windowed.weir.fst',header=T)
ru.gu <- read.table('rustica_gutturalis.100kb.windowed.weir.fst',header=T)
ru.sa <- read.table('rustica_savignii.100kb.windowed.weir.fst',header=T)
ru.tr <- read.table('rustica_transitiva.100kb.windowed.weir.fst',header=T)
ru.ty <- read.table('rustica_tytleri.100kb.windowed.weir.fst',header=T)
sa.tr <- read.table('savignii_transitiva.100kb.windowed.weir.fst',header=T)
sa.ty <- read.table('savignii_tytleri.100kb.windowed.weir.fst',header=T)
tr.ty <- read.table('transitiva_tytleri.100kb.windowed.weir.fst',header=T)

gu.guty <- read.table('gutturalis_tytleri-gutturalis.100kb.windowed.weir.fst',header=T)
ty.guty <- read.table('tytleri_tytleri-gutturalis.100kb.windowed.weir.fst',header=T)
ru.rugu <- read.table('rustica_rustica-gutturalis.100kb.windowed.weir.fst',header=T)
gu.rugu <- read.table('gutturalis_rustica-gutturalis.100kb.windowed.weir.fst',header=T)
ru.ruty <- read.table('rustica_rustica-tytleri.100kb.windowed.weir.fst',header=T)
ty.ruty <- read.table('tytleri_rustica-tytleri.100kb.windowed.weir.fst',header=T)

### Subset data by Z and autosome lists--------------------------------------

er.gu.z <- setDT(er.gu)[CHROM %chin% list.z$V1]
er.gu.a <- setDT(er.gu)[CHROM %chin% list.a$V1]
er.ru.z <- setDT(er.ru)[CHROM %chin% list.z$V1]
er.ru.a <- setDT(er.ru)[CHROM %chin% list.a$V1]
er.sa.z <- setDT(er.sa)[CHROM %chin% list.z$V1]
er.sa.a <- setDT(er.sa)[CHROM %chin% list.a$V1]
er.tr.z <- setDT(er.tr)[CHROM %chin% list.z$V1]
er.tr.a <- setDT(er.tr)[CHROM %chin% list.a$V1]
er.ty.z <- setDT(er.ty)[CHROM %chin% list.z$V1]
er.ty.a <- setDT(er.ty)[CHROM %chin% list.a$V1]
gu.sa.z <- setDT(gu.sa)[CHROM %chin% list.z$V1]
gu.sa.a <- setDT(gu.sa)[CHROM %chin% list.a$V1]
gu.tr.z <- setDT(gu.tr)[CHROM %chin% list.z$V1]
gu.tr.a <- setDT(gu.tr)[CHROM %chin% list.a$V1]
gu.ty.z <- setDT(gu.ty)[CHROM %chin% list.z$V1]
gu.ty.a <- setDT(gu.ty)[CHROM %chin% list.a$V1]
ru.gu.z <- setDT(ru.gu)[CHROM %chin% list.z$V1]
ru.gu.a <- setDT(ru.gu)[CHROM %chin% list.a$V1]
ru.sa.z <- setDT(ru.sa)[CHROM %chin% list.z$V1]
ru.sa.a <- setDT(ru.sa)[CHROM %chin% list.a$V1]
ru.tr.z <- setDT(ru.tr)[CHROM %chin% list.z$V1]
ru.tr.a <- setDT(ru.tr)[CHROM %chin% list.a$V1]
ru.ty.z <- setDT(ru.ty)[CHROM %chin% list.z$V1]
ru.ty.a <- setDT(ru.ty)[CHROM %chin% list.a$V1]
sa.tr.z <- setDT(sa.tr)[CHROM %chin% list.z$V1]
sa.tr.a <- setDT(sa.tr)[CHROM %chin% list.a$V1]
sa.ty.z <- setDT(sa.ty)[CHROM %chin% list.z$V1]
sa.ty.a <- setDT(sa.ty)[CHROM %chin% list.a$V1]
tr.ty.z <- setDT(tr.ty)[CHROM %chin% list.z$V1]
tr.ty.a <- setDT(tr.ty)[CHROM %chin% list.a$V1]

gu.guty.z <- setDT(gu.guty)[CHROM %chin% list.z$V1]
gu.guty.a <- setDT(gu.guty)[CHROM %chin% list.a$V1]
ty.guty.z <- setDT(ty.guty)[CHROM %chin% list.z$V1]
ty.guty.a <- setDT(ty.guty)[CHROM %chin% list.a$V1]
ru.rugu.z <- setDT(ru.rugu)[CHROM %chin% list.z$V1]
ru.rugu.a <- setDT(ru.rugu)[CHROM %chin% list.a$V1]
gu.rugu.z <- setDT(gu.rugu)[CHROM %chin% list.z$V1]
gu.rugu.a <- setDT(gu.rugu)[CHROM %chin% list.a$V1]
ru.ruty.z <- setDT(ru.ruty)[CHROM %chin% list.z$V1]
ru.ruty.a <- setDT(ru.ruty)[CHROM %chin% list.a$V1]
ty.ruty.z <- setDT(ty.ruty)[CHROM %chin% list.z$V1]
ty.ruty.a <- setDT(ty.ruty)[CHROM %chin% list.a$V1]

### Summary stats-----------------------------------------------------------

# Calculate means for whole genome
mean(er.gu$MEAN_FST)
mean(er.ru$MEAN_FST)
mean(er.sa$MEAN_FST)
mean(er.tr$MEAN_FST)
mean(er.ty$MEAN_FST)
mean(ru.gu$MEAN_FST)
mean(gu.sa$MEAN_FST)
mean(gu.tr$MEAN_FST)
mean(gu.ty$MEAN_FST)
mean(ru.sa$MEAN_FST)
mean(ru.tr$MEAN_FST)
mean(ru.ty$MEAN_FST)
mean(sa.tr$MEAN_FST)
mean(sa.ty$MEAN_FST)
mean(tr.ty$MEAN_FST)
mean(gu.rugu$MEAN_FST)
mean(ru.rugu$MEAN_FST)
mean(gu.guty$MEAN_FST)
mean(ty.guty$MEAN_FST)
mean(ru.ruty$MEAN_FST)
mean(ty.ruty$MEAN_FST)

# Calculate means for autosomes and Z
mean(er.gu.a$MEAN_FST)
mean(er.ru.a$MEAN_FST)
mean(er.sa.a$MEAN_FST)
mean(er.tr.a$MEAN_FST)
mean(er.ty.a$MEAN_FST)
mean(ru.gu.a$MEAN_FST)
mean(gu.sa.a$MEAN_FST)
mean(gu.tr.a$MEAN_FST)
mean(gu.ty.a$MEAN_FST)
mean(ru.sa.a$MEAN_FST)
mean(ru.tr.a$MEAN_FST)
mean(ru.ty.a$MEAN_FST)
mean(sa.tr.a$MEAN_FST)
mean(sa.ty.a$MEAN_FST)
mean(tr.ty.a$MEAN_FST)
mean(gu.rugu.a$MEAN_FST)
mean(ru.rugu.a$MEAN_FST)
mean(gu.guty.a$MEAN_FST)
mean(ty.guty.a$MEAN_FST)
mean(ru.ruty.a$MEAN_FST)
mean(ty.ruty.a$MEAN_FST)

mean(er.gu.z$MEAN_FST)
mean(er.ru.z$MEAN_FST)
mean(er.sa.z$MEAN_FST)
mean(er.tr.z$MEAN_FST)
mean(er.ty.z$MEAN_FST)
mean(ru.gu.z$MEAN_FST)
mean(gu.sa.z$MEAN_FST)
mean(gu.tr.z$MEAN_FST)
mean(gu.ty.z$MEAN_FST)
mean(ru.sa.z$MEAN_FST)
mean(ru.tr.z$MEAN_FST)
mean(ru.ty.z$MEAN_FST)
mean(sa.tr.z$MEAN_FST)
mean(sa.ty.z$MEAN_FST)
mean(tr.ty.z$MEAN_FST)
mean(gu.rugu.z$MEAN_FST)
mean(ru.rugu.z$MEAN_FST)
mean(gu.guty.z$MEAN_FST)
mean(ty.guty.z$MEAN_FST)
mean(ru.ruty.z$MEAN_FST)
mean(ty.ruty.z$MEAN_FST)

# Two-sample Mann-Whitney U tests!
wilcox.test(er.gu.a$MEAN_FST,er.gu.z$MEAN_FST)
wilcox.test(er.ru.a$MEAN_FST,er.ru.z$MEAN_FST)
wilcox.test(er.sa.a$MEAN_FST,er.sa.z$MEAN_FST)
wilcox.test(er.tr.a$MEAN_FST,er.tr.z$MEAN_FST)
wilcox.test(er.ty.a$MEAN_FST,er.ty.z$MEAN_FST)
wilcox.test(ru.gu.a$MEAN_FST,ru.gu.z$MEAN_FST)
wilcox.test(gu.sa.a$MEAN_FST,gu.sa.z$MEAN_FST)
wilcox.test(gu.tr.a$MEAN_FST,gu.tr.z$MEAN_FST)
wilcox.test(gu.ty.a$MEAN_FST,gu.ty.z$MEAN_FST)
wilcox.test(ru.sa.a$MEAN_FST,ru.sa.z$MEAN_FST)
wilcox.test(ru.tr.a$MEAN_FST,ru.tr.z$MEAN_FST)
wilcox.test(ru.ty.a$MEAN_FST,ru.ty.z$MEAN_FST)
wilcox.test(sa.tr.a$MEAN_FST,sa.tr.z$MEAN_FST)
wilcox.test(sa.ty.a$MEAN_FST,sa.ty.z$MEAN_FST)
wilcox.test(tr.ty.a$MEAN_FST,tr.ty.z$MEAN_FST)

wilcox.test(gu.rugu.a$MEAN_FST,gu.rugu.z$MEAN_FST)
wilcox.test(ru.rugu.a$MEAN_FST,ru.rugu.z$MEAN_FST)
wilcox.test(gu.guty.a$MEAN_FST,gu.guty.z$MEAN_FST)
wilcox.test(ty.guty.a$MEAN_FST,ty.guty.z$MEAN_FST)
wilcox.test(ru.ruty.a$MEAN_FST,ru.ruty.z$MEAN_FST)
wilcox.test(ty.ruty.a$MEAN_FST,ty.ruty.z$MEAN_FST)

# Are FstZA ratios different for pairwise subspecies versus parental/hybrid comparisons?
fstZA.sub <- c(3.90,7.38,5.01,5.72,1.88,4.90,4.02,4.41,3.69,1.74,2.86,9.15,2.22,5.62,6.46)
fstZA.hyb <- c(7.98,3.45,4.03,2.79,10.51)

t.test(fstZA.sub,fstZA.hyb)

# No. Probably have too low of power to detect any difference.


### Plot scatterplot of ratios----------------------------------------------

mean.fst <- read.table('./fst/mean_Fst_summary.txt',header=T)
mean(mean.fst$mean.fst)
sd(mean.fst$mean.fst)

par(mfrow=c(1,2))
scatt <- palette(c('grey','black'))
palette()
color = rep(NA, length=length(mean.fst$category))
color[which(mean.fst$category=="subspecies")] = "black"
color[which(mean.fst$category=="parental-hybrid")] = "grey"

plot(mean.fst$mean.fst.auto,mean.fst$mean.fst.z,pch=20,xlim=c(0,0.2),ylim=c(0,0.2),col=color,xlab='Autosomal Fst',ylab='Z Chromosome Fst')
xax <- c(0,0.05,0.1,0.15,0.2)
yax <- c(0,0.05,0.1,0.15,0.2)
abline(lm(yax~xax),lty=1)

plot(mean.fst$mean.fst.z,mean.fst$ratio.fst.ZA,pch=20,col=color,xlab='Z Chromosome Fst',ylab='Fst ZA',ylim=c(0,10))
#points(mean.fst$mean.fst.auto,mean.fst$ratio.fst.ZA,pch=20,col='blue')
abline(lm(mean.fst$ratio.fst.ZA~mean.fst$mean.fst.z),lty=2)
cor.test(mean.fst$mean.fst.z,mean.fst$ratio.fst.ZA,method='spearman')

hyb <- mean.fst[which(mean.fst$category=='parental-hybrid'),]
cor.test(hyb$mean.fst.z,hyb$ratio.fst.ZA,method='spearman')

cor.test(mean.fst$mean.fst.auto,mean.fst$ratio.fst.ZA)

### Plot density distributions----------------------------------------------

par(mfrow=c(3,3))
plot(density(er.ty.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst erythrogaster-tytleri")
polygon(density(er.ty.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(er.ty.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(er.ty.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(gu.ty.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst gutturalis-tytleri")
polygon(density(gu.ty.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(gu.ty.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(gu.ty.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(ru.gu.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst gutturalis-rustica")
polygon(density(ru.gu.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(ru.gu.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(ru.gu.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(ru.ty.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst rustica-tytleri")
polygon(density(ru.ty.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(ru.ty.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(ru.ty.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(ru.sa.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst rustica-savignii")
polygon(density(ru.sa.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(ru.sa.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(ru.sa.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(ru.tr.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst rustica-transitiva")
polygon(density(ru.tr.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(ru.tr.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(ru.tr.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(sa.tr.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst savignii-transitiva")
polygon(density(sa.tr.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(sa.tr.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(sa.tr.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(gu.guty.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst gutturalis-gutturalis x tytleri")
polygon(density(gu.guty.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(gu.guty.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(gu.guty.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(ty.guty.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst tytleri-gutturalis x tytleri")
polygon(density(ty.guty.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(ty.guty.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(ty.guty.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(gu.rugu.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst gutturalis-gutturalis x rustica")
polygon(density(gu.rugu.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(gu.rugu.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(gu.rugu.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(ru.rugu.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst rustica-gutturalis x rustica")
polygon(density(ru.rugu.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(ru.rugu.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(ru.rugu.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(ru.ruty.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst rustica-rustica x tytleri")
polygon(density(ru.ruty.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(ru.ruty.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(ru.ruty.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(ty.ruty.a$MEAN_FST,na.rm=T),col='lightgrey',main=NA,xlab="Fst tytleri-rustica x tytleri")
polygon(density(ty.ruty.a$MEAN_FST,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(ty.ruty.z$MEAN_FST,na.rm=T),col='seagreen4')
polygon(density(ty.ruty.z$MEAN_FST,na.rm=T),col=alpha('seagreen4',0.5))


### ------------------------------------------------------------------------
### Population branch statistics (PBS) per subspecies
### ------------------------------------------------------------------------

# Prior to performing analyses, need to merge sets of Fst results:

pbs.er.data <-merge(er.ty,er.sa,by=c("CHROM","BIN_START")) 
pbs.er.data <-merge(pbs.er.data,sa.ty,by=c("CHROM","BIN_START"))
pbs.er <- ((pbs.er.data$WEIGHTED_FST.x + pbs.er.data$WEIGHTED_FST.y)-pbs.er.data$MEAN_FST)/2
pbs.er.data$pbs.er <- pbs.er

pbs.gu.data <-merge(gu.ty,gu.sa,by=c("CHROM","BIN_START")) 
pbs.gu.data <-merge(pbs.gu.data,sa.ty,by=c("CHROM","BIN_START"))
pbs.gu <- ((pbs.gu.data$WEIGHTED_FST.x + pbs.gu.data$WEIGHTED_FST.y)-pbs.gu.data$MEAN_FST)/2
pbs.gu.data$pbs.gu <- pbs.gu

pbs.ru.data <-merge(ru.gu,er.ru,by=c("CHROM","BIN_START")) 
pbs.ru.data <-merge(pbs.ru.data,er.gu,by=c("CHROM","BIN_START"))
pbs.ru <- ((pbs.ru.data$WEIGHTED_FST.x + pbs.ru.data$WEIGHTED_FST.y)-pbs.ru.data$MEAN_FST)/2
pbs.ru.data$pbs.ru <- pbs.ru

pbs.sa.data <-merge(ru.sa,er.sa,by=c("CHROM","BIN_START")) 
pbs.sa.data <-merge(pbs.sa.data,er.ru,by=c("CHROM","BIN_START"))
pbs.sa <- ((pbs.sa.data$WEIGHTED_FST.x + pbs.sa.data$WEIGHTED_FST.y)-pbs.sa.data$MEAN_FST)/2
pbs.sa.data$pbs.sa <- pbs.sa

pbs.tr.data <-merge(ru.tr,er.tr,by=c("CHROM","BIN_START")) 
pbs.tr.data <-merge(pbs.tr.data,er.ru,by=c("CHROM","BIN_START"))
pbs.tr <- ((pbs.tr.data$WEIGHTED_FST.x + pbs.tr.data$WEIGHTED_FST.y)-pbs.tr.data$MEAN_FST)/2
pbs.tr.data$pbs.tr <- pbs.tr

pbs.ty.data <-merge(ru.ty,gu.ty,by=c("CHROM","BIN_START")) 
pbs.ty.data <-merge(pbs.ty.data,ru.gu,by=c("CHROM","BIN_START"))
pbs.ty <- ((pbs.ty.data$WEIGHTED_FST.x + pbs.ty.data$WEIGHTED_FST.y)-pbs.ty.data$MEAN_FST)/2
pbs.ty.data$pbs.ty <- pbs.ty

# Write tables!
write.table(pbs.er.data,file='./pbs/pbs.er.data.txt',quote=F,row.names=F,sep = "\t")
write.table(pbs.gu.data,file='./pbs/pbs.gu.data.txt',quote=F,row.names=F,sep = "\t")
write.table(pbs.ru.data,file='./pbs/pbs.ru.data.txt',quote=F,row.names=F,sep = "\t")
write.table(pbs.sa.data,file='./pbs/pbs.sa.data.txt',quote=F,row.names=F,sep = "\t")
write.table(pbs.tr.data,file='./pbs/pbs.tr.data.txt',quote=F,row.names=F,sep = "\t")
write.table(pbs.ty.data,file='./pbs/pbs.ty.data.txt',quote=F,row.names=F,sep = "\t")

# Subset by autosomes and Z
pbs.er.data.z <- setDT(pbs.er.data)[CHROM %chin% list.z$V1]
pbs.er.data.a <- setDT(pbs.er.data)[CHROM %chin% list.a$V1]
pbs.gu.data.z <- setDT(pbs.gu.data)[CHROM %chin% list.z$V1]
pbs.gu.data.a <- setDT(pbs.gu.data)[CHROM %chin% list.a$V1]
pbs.ru.data.z <- setDT(pbs.ru.data)[CHROM %chin% list.z$V1]
pbs.ru.data.a <- setDT(pbs.ru.data)[CHROM %chin% list.a$V1]
pbs.sa.data.z <- setDT(pbs.sa.data)[CHROM %chin% list.z$V1]
pbs.sa.data.a <- setDT(pbs.sa.data)[CHROM %chin% list.a$V1]
pbs.tr.data.z <- setDT(pbs.tr.data)[CHROM %chin% list.z$V1]
pbs.tr.data.a <- setDT(pbs.tr.data)[CHROM %chin% list.a$V1]
pbs.ty.data.z <- setDT(pbs.ty.data)[CHROM %chin% list.z$V1]
pbs.ty.data.a <- setDT(pbs.ty.data)[CHROM %chin% list.a$V1]

# Boxplot!
par(mfrow=c(1,1))
boxplot(pbs.er.data.a$pbs.er,pbs.er.data.z$pbs.er,
        pbs.gu.data.a$pbs.gu,pbs.gu.data.z$pbs.gu,
        pbs.ru.data.a$pbs.ru,pbs.ru.data.z$pbs.ru,
        pbs.sa.data.a$pbs.sa,pbs.sa.data.z$pbs.sa,
        pbs.tr.data.a$pbs.tr,pbs.tr.data.z$pbs.tr,
        pbs.ty.data.a$pbs.ty,pbs.ty.data.z$pbs.ty,
        col=c('grey','seagreen4'))

# Z scans!
par(mfrow=c(6,1))
plot(pbs.er.data.z$pbs.er,pch=20,col=alpha('seagreen4',0.5),ylim=c(0,1))
plot(pbs.gu.data.z$pbs.gu,pch=20,col=alpha('seagreen4',0.5),ylim=c(0,1))
plot(pbs.ru.data.z$pbs.ru,pch=20,col=alpha('seagreen4',0.5),ylim=c(0,1))
plot(pbs.sa.data.z$pbs.sa,pch=20,col=alpha('seagreen4',0.5),ylim=c(0,1))
plot(pbs.tr.data.z$pbs.tr,pch=20,col=alpha('seagreen4',0.5),ylim=c(0,1))
plot(pbs.ty.data.z$pbs.ty,pch=20,col=alpha('seagreen4',0.5),ylim=c(0,1))

# Wrote Python script to order scans according to Flycatcher

# python order_scans_pbs.py ChrZ_order.txt pbs.er.data.txt pbs.er.chrZ.txt
# python order_scans_pbs.py ChrZ_order.txt pbs.gu.data.txt pbs.gu.chrZ.txt
# python order_scans_pbs.py ChrZ_order.txt pbs.ru.data.txt pbs.ru.chrZ.txt
# python order_scans_pbs.py ChrZ_order.txt pbs.sa.data.txt pbs.sa.chrZ.txt
# python order_scans_pbs.py ChrZ_order.txt pbs.tr.data.txt pbs.tr.chrZ.txt
# python order_scans_pbs.py ChrZ_order.txt pbs.ty.data.txt pbs.ty.chrZ.txt

# python order_scans_pbs.py Chr4_order.txt pbs.er.data.txt pbs.er.chr4.txt
# python order_scans_pbs.py Chr4_order.txt pbs.gu.data.txt pbs.gu.chr4.txt
# python order_scans_pbs.py Chr4_order.txt pbs.ru.data.txt pbs.ru.chr4.txt
# python order_scans_pbs.py Chr4_order.txt pbs.sa.data.txt pbs.sa.chr4.txt
# python order_scans_pbs.py Chr4_order.txt pbs.tr.data.txt pbs.tr.chr4.txt
# python order_scans_pbs.py Chr4_order.txt pbs.ty.data.txt pbs.ty.chr4.txt

# Also wrote shell script to kick out an ordered PBS scan per chromosome.

# Then concatenated the output like so:

# head -n 1 pbs.er.chr1.txt > pbs.er.data.chrom_order.txt; for chrom in `cat ../chrom.list`; do tail -n +2 pbs.er.$chrom.txt >> pbs.er.data.chrom_order.txt; done
# head -n 1 pbs.gu.chr1.txt > pbs.gu.data.chrom_order.txt; for chrom in `cat ../chrom.list`; do tail -n +2 pbs.gu.$chrom.txt >> pbs.gu.data.chrom_order.txt; done
# head -n 1 pbs.ru.chr1.txt > pbs.ru.data.chrom_order.txt; for chrom in `cat ../chrom.list`; do tail -n +2 pbs.ru.$chrom.txt >> pbs.ru.data.chrom_order.txt; done
# head -n 1 pbs.sa.chr1.txt > pbs.sa.data.chrom_order.txt; for chrom in `cat ../chrom.list`; do tail -n +2 pbs.sa.$chrom.txt >> pbs.sa.data.chrom_order.txt; done
# head -n 1 pbs.tr.chr1.txt > pbs.tr.data.chrom_order.txt; for chrom in `cat ../chrom.list`; do tail -n +2 pbs.tr.$chrom.txt >> pbs.tr.data.chrom_order.txt; done
# head -n 1 pbs.ty.chr1.txt > pbs.ty.data.chrom_order.txt; for chrom in `cat ../chrom.list`; do tail -n +2 pbs.ty.$chrom.txt >> pbs.ty.data.chrom_order.txt; done

# Read in ordered scans:

pbs.er.data.4 <- read.table('./pbs/pbs.er.chr4.txt',header=T)
pbs.gu.data.4 <- read.table('./pbs/pbs.gu.chr4.txt',header=T)
pbs.ru.data.4 <- read.table('./pbs/pbs.ru.chr4.txt',header=T)
pbs.sa.data.4 <- read.table('./pbs/pbs.sa.chr4.txt',header=T)
pbs.tr.data.4 <- read.table('./pbs/pbs.tr.chr4.txt',header=T)
pbs.ty.data.4 <- read.table('./pbs/pbs.ty.chr4.txt',header=T)

pbs.er.data.z <- read.table('./pbs/pbs.er.Z.txt',header=T)
pbs.gu.data.z <- read.table('./pbs/pbs.gu.Z.txt',header=T)
pbs.ru.data.z <- read.table('./pbs/pbs.ru.Z.txt',header=T)
pbs.sa.data.z <- read.table('./pbs/pbs.sa.Z.txt',header=T)
pbs.tr.data.z <- read.table('./pbs/pbs.tr.Z.txt',header=T)
pbs.ty.data.z <- read.table('./pbs/pbs.ty.Z.txt',header=T)

# Together, in Flycatcher order!

par(mfrow=c(6,2))
plot(pbs.er.data.4$pbs,type='l',lwd=1.5,col=alpha('grey',1),ylim=c(0,1),ylab='pbs')
plot(pbs.er.data.z$pbs,type='l',lwd=1.5,col=alpha('seagreen4',1),ylim=c(0,1),ylab='pbs')
plot(pbs.gu.data.4$pbs,type='l',lwd=1.5,col=alpha('grey',1),ylim=c(0,1),ylab='pbs')
plot(pbs.gu.data.z$pbs,type='l',lwd=1.5,col=alpha('seagreen4',1),ylim=c(0,1),ylab='pbs')
plot(pbs.ru.data.4$pbs,type='l',lwd=1.5,col=alpha('grey',1),ylim=c(0,1),ylab='pbs')
plot(pbs.ru.data.z$pbs,type='l',lwd=1.5,col=alpha('seagreen4',1),ylim=c(0,1),ylab='pbs')
plot(pbs.sa.data.4$pbs,type='l',lwd=1.5,col=alpha('grey',1),ylim=c(0,1),ylab='pbs')
plot(pbs.sa.data.z$pbs,type='l',lwd=1.5,col=alpha('seagreen4',1),ylim=c(0,1),ylab='pbs')
plot(pbs.tr.data.4$pbs,type='l',lwd=1.5,col=alpha('grey',1),ylim=c(0,1),ylab='pbs')
plot(pbs.tr.data.z$pbs,type='l',lwd=1.5,col=alpha('seagreen4',1),ylim=c(0,1),ylab='pbs')
plot(pbs.ty.data.4$pbs,type='l',lwd=1.5,col=alpha('grey',1),ylim=c(0,1),ylab='pbs')
plot(pbs.ty.data.z$pbs,type='l',lwd=1.5,col=alpha('seagreen4',1),ylim=c(0,1),ylab='pbs')

# Summary statistics!
mean(pbs.er.data$pbs.er)
mean(pbs.gu.data$pbs.gu)
mean(pbs.ru.data$pbs.ru)
mean(pbs.sa.data$pbs.sa)
mean(pbs.tr.data$pbs.tr)
mean(pbs.ty.data$pbs.ty)

mean(pbs.er.data.a$pbs.er)
mean(pbs.gu.data.a$pbs.gu)
mean(pbs.ru.data.a$pbs.ru)
mean(pbs.sa.data.a$pbs.sa)
mean(pbs.tr.data.a$pbs.tr)
mean(pbs.ty.data.a$pbs.ty)

mean(pbs.er.data.z$pbs.er)
mean(pbs.gu.data.z$pbs.gu)
mean(pbs.ru.data.z$pbs.ru)
mean(pbs.sa.data.z$pbs.sa)
mean(pbs.tr.data.z$pbs.tr)
mean(pbs.ty.data.z$pbs.ty)

# Run two-sample Mann-Whitney U tests to compare Z and auto distributions!
wilcox.test(pbs.er.data.a$pbs.er,pbs.er.data.z$pbs)
wilcox.test(pbs.gu.data.a$pbs.gu,pbs.gu.data.z$pbs)
wilcox.test(pbs.ru.data.a$pbs.ru,pbs.ru.data.z$pbs)
wilcox.test(pbs.sa.data.a$pbs.sa,pbs.sa.data.z$pbs)
wilcox.test(pbs.tr.data.a$pbs.tr,pbs.tr.data.z$pbs)
wilcox.test(pbs.ty.data.a$pbs.ty,pbs.ty.data.z$pbs)

# Now plot scans for each chromosome.

pbs.er.chrom <- read.table('./pbs/pbs.er.data.chrom_order.txt',header=T)
pbs.gu.chrom <- read.table('./pbs/pbs.gu.data.chrom_order.txt',header=T)
pbs.ru.chrom <- read.table('./pbs/pbs.ru.data.chrom_order.txt',header=T)
pbs.sa.chrom <- read.table('./pbs/pbs.sa.data.chrom_order.txt',header=T)
pbs.tr.chrom <- read.table('./pbs/pbs.tr.data.chrom_order.txt',header=T)
pbs.ty.chrom <- read.table('./pbs/pbs.ty.data.chrom_order.txt',header=T)

pbs.er.chrom$pbs[pbs.er.chrom$pbs<=0] <- 0
pbs.gu.chrom$pbs[pbs.gu.chrom$pbs<=0] <- 0
pbs.ru.chrom$pbs[pbs.ru.chrom$pbs<=0] <- 0
pbs.sa.chrom$pbs[pbs.sa.chrom$pbs<=0] <- 0
pbs.tr.chrom$pbs[pbs.tr.chrom$pbs<=0] <- 0
pbs.ty.chrom$pbs[pbs.ty.chrom$pbs<=0] <- 0

pbs.pal <- palette(c('steelblue','lightblue','orange','seagreen'))
length(unique(pbs.er.chrom$CHROM))
color = rep(NA, length=length(pbs.er.chrom$CHROM))
color[which(pbs.er.chrom$CHROM=="chr1")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr1A")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr2")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr3")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr4")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr4A")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr5")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr6")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr7")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr8")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr9")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr10")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr11")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr12")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr13")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr14")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr15")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr17")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr18")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr19")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr20")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr21")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr22")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr23")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr24")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr25")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr26")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="chr27")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="chr28")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="LG34")] = "seagreen"
color[which(pbs.er.chrom$CHROM=="LG35")] = "skyblue"
color[which(pbs.er.chrom$CHROM=="Z")] = "seagreen"

par(mfrow=c(6,1))
plot(pbs.er.chrom$pbs,pch=20,col=color)
plot(pbs.gu.chrom$pbs,pch=20,col=color)
plot(pbs.ru.chrom$pbs,pch=20,col=color)
plot(pbs.sa.chrom$pbs,pch=20,col=color)
plot(pbs.tr.chrom$pbs,pch=20,col=color)
plot(pbs.ty.chrom$pbs,pch=20,col=color)


### ------------------------------------------------------------------------
### Relative population differentiation correlations between subspecies
### ------------------------------------------------------------------------

# Merge PBS datasets for all pairwise subspecies

pbs.er.gu <- merge(pbs.er.data,pbs.gu.data,by=c("CHROM","BIN_START"))
pbs.er.ru <- merge(pbs.er.data,pbs.ru.data,by=c("CHROM","BIN_START"))
pbs.er.sa <- merge(pbs.er.data,pbs.sa.data,by=c("CHROM","BIN_START"))
pbs.er.tr <- merge(pbs.er.data,pbs.tr.data,by=c("CHROM","BIN_START"))
pbs.er.ty <- merge(pbs.er.data,pbs.ty.data,by=c("CHROM","BIN_START"))
pbs.gu.ru <- merge(pbs.gu.data,pbs.ru.data,by=c("CHROM","BIN_START"))
pbs.gu.sa <- merge(pbs.gu.data,pbs.sa.data,by=c("CHROM","BIN_START"))
pbs.gu.tr <- merge(pbs.gu.data,pbs.tr.data,by=c("CHROM","BIN_START"))
pbs.gu.ty <- merge(pbs.gu.data,pbs.ty.data,by=c("CHROM","BIN_START"))
pbs.ru.sa <- merge(pbs.ru.data,pbs.sa.data,by=c("CHROM","BIN_START"))
pbs.ru.tr <- merge(pbs.ru.data,pbs.tr.data,by=c("CHROM","BIN_START"))
pbs.ru.ty <- merge(pbs.ru.data,pbs.ty.data,by=c("CHROM","BIN_START"))
pbs.sa.tr <- merge(pbs.sa.data,pbs.tr.data,by=c("CHROM","BIN_START"))
pbs.sa.ty <- merge(pbs.sa.data,pbs.ty.data,by=c("CHROM","BIN_START"))
pbs.tr.ty <- merge(pbs.tr.data,pbs.ty.data,by=c("CHROM","BIN_START"))

# Scatterplots

par(mfrow=c(5,3))
plot(pbs.er.gu$pbs.er,pbs.er.gu$pbs.gu,pch=20,col=alpha('black',0.5),xlab='pbs erythrogaster',ylab='pbs gutturalis')
plot(pbs.er.ru$pbs.er,pbs.er.ru$pbs.ru,pch=20,col=alpha('black',0.5),xlab='pbs erythrogaster',ylab='pbs rustica')
plot(pbs.er.sa$pbs.er,pbs.er.sa$pbs.sa,pch=20,col=alpha('black',0.5),xlab='pbs erythrogaster',ylab='pbs savignii')
plot(pbs.er.tr$pbs.er,pbs.er.tr$pbs.tr,pch=20,col=alpha('black',0.5),xlab='pbs erythrogaster',ylab='pbs transitiva')
plot(pbs.er.ty$pbs.er,pbs.er.ty$pbs.ty,pch=20,col=alpha('black',0.5),xlab='pbs erythrogaster',ylab='pbs tytleri')
plot(pbs.gu.ru$pbs.gu,pbs.gu.ru$pbs.ru,pch=20,col=alpha('black',0.5),xlab='pbs gutturalis',ylab='pbs rustica')
plot(pbs.gu.sa$pbs.gu,pbs.gu.sa$pbs.sa,pch=20,col=alpha('black',0.5),xlab='pbs gutturalis',ylab='pbs savignii')
plot(pbs.gu.tr$pbs.gu,pbs.gu.tr$pbs.tr,pch=20,col=alpha('black',0.5),xlab='pbs gutturalis',ylab='pbs transitiva')
plot(pbs.gu.ty$pbs.gu,pbs.gu.ty$pbs.ty,pch=20,col=alpha('black',0.5),xlab='pbs gutturalis',ylab='pbs tytleri')
plot(pbs.ru.sa$pbs.ru,pbs.ru.sa$pbs.sa,pch=20,col=alpha('black',0.5),xlab='pbs rustica',ylab='pbs savignii')
plot(pbs.ru.tr$pbs.ru,pbs.ru.tr$pbs.tr,pch=20,col=alpha('black',0.5),xlab='pbs rustica',ylab='pbs transitiva')
plot(pbs.ru.ty$pbs.ru,pbs.ru.ty$pbs.ty,pch=20,col=alpha('black',0.5),xlab='pbs rustica',ylab='pbs tytleri')
plot(pbs.sa.tr$pbs.sa,pbs.sa.tr$pbs.tr,pch=20,col=alpha('black',0.5),xlab='pbs savignii',ylab='pbs transitiva')
plot(pbs.sa.ty$pbs.sa,pbs.sa.ty$pbs.ty,pch=20,col=alpha('black',0.5),xlab='pbs savignii',ylab='pbs tytleri')
plot(pbs.tr.ty$pbs.tr,pbs.tr.ty$pbs.ty,pch=20,col=alpha('black',0.5),xlab='pbs transitiva',ylab='pbs tytleri')

# Spearman's rank order correlations

cor.test(pbs.er.gu$pbs.er,pbs.er.gu$pbs.gu,method="spearman")
cor.test(pbs.er.ru$pbs.er,pbs.er.ru$pbs.ru,method="spearman")
cor.test(pbs.er.sa$pbs.er,pbs.er.sa$pbs.sa,method="spearman")
cor.test(pbs.er.tr$pbs.er,pbs.er.tr$pbs.tr,method="spearman")
cor.test(pbs.er.ty$pbs.er,pbs.er.ty$pbs.ty,method="spearman")
cor.test(pbs.gu.ru$pbs.gu,pbs.gu.ru$pbs.ru,method="spearman")
cor.test(pbs.gu.sa$pbs.gu,pbs.gu.sa$pbs.sa,method="spearman")
cor.test(pbs.gu.tr$pbs.gu,pbs.gu.tr$pbs.tr,method="spearman")
cor.test(pbs.gu.ty$pbs.gu,pbs.gu.ty$pbs.ty,method="spearman")
cor.test(pbs.ru.sa$pbs.ru,pbs.ru.sa$pbs.sa,method="spearman")
cor.test(pbs.ru.tr$pbs.ru,pbs.ru.tr$pbs.tr,method="spearman")
cor.test(pbs.ru.ty$pbs.ru,pbs.ru.ty$pbs.ty,method="spearman")
cor.test(pbs.sa.tr$pbs.sa,pbs.sa.tr$pbs.tr,method="spearman")
cor.test(pbs.sa.ty$pbs.sa,pbs.sa.ty$pbs.ty,method="spearman")
cor.test(pbs.tr.ty$pbs.tr,pbs.tr.ty$pbs.ty,method="spearman")

# All correlations are positive and significant at p < 2.2 x 10-16.
# r values range from 0.23 - 0.59.

# Specific Z and autosome relationships

pbs.er.gu.z <- setDT(pbs.er.gu)[CHROM %chin% list.z$V1]
pbs.er.ru.z <- setDT(pbs.er.ru)[CHROM %chin% list.z$V1]
pbs.er.sa.z <- setDT(pbs.er.sa)[CHROM %chin% list.z$V1]
pbs.er.tr.z <- setDT(pbs.er.tr)[CHROM %chin% list.z$V1]
pbs.er.ty.z <- setDT(pbs.er.ty)[CHROM %chin% list.z$V1]
pbs.gu.ru.z <- setDT(pbs.gu.ru)[CHROM %chin% list.z$V1]
pbs.gu.sa.z <- setDT(pbs.gu.sa)[CHROM %chin% list.z$V1]
pbs.gu.tr.z <- setDT(pbs.gu.tr)[CHROM %chin% list.z$V1]
pbs.gu.ty.z <- setDT(pbs.gu.ty)[CHROM %chin% list.z$V1]
pbs.ru.sa.z <- setDT(pbs.ru.sa)[CHROM %chin% list.z$V1]
pbs.ru.tr.z <- setDT(pbs.ru.tr)[CHROM %chin% list.z$V1]
pbs.ru.ty.z <- setDT(pbs.ru.ty)[CHROM %chin% list.z$V1]
pbs.sa.tr.z <- setDT(pbs.sa.tr)[CHROM %chin% list.z$V1]
pbs.sa.ty.z <- setDT(pbs.sa.ty)[CHROM %chin% list.z$V1]
pbs.tr.ty.z <- setDT(pbs.tr.ty)[CHROM %chin% list.z$V1]

cor.test(pbs.er.gu.z$pbs.er,pbs.er.gu.z$pbs.gu,method="spearman")
cor.test(pbs.er.ru.z$pbs.er,pbs.er.ru.z$pbs.ru,method="spearman")
cor.test(pbs.er.sa.z$pbs.er,pbs.er.sa.z$pbs.sa,method="spearman")
cor.test(pbs.er.tr.z$pbs.er,pbs.er.tr.z$pbs.tr,method="spearman")
cor.test(pbs.er.ty.z$pbs.er,pbs.er.ty.z$pbs.ty,method="spearman")
cor.test(pbs.gu.ru.z$pbs.gu,pbs.gu.ru.z$pbs.ru,method="spearman")
cor.test(pbs.gu.sa.z$pbs.gu,pbs.gu.sa.z$pbs.sa,method="spearman")
cor.test(pbs.gu.tr.z$pbs.gu,pbs.gu.tr.z$pbs.tr,method="spearman")
cor.test(pbs.gu.ty.z$pbs.gu,pbs.gu.ty.z$pbs.ty,method="spearman")
cor.test(pbs.ru.sa.z$pbs.ru,pbs.ru.sa.z$pbs.sa,method="spearman")
cor.test(pbs.ru.tr.z$pbs.ru,pbs.ru.tr.z$pbs.tr,method="spearman")
cor.test(pbs.ru.ty.z$pbs.ru,pbs.ru.ty.z$pbs.ty,method="spearman")
cor.test(pbs.sa.tr.z$pbs.sa,pbs.sa.tr.z$pbs.tr,method="spearman")
cor.test(pbs.sa.ty.z$pbs.sa,pbs.sa.ty.z$pbs.ty,method="spearman")
cor.test(pbs.tr.ty.z$pbs.tr,pbs.tr.ty.z$pbs.ty,method="spearman")

# All correlations are positive with p-values < 2.2 x 10-16.
# r values range from 0.419 - 0.78 (STRONGER than genome-wide).

pbs.er.gu.a <- setDT(pbs.er.gu)[CHROM %chin% list.a$V1]
pbs.er.ru.a <- setDT(pbs.er.ru)[CHROM %chin% list.a$V1]
pbs.er.sa.a <- setDT(pbs.er.sa)[CHROM %chin% list.a$V1]
pbs.er.tr.a <- setDT(pbs.er.tr)[CHROM %chin% list.a$V1]
pbs.er.ty.a <- setDT(pbs.er.ty)[CHROM %chin% list.a$V1]
pbs.gu.ru.a <- setDT(pbs.gu.ru)[CHROM %chin% list.a$V1]
pbs.gu.sa.a <- setDT(pbs.gu.sa)[CHROM %chin% list.a$V1]
pbs.gu.tr.a <- setDT(pbs.gu.tr)[CHROM %chin% list.a$V1]
pbs.gu.ty.a <- setDT(pbs.gu.ty)[CHROM %chin% list.a$V1]
pbs.ru.sa.a <- setDT(pbs.ru.sa)[CHROM %chin% list.a$V1]
pbs.ru.tr.a <- setDT(pbs.ru.tr)[CHROM %chin% list.a$V1]
pbs.ru.ty.a <- setDT(pbs.ru.ty)[CHROM %chin% list.a$V1]
pbs.sa.tr.a <- setDT(pbs.sa.tr)[CHROM %chin% list.a$V1]
pbs.sa.ty.a <- setDT(pbs.sa.ty)[CHROM %chin% list.a$V1]
pbs.tr.ty.a <- setDT(pbs.tr.ty)[CHROM %chin% list.a$V1]

cor.test(pbs.er.gu.a$pbs.er,pbs.er.gu.a$pbs.gu,method="spearman")
cor.test(pbs.er.ru.a$pbs.er,pbs.er.ru.a$pbs.ru,method="spearman")
cor.test(pbs.er.sa.a$pbs.er,pbs.er.sa.a$pbs.sa,method="spearman")
cor.test(pbs.er.tr.a$pbs.er,pbs.er.tr.a$pbs.tr,method="spearman")
cor.test(pbs.er.ty.a$pbs.er,pbs.er.ty.a$pbs.ty,method="spearman")
cor.test(pbs.gu.ru.a$pbs.gu,pbs.gu.ru.a$pbs.ru,method="spearman")
cor.test(pbs.gu.sa.a$pbs.gu,pbs.gu.sa.a$pbs.sa,method="spearman")
cor.test(pbs.gu.tr.a$pbs.gu,pbs.gu.tr.a$pbs.tr,method="spearman")
cor.test(pbs.gu.ty.a$pbs.gu,pbs.gu.ty.a$pbs.ty,method="spearman")
cor.test(pbs.ru.sa.a$pbs.ru,pbs.ru.sa.a$pbs.sa,method="spearman")
cor.test(pbs.ru.tr.a$pbs.ru,pbs.ru.tr.a$pbs.tr,method="spearman")
cor.test(pbs.ru.ty.a$pbs.ru,pbs.ru.ty.a$pbs.ty,method="spearman")
cor.test(pbs.sa.tr.a$pbs.sa,pbs.sa.tr.a$pbs.tr,method="spearman")
cor.test(pbs.sa.ty.a$pbs.sa,pbs.sa.ty.a$pbs.ty,method="spearman")
cor.test(pbs.tr.ty.a$pbs.tr,pbs.tr.ty.a$pbs.ty,method="spearman")

# All correlations are positive with p-values =< 2.29 x 10-5.
# r values range from 0.042 - 0.476 (Weaker than when Z is included).


### ------------------------------------------------------------------------
### Comparison with nucleotide diversity
### ------------------------------------------------------------------------

# Read in pixy pi data
pixy.pi <- read.table('../pixy_results/pixy_100kb_all_pi.txt',header=T)

# Subset subspecies pi data
er <- pixy.pi[which(pixy.pi$pop=='er'),]
gu <- pixy.pi[which(pixy.pi$pop=='gu'),]
ru <- pixy.pi[which(pixy.pi$pop=='ru'),]
sa <- pixy.pi[which(pixy.pi$pop=='sa'),]
tr <- pixy.pi[which(pixy.pi$pop=='tr'),]
ty <- pixy.pi[which(pixy.pi$pop=='ty'),]

# Merge pi and pbs data
er.pi.pbs <-merge(pbs.er.data,er,by.x=c("CHROM","BIN_START"),by.y=c("chromosome","window_pos_1")) 
gu.pi.pbs <-merge(pbs.gu.data,gu,by.x=c("CHROM","BIN_START"),by.y=c("chromosome","window_pos_1")) 
ru.pi.pbs <-merge(pbs.ru.data,ru,by.x=c("CHROM","BIN_START"),by.y=c("chromosome","window_pos_1")) 
sa.pi.pbs <-merge(pbs.sa.data,sa,by.x=c("CHROM","BIN_START"),by.y=c("chromosome","window_pos_1")) 
tr.pi.pbs <-merge(pbs.tr.data,tr,by.x=c("CHROM","BIN_START"),by.y=c("chromosome","window_pos_1")) 
ty.pi.pbs <-merge(pbs.ty.data,ty,by.x=c("CHROM","BIN_START"),by.y=c("chromosome","window_pos_1")) 

# Scatterplots
par(mfrow=c(2,3))
plot(er.pi.pbs$avg_pi,er.pi.pbs$pbs.er,pch=20,col=alpha('black',0.5))
plot(gu.pi.pbs$avg_pi,gu.pi.pbs$pbs.gu,pch=20,col=alpha('black',0.5))
plot(ru.pi.pbs$avg_pi,ru.pi.pbs$pbs.ru,pch=20,col=alpha('black',0.5))
plot(sa.pi.pbs$avg_pi,sa.pi.pbs$pbs.sa,pch=20,col=alpha('black',0.5))
plot(tr.pi.pbs$avg_pi,tr.pi.pbs$pbs.tr,pch=20,col=alpha('black',0.5))
plot(ty.pi.pbs$avg_pi,ty.pi.pbs$pbs.ty,pch=20,col=alpha('black',0.5))

# Spearman's rank order correlations
cor.test(er.pi.pbs$avg_pi,er.pi.pbs$pbs.er,method="spearman")
cor.test(gu.pi.pbs$avg_pi,gu.pi.pbs$pbs.gu,method="spearman")
cor.test(ru.pi.pbs$avg_pi,ru.pi.pbs$pbs.ru,method="spearman")
cor.test(sa.pi.pbs$avg_pi,sa.pi.pbs$pbs.sa,method="spearman")
cor.test(tr.pi.pbs$avg_pi,tr.pi.pbs$pbs.tr,method="spearman")
cor.test(ty.pi.pbs$avg_pi,ty.pi.pbs$pbs.ty,method="spearman")

# All correlations are negative with p-values < 2.2 x 10-16.
# r values range from -0.26 - -0.54.

# Specific Z and autosome relationships

er.pi.pbs.z <- setDT(er.pi.pbs)[CHROM %chin% list.z$V1]
gu.pi.pbs.z <- setDT(gu.pi.pbs)[CHROM %chin% list.z$V1]
ru.pi.pbs.z <- setDT(ru.pi.pbs)[CHROM %chin% list.z$V1]
sa.pi.pbs.z <- setDT(sa.pi.pbs)[CHROM %chin% list.z$V1]
tr.pi.pbs.z <- setDT(tr.pi.pbs)[CHROM %chin% list.z$V1]
ty.pi.pbs.z <- setDT(ty.pi.pbs)[CHROM %chin% list.z$V1]

par(mfrow=c(2,3))
plot(er.pi.pbs.z$PI,er.pi.pbs.z$pbs.er,pch=20,col=alpha('black',0.5))
plot(gu.pi.pbs.z$PI,gu.pi.pbs.z$pbs.gu,pch=20,col=alpha('black',0.5))
plot(ru.pi.pbs.z$PI,ru.pi.pbs.z$pbs.ru,pch=20,col=alpha('black',0.5))
plot(sa.pi.pbs.z$PI,sa.pi.pbs.z$pbs.sa,pch=20,col=alpha('black',0.5))
plot(tr.pi.pbs.z$PI,tr.pi.pbs.z$pbs.tr,pch=20,col=alpha('black',0.5))
plot(ty.pi.pbs.z$PI,ty.pi.pbs.z$pbs.ty,pch=20,col=alpha('black',0.5))

cor.test(er.pi.pbs.z$PI,er.pi.pbs.z$pbs.er,method="spearman")
cor.test(gu.pi.pbs.z$PI,gu.pi.pbs.z$pbs.gu,method="spearman")
cor.test(ru.pi.pbs.z$PI,ru.pi.pbs.z$pbs.ru,method="spearman")
cor.test(sa.pi.pbs.z$PI,sa.pi.pbs.z$pbs.sa,method="spearman")
cor.test(tr.pi.pbs.z$PI,tr.pi.pbs.z$pbs.tr,method="spearman")
cor.test(ty.pi.pbs.z$PI,ty.pi.pbs.z$pbs.ty,method="spearman")

# All Z-specific correlations negative with p-values < 2.2 x 10-16.
# r values range from -0.425 - -0.66 (STRONGER relationship than genome-wide).

er.pi.pbs.a <- setDT(er.pi.pbs)[CHROM %chin% list.a$V1]
gu.pi.pbs.a <- setDT(gu.pi.pbs)[CHROM %chin% list.a$V1]
ru.pi.pbs.a <- setDT(ru.pi.pbs)[CHROM %chin% list.a$V1]
sa.pi.pbs.a <- setDT(sa.pi.pbs)[CHROM %chin% list.a$V1]
tr.pi.pbs.a <- setDT(tr.pi.pbs)[CHROM %chin% list.a$V1]
ty.pi.pbs.a <- setDT(ty.pi.pbs)[CHROM %chin% list.a$V1]

par(mfrow=c(2,3))
plot(er.pi.pbs.a$PI,er.pi.pbs.a$pbs.er,pch=20,col=alpha('black',0.5))
plot(gu.pi.pbs.a$PI,gu.pi.pbs.a$pbs.gu,pch=20,col=alpha('black',0.5))
plot(ru.pi.pbs.a$PI,ru.pi.pbs.a$pbs.ru,pch=20,col=alpha('black',0.5))
plot(sa.pi.pbs.a$PI,sa.pi.pbs.a$pbs.sa,pch=20,col=alpha('black',0.5))
plot(tr.pi.pbs.a$PI,tr.pi.pbs.a$pbs.tr,pch=20,col=alpha('black',0.5))
plot(ty.pi.pbs.a$PI,ty.pi.pbs.a$pbs.ty,pch=20,col=alpha('black',0.5))

cor.test(er.pi.pbs.a$PI,er.pi.pbs.a$pbs.er,method="spearman")
cor.test(gu.pi.pbs.a$PI,gu.pi.pbs.a$pbs.gu,method="spearman")
cor.test(ru.pi.pbs.a$PI,ru.pi.pbs.a$pbs.ru,method="spearman")
cor.test(sa.pi.pbs.a$PI,sa.pi.pbs.a$pbs.sa,method="spearman")
cor.test(tr.pi.pbs.a$PI,tr.pi.pbs.a$pbs.tr,method="spearman")
cor.test(ty.pi.pbs.a$PI,ty.pi.pbs.a$pbs.ty,method="spearman")

# All Z-specific correlations negative with p-values < 0.001.
# r values range from -0.034 - -0.35 (WEAKER relationship than when Z is included).

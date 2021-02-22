############################################################################
# Barn Swallow introgression statistics
############################################################################

### Goal: quantify introgression between Barn Swallow subspecies across the
### genome and test the hypothesis that the Z Chromosome shows evidence of
### less introgression than autosomes.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('../abba-baba_results/')

library(data.table)
library(dplyr)
library(zoo)

### Read in Z-linked and autosomal scaffold lists---------------------------

list.z <- read.table('../processing_files/hirundo_rustica_scaffold_list.Z.txt',header=F)
list.a <- read.table('../processing_files/hirundo_rustica_scaffold_list.auto.txt',header=F)

### Read in data------------------------------------------------------------

rurt <- read.csv('abbababa.output.sa-ru-rt-sm.w2000.csv',header=T)
rurg <- read.csv('abbababa.output.sa-ru-rg-sm.w2000.csv',header=T)
gutg <- read.csv('abbababa.output.er-gu-tg-sm.w2000.csv',header=T)

### Split by autosome/Z-----------------------------------------------------

rurt.a <- setDT(rurt)[scaffold %chin% list.a$V1]
rurt.z <- setDT(rurt)[scaffold %chin% list.z$V1]

rurg.a <- setDT(rurg)[scaffold %chin% list.a$V1]
rurg.z <- setDT(rurg)[scaffold %chin% list.z$V1]

gutg.a <- setDT(gutg)[scaffold %chin% list.a$V1]
gutg.z <- setDT(gutg)[scaffold %chin% list.z$V1]

### Remove fd data where D < 0----------------------------------------------

# fd is meaningless between P2 and P3 if D is negative!

rurt.a$fd[rurt.a$D<0] <- NA
rurt.z$fd[rurt.z$D<0] <- NA

rurg.a$fd[rurg.a$D<0] <- NA
rurg.z$fd[rurg.z$D<0] <- NA

gutg.a$fd[gutg.a$D<0] <- NA
gutg.z$fd[gutg.z$D<0] <- NA

### Calculate means and SD--------------------------------------------------

mean(rurt.a$fd,na.rm=T)
mean(rurt.z$fd,na.rm=T)
sd(rurt.a$fd,na.rm=T)
sd(rurt.z$fd,na.rm=T)

mean(rurg.a$fd,na.rm=T)
mean(rurg.z$fd,na.rm=T)
sd(rurg.a$fd,na.rm=T)
sd(rurg.z$fd,na.rm=T)

mean(gutg.a$fd,na.rm=T)
mean(gutg.z$fd,na.rm=T)
sd(gutg.a$fd,na.rm=T)
sd(gutg.z$fd,na.rm=T)

### Plot distributions of fd------------------------------------------------

par(mfrow=c(1,3))
boxplot(rurt.a$fd,rurt.z$fd,pch=20,col=c('grey','seagreen'),ylab='fd',names=c('Auto','Z'))
boxplot(rurg.a$fd,rurg.z$fd,pch=20,col=c('grey','seagreen'),ylab='fd',names=c('Auto','Z'))
boxplot(gutg.a$fd,gutg.z$fd,pch=20,col=c('grey','seagreen'),ylab='fd',names=c('Auto','Z'))

### Mann-Whitney U tests----------------------------------------------------

wilcox.test(rurt.a$fd,rurt.z$fd)
wilcox.test(rurg.a$fd,rurg.z$fd)
wilcox.test(gutg.a$fd,gutg.z$fd)


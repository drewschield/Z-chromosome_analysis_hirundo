############################################################################
# Barn Swallow Z-linked and autosomal Tajima's D
############################################################################

### Goal: compare distributions of Tajima's D statistic between autosomes
### and the Z Chromosome for subspecies and hybrid zones in Barn Swallows.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('../tajimas_d_results/')

library(data.table)
library(dplyr)
library(scales)

### Read in Z-linked and autosomal scaffold lists---------------------------

list.z <- read.table('../processing_files/hirundo_rustica_scaffold_list.Z.txt',header=F)
list.a <- read.table('../processing_files/hirundo_rustica_scaffold_list.auto.txt',header=F)

### Read in data------------------------------------------------------------

## Male-only
# Subspecies
er.taj.male <- read.table('erythrogaster.male-only.100kb.Tajima.D',header=T)
gu.taj.male <- read.table('gutturalis.male-only.100kb.Tajima.D',header=T)
ru.taj.male <- read.table('rustica.male-only.100kb.Tajima.D',header=T)
sa.taj.male <- read.table('savignii.male-only.100kb.Tajima.D',header=T)
tr.taj.male <- read.table('transitiva.male-only.100kb.Tajima.D',header=T)
ty.taj.male <- read.table('tytleri.male-only.100kb.Tajima.D',header=T)
# Hybrids
rugu.taj.male <- read.table('rustica-gutturalis.male-only.100kb.Tajima.D',header=T)
ruty.taj.male <- read.table('rustica-tytleri.male-only.100kb.Tajima.D',header=T)
guty.taj.male <- read.table('tytleri-gutturalis.male-only.100kb.Tajima.D',header=T)

## Max-missing = 0.2
# Subspecies
er.taj.miss02 <- read.table('erythrogaster.miss02.100kb.Tajima.D',header=T)
gu.taj.miss02 <- read.table('gutturalis.miss02.100kb.Tajima.D',header=T)
ru.taj.miss02 <- read.table('rustica.miss02.100kb.Tajima.D',header=T)
sa.taj.miss02 <- read.table('savignii.miss02.100kb.Tajima.D',header=T)
tr.taj.miss02 <- read.table('transitiva.miss02.100kb.Tajima.D',header=T)
ty.taj.miss02 <- read.table('tytleri.miss02.100kb.Tajima.D',header=T)
# Hybrids
rugu.taj.miss02 <- read.table('rustica-gutturalis.miss02.100kb.Tajima.D',header=T)
ruty.taj.miss02 <- read.table('rustica-tytleri.miss02.100kb.Tajima.D',header=T)
guty.taj.miss02 <- read.table('tytleri-gutturalis.miss02.100kb.Tajima.D',header=T)

## Max-missing = 0.4
# Subspecies
er.taj.miss04 <- read.table('erythrogaster.miss04.100kb.Tajima.D',header=T)
gu.taj.miss04 <- read.table('gutturalis.miss04.100kb.Tajima.D',header=T)
ru.taj.miss04 <- read.table('rustica.miss04.100kb.Tajima.D',header=T)
sa.taj.miss04 <- read.table('savignii.miss04.100kb.Tajima.D',header=T)
tr.taj.miss04 <- read.table('transitiva.miss04.100kb.Tajima.D',header=T)
ty.taj.miss04 <- read.table('tytleri.miss04.100kb.Tajima.D',header=T)
# Hybrids
rugu.taj.miss04 <- read.table('rustica-gutturalis.miss04.100kb.Tajima.D',header=T)
ruty.taj.miss04 <- read.table('rustica-tytleri.miss04.100kb.Tajima.D',header=T)
guty.taj.miss04 <- read.table('tytleri-gutturalis.miss04.100kb.Tajima.D',header=T)

## Max-missing = 0.6
# Subspecies
er.taj.miss06 <- read.table('erythrogaster.miss06.100kb.Tajima.D',header=T)
gu.taj.miss06 <- read.table('gutturalis.miss06.100kb.Tajima.D',header=T)
ru.taj.miss06 <- read.table('rustica.miss06.100kb.Tajima.D',header=T)
sa.taj.miss06 <- read.table('savignii.miss06.100kb.Tajima.D',header=T)
tr.taj.miss06 <- read.table('transitiva.miss06.100kb.Tajima.D',header=T)
ty.taj.miss06 <- read.table('tytleri.miss06.100kb.Tajima.D',header=T)
# Hybrids
rugu.taj.miss06 <- read.table('rustica-gutturalis.miss06.100kb.Tajima.D',header=T)
ruty.taj.miss06 <- read.table('rustica-tytleri.miss06.100kb.Tajima.D',header=T)
guty.taj.miss06 <- read.table('tytleri-gutturalis.miss06.100kb.Tajima.D',header=T)

### Subset data by Z and autosome lists--------------------------------------

## Male-only
# Subspecies
er.taj.male.z <- setDT(er.taj.male)[CHROM %chin% list.z$V1]
er.taj.male.a <- setDT(er.taj.male)[CHROM %chin% list.a$V1]
gu.taj.male.z <- setDT(gu.taj.male)[CHROM %chin% list.z$V1]
gu.taj.male.a <- setDT(gu.taj.male)[CHROM %chin% list.a$V1]
ru.taj.male.z <- setDT(ru.taj.male)[CHROM %chin% list.z$V1]
ru.taj.male.a <- setDT(ru.taj.male)[CHROM %chin% list.a$V1]
sa.taj.male.z <- setDT(sa.taj.male)[CHROM %chin% list.z$V1]
sa.taj.male.a <- setDT(sa.taj.male)[CHROM %chin% list.a$V1]
tr.taj.male.z <- setDT(tr.taj.male)[CHROM %chin% list.z$V1]
tr.taj.male.a <- setDT(tr.taj.male)[CHROM %chin% list.a$V1]
ty.taj.male.z <- setDT(ty.taj.male)[CHROM %chin% list.z$V1]
ty.taj.male.a <- setDT(ty.taj.male)[CHROM %chin% list.a$V1]
# Hybrids
rugu.taj.male.z <- setDT(rugu.taj.male)[CHROM %chin% list.z$V1]
rugu.taj.male.a <- setDT(rugu.taj.male)[CHROM %chin% list.a$V1]
ruty.taj.male.z <- setDT(ruty.taj.male)[CHROM %chin% list.z$V1]
ruty.taj.male.a <- setDT(ruty.taj.male)[CHROM %chin% list.a$V1]
guty.taj.male.z <- setDT(guty.taj.male)[CHROM %chin% list.z$V1]
guty.taj.male.a <- setDT(guty.taj.male)[CHROM %chin% list.a$V1]

## Max-missing = 0.2
# Subspecies
er.taj.miss02.z <- setDT(er.taj.miss02)[CHROM %chin% list.z$V1]
er.taj.miss02.a <- setDT(er.taj.miss02)[CHROM %chin% list.a$V1]
gu.taj.miss02.z <- setDT(gu.taj.miss02)[CHROM %chin% list.z$V1]
gu.taj.miss02.a <- setDT(gu.taj.miss02)[CHROM %chin% list.a$V1]
ru.taj.miss02.z <- setDT(ru.taj.miss02)[CHROM %chin% list.z$V1]
ru.taj.miss02.a <- setDT(ru.taj.miss02)[CHROM %chin% list.a$V1]
sa.taj.miss02.z <- setDT(sa.taj.miss02)[CHROM %chin% list.z$V1]
sa.taj.miss02.a <- setDT(sa.taj.miss02)[CHROM %chin% list.a$V1]
tr.taj.miss02.z <- setDT(tr.taj.miss02)[CHROM %chin% list.z$V1]
tr.taj.miss02.a <- setDT(tr.taj.miss02)[CHROM %chin% list.a$V1]
ty.taj.miss02.z <- setDT(ty.taj.miss02)[CHROM %chin% list.z$V1]
ty.taj.miss02.a <- setDT(ty.taj.miss02)[CHROM %chin% list.a$V1]
# Hybrids
rugu.taj.miss02.z <- setDT(rugu.taj.miss02)[CHROM %chin% list.z$V1]
rugu.taj.miss02.a <- setDT(rugu.taj.miss02)[CHROM %chin% list.a$V1]
ruty.taj.miss02.z <- setDT(ruty.taj.miss02)[CHROM %chin% list.z$V1]
ruty.taj.miss02.a <- setDT(ruty.taj.miss02)[CHROM %chin% list.a$V1]
guty.taj.miss02.z <- setDT(guty.taj.miss02)[CHROM %chin% list.z$V1]
guty.taj.miss02.a <- setDT(guty.taj.miss02)[CHROM %chin% list.a$V1]

## Max-missing = 0.4
# Subspecies
er.taj.miss04.z <- setDT(er.taj.miss04)[CHROM %chin% list.z$V1]
er.taj.miss04.a <- setDT(er.taj.miss04)[CHROM %chin% list.a$V1]
gu.taj.miss04.z <- setDT(gu.taj.miss04)[CHROM %chin% list.z$V1]
gu.taj.miss04.a <- setDT(gu.taj.miss04)[CHROM %chin% list.a$V1]
ru.taj.miss04.z <- setDT(ru.taj.miss04)[CHROM %chin% list.z$V1]
ru.taj.miss04.a <- setDT(ru.taj.miss04)[CHROM %chin% list.a$V1]
sa.taj.miss04.z <- setDT(sa.taj.miss04)[CHROM %chin% list.z$V1]
sa.taj.miss04.a <- setDT(sa.taj.miss04)[CHROM %chin% list.a$V1]
tr.taj.miss04.z <- setDT(tr.taj.miss04)[CHROM %chin% list.z$V1]
tr.taj.miss04.a <- setDT(tr.taj.miss04)[CHROM %chin% list.a$V1]
ty.taj.miss04.z <- setDT(ty.taj.miss04)[CHROM %chin% list.z$V1]
ty.taj.miss04.a <- setDT(ty.taj.miss04)[CHROM %chin% list.a$V1]
# Hybrids
rugu.taj.miss04.z <- setDT(rugu.taj.miss04)[CHROM %chin% list.z$V1]
rugu.taj.miss04.a <- setDT(rugu.taj.miss04)[CHROM %chin% list.a$V1]
ruty.taj.miss04.z <- setDT(ruty.taj.miss04)[CHROM %chin% list.z$V1]
ruty.taj.miss04.a <- setDT(ruty.taj.miss04)[CHROM %chin% list.a$V1]
guty.taj.miss04.z <- setDT(guty.taj.miss04)[CHROM %chin% list.z$V1]
guty.taj.miss04.a <- setDT(guty.taj.miss04)[CHROM %chin% list.a$V1]

## Max-missing = 0.6
# Subspecies
er.taj.miss06.z <- setDT(er.taj.miss06)[CHROM %chin% list.z$V1]
er.taj.miss06.a <- setDT(er.taj.miss06)[CHROM %chin% list.a$V1]
gu.taj.miss06.z <- setDT(gu.taj.miss06)[CHROM %chin% list.z$V1]
gu.taj.miss06.a <- setDT(gu.taj.miss06)[CHROM %chin% list.a$V1]
ru.taj.miss06.z <- setDT(ru.taj.miss06)[CHROM %chin% list.z$V1]
ru.taj.miss06.a <- setDT(ru.taj.miss06)[CHROM %chin% list.a$V1]
sa.taj.miss06.z <- setDT(sa.taj.miss06)[CHROM %chin% list.z$V1]
sa.taj.miss06.a <- setDT(sa.taj.miss06)[CHROM %chin% list.a$V1]
tr.taj.miss06.z <- setDT(tr.taj.miss06)[CHROM %chin% list.z$V1]
tr.taj.miss06.a <- setDT(tr.taj.miss06)[CHROM %chin% list.a$V1]
ty.taj.miss06.z <- setDT(ty.taj.miss06)[CHROM %chin% list.z$V1]
ty.taj.miss06.a <- setDT(ty.taj.miss06)[CHROM %chin% list.a$V1]
# Hybrids
rugu.taj.miss06.z <- setDT(rugu.taj.miss06)[CHROM %chin% list.z$V1]
rugu.taj.miss06.a <- setDT(rugu.taj.miss06)[CHROM %chin% list.a$V1]
ruty.taj.miss06.z <- setDT(ruty.taj.miss06)[CHROM %chin% list.z$V1]
ruty.taj.miss06.a <- setDT(ruty.taj.miss06)[CHROM %chin% list.a$V1]
guty.taj.miss06.z <- setDT(guty.taj.miss06)[CHROM %chin% list.z$V1]
guty.taj.miss06.a <- setDT(guty.taj.miss06)[CHROM %chin% list.a$V1]

### Statistical summary------------------------------------------------------

## Miss04 (means of all)
mean(c(er.taj.miss04.a$TajimaD,
       gu.taj.miss04.a$TajimaD,
       ru.taj.miss04.a$TajimaD,
       sa.taj.miss04.a$TajimaD,
       tr.taj.miss04.a$TajimaD,
       ty.taj.miss04.a$TajimaD,
       rugu.taj.miss04.a$TajimaD,
       ruty.taj.miss04.a$TajimaD,
       guty.taj.miss04.a$TajimaD),
     na.rm=T)

sd(c(er.taj.miss04.a$TajimaD,
       gu.taj.miss04.a$TajimaD,
       ru.taj.miss04.a$TajimaD,
       sa.taj.miss04.a$TajimaD,
       tr.taj.miss04.a$TajimaD,
       ty.taj.miss04.a$TajimaD,
       rugu.taj.miss04.a$TajimaD,
       ruty.taj.miss04.a$TajimaD,
       guty.taj.miss04.a$TajimaD),
     na.rm=T)

mean(c(er.taj.miss04.z$TajimaD,
       gu.taj.miss04.z$TajimaD,
       ru.taj.miss04.z$TajimaD,
       sa.taj.miss04.z$TajimaD,
       tr.taj.miss04.z$TajimaD,
       ty.taj.miss04.z$TajimaD,
       rugu.taj.miss04.z$TajimaD,
       ruty.taj.miss04.z$TajimaD,
       guty.taj.miss04.z$TajimaD),
     na.rm=T)

sd(c(er.taj.miss04.z$TajimaD,
     gu.taj.miss04.z$TajimaD,
     ru.taj.miss04.z$TajimaD,
     sa.taj.miss04.z$TajimaD,
     tr.taj.miss04.z$TajimaD,
     ty.taj.miss04.z$TajimaD,
     rugu.taj.miss04.z$TajimaD,
     ruty.taj.miss04.z$TajimaD,
     guty.taj.miss04.z$TajimaD),
   na.rm=T)


## Male-only
# autosomes
mean(er.taj.male.a$TajimaD,na.rm=T)
mean(gu.taj.male.a$TajimaD,na.rm=T)
mean(ru.taj.male.a$TajimaD,na.rm=T)
mean(sa.taj.male.a$TajimaD,na.rm=T)
mean(tr.taj.male.a$TajimaD,na.rm=T)
mean(ty.taj.male.a$TajimaD,na.rm=T)
mean(rugu.taj.male.a$TajimaD,na.rm=T)
mean(ruty.taj.male.a$TajimaD,na.rm=T)
mean(guty.taj.male.a$TajimaD,na.rm=T)

sd(er.taj.male.a$TajimaD,na.rm=T)
sd(gu.taj.male.a$TajimaD,na.rm=T)
sd(ru.taj.male.a$TajimaD,na.rm=T)
sd(sa.taj.male.a$TajimaD,na.rm=T)
sd(tr.taj.male.a$TajimaD,na.rm=T)
sd(ty.taj.male.a$TajimaD,na.rm=T)
sd(rugu.taj.male.a$TajimaD,na.rm=T)
sd(ruty.taj.male.a$TajimaD,na.rm=T)
sd(guty.taj.male.a$TajimaD,na.rm=T)

# Z chromosome
mean(er.taj.male.z$TajimaD,na.rm=T)
mean(gu.taj.male.z$TajimaD,na.rm=T)
mean(ru.taj.male.z$TajimaD,na.rm=T)
mean(sa.taj.male.z$TajimaD,na.rm=T)
mean(tr.taj.male.z$TajimaD,na.rm=T)
mean(ty.taj.male.z$TajimaD,na.rm=T)
mean(rugu.taj.male.z$TajimaD,na.rm=T)
mean(ruty.taj.male.z$TajimaD,na.rm=T)
mean(guty.taj.male.z$TajimaD,na.rm=T)

sd(er.taj.male.z$TajimaD,na.rm=T)
sd(gu.taj.male.z$TajimaD,na.rm=T)
sd(ru.taj.male.z$TajimaD,na.rm=T)
sd(sa.taj.male.z$TajimaD,na.rm=T)
sd(tr.taj.male.z$TajimaD,na.rm=T)
sd(ty.taj.male.z$TajimaD,na.rm=T)
sd(rugu.taj.male.z$TajimaD,na.rm=T)
sd(ruty.taj.male.z$TajimaD,na.rm=T)
sd(guty.taj.male.z$TajimaD,na.rm=T)

## Max-missing = 0.2
# autosomes
mean(er.taj.miss02.a$TajimaD,na.rm=T)
mean(gu.taj.miss02.a$TajimaD,na.rm=T)
mean(ru.taj.miss02.a$TajimaD,na.rm=T)
mean(sa.taj.miss02.a$TajimaD,na.rm=T)
mean(tr.taj.miss02.a$TajimaD,na.rm=T)
mean(ty.taj.miss02.a$TajimaD,na.rm=T)
mean(rugu.taj.miss02.a$TajimaD,na.rm=T)
mean(ruty.taj.miss02.a$TajimaD,na.rm=T)
mean(guty.taj.miss02.a$TajimaD,na.rm=T)

sd(er.taj.miss02.a$TajimaD,na.rm=T)
sd(gu.taj.miss02.a$TajimaD,na.rm=T)
sd(ru.taj.miss02.a$TajimaD,na.rm=T)
sd(sa.taj.miss02.a$TajimaD,na.rm=T)
sd(tr.taj.miss02.a$TajimaD,na.rm=T)
sd(ty.taj.miss02.a$TajimaD,na.rm=T)
sd(rugu.taj.miss02.a$TajimaD,na.rm=T)
sd(ruty.taj.miss02.a$TajimaD,na.rm=T)
sd(guty.taj.miss02.a$TajimaD,na.rm=T)

# Z chromosome
mean(er.taj.miss02.z$TajimaD,na.rm=T)
mean(gu.taj.miss02.z$TajimaD,na.rm=T)
mean(ru.taj.miss02.z$TajimaD,na.rm=T)
mean(sa.taj.miss02.z$TajimaD,na.rm=T)
mean(tr.taj.miss02.z$TajimaD,na.rm=T)
mean(ty.taj.miss02.z$TajimaD,na.rm=T)
mean(rugu.taj.miss02.z$TajimaD,na.rm=T)
mean(ruty.taj.miss02.z$TajimaD,na.rm=T)
mean(guty.taj.miss02.z$TajimaD,na.rm=T)

sd(er.taj.miss02.z$TajimaD,na.rm=T)
sd(gu.taj.miss02.z$TajimaD,na.rm=T)
sd(ru.taj.miss02.z$TajimaD,na.rm=T)
sd(sa.taj.miss02.z$TajimaD,na.rm=T)
sd(tr.taj.miss02.z$TajimaD,na.rm=T)
sd(ty.taj.miss02.z$TajimaD,na.rm=T)
sd(rugu.taj.miss02.z$TajimaD,na.rm=T)
sd(ruty.taj.miss02.z$TajimaD,na.rm=T)
sd(guty.taj.miss02.z$TajimaD,na.rm=T)

## Max-missing = 0.4
# autosomes
mean(er.taj.miss04.a$TajimaD,na.rm=T)
mean(gu.taj.miss04.a$TajimaD,na.rm=T)
mean(ru.taj.miss04.a$TajimaD,na.rm=T)
mean(sa.taj.miss04.a$TajimaD,na.rm=T)
mean(tr.taj.miss04.a$TajimaD,na.rm=T)
mean(ty.taj.miss04.a$TajimaD,na.rm=T)
mean(rugu.taj.miss04.a$TajimaD,na.rm=T)
mean(ruty.taj.miss04.a$TajimaD,na.rm=T)
mean(guty.taj.miss04.a$TajimaD,na.rm=T)

sd(er.taj.miss04.a$TajimaD,na.rm=T)
sd(gu.taj.miss04.a$TajimaD,na.rm=T)
sd(ru.taj.miss04.a$TajimaD,na.rm=T)
sd(sa.taj.miss04.a$TajimaD,na.rm=T)
sd(tr.taj.miss04.a$TajimaD,na.rm=T)
sd(ty.taj.miss04.a$TajimaD,na.rm=T)
sd(rugu.taj.miss04.a$TajimaD,na.rm=T)
sd(ruty.taj.miss04.a$TajimaD,na.rm=T)
sd(guty.taj.miss04.a$TajimaD,na.rm=T)

# Z chromosome
mean(er.taj.miss04.z$TajimaD,na.rm=T)
mean(gu.taj.miss04.z$TajimaD,na.rm=T)
mean(ru.taj.miss04.z$TajimaD,na.rm=T)
mean(sa.taj.miss04.z$TajimaD,na.rm=T)
mean(tr.taj.miss04.z$TajimaD,na.rm=T)
mean(ty.taj.miss04.z$TajimaD,na.rm=T)
mean(rugu.taj.miss04.z$TajimaD,na.rm=T)
mean(ruty.taj.miss04.z$TajimaD,na.rm=T)
mean(guty.taj.miss04.z$TajimaD,na.rm=T)

sd(er.taj.miss04.z$TajimaD,na.rm=T)
sd(gu.taj.miss04.z$TajimaD,na.rm=T)
sd(ru.taj.miss04.z$TajimaD,na.rm=T)
sd(sa.taj.miss04.z$TajimaD,na.rm=T)
sd(tr.taj.miss04.z$TajimaD,na.rm=T)
sd(ty.taj.miss04.z$TajimaD,na.rm=T)
sd(rugu.taj.miss04.z$TajimaD,na.rm=T)
sd(ruty.taj.miss04.z$TajimaD,na.rm=T)
sd(guty.taj.miss04.z$TajimaD,na.rm=T)

## Max-missing = 0.6
# autosomes
mean(er.taj.miss06.a$TajimaD,na.rm=T)
mean(gu.taj.miss06.a$TajimaD,na.rm=T)
mean(ru.taj.miss06.a$TajimaD,na.rm=T)
mean(sa.taj.miss06.a$TajimaD,na.rm=T)
mean(tr.taj.miss06.a$TajimaD,na.rm=T)
mean(ty.taj.miss06.a$TajimaD,na.rm=T)
mean(rugu.taj.miss06.a$TajimaD,na.rm=T)
mean(ruty.taj.miss06.a$TajimaD,na.rm=T)
mean(guty.taj.miss06.a$TajimaD,na.rm=T)

sd(er.taj.miss06.a$TajimaD,na.rm=T)
sd(gu.taj.miss06.a$TajimaD,na.rm=T)
sd(ru.taj.miss06.a$TajimaD,na.rm=T)
sd(sa.taj.miss06.a$TajimaD,na.rm=T)
sd(tr.taj.miss06.a$TajimaD,na.rm=T)
sd(ty.taj.miss06.a$TajimaD,na.rm=T)
sd(rugu.taj.miss06.a$TajimaD,na.rm=T)
sd(ruty.taj.miss06.a$TajimaD,na.rm=T)
sd(guty.taj.miss06.a$TajimaD,na.rm=T)

# Z chromosome
mean(er.taj.miss06.z$TajimaD,na.rm=T)
mean(gu.taj.miss06.z$TajimaD,na.rm=T)
mean(ru.taj.miss06.z$TajimaD,na.rm=T)
mean(sa.taj.miss06.z$TajimaD,na.rm=T)
mean(tr.taj.miss06.z$TajimaD,na.rm=T)
mean(ty.taj.miss06.z$TajimaD,na.rm=T)
mean(rugu.taj.miss06.z$TajimaD,na.rm=T)
mean(ruty.taj.miss06.z$TajimaD,na.rm=T)
mean(guty.taj.miss06.z$TajimaD,na.rm=T)

sd(er.taj.miss06.z$TajimaD,na.rm=T)
sd(gu.taj.miss06.z$TajimaD,na.rm=T)
sd(ru.taj.miss06.z$TajimaD,na.rm=T)
sd(sa.taj.miss06.z$TajimaD,na.rm=T)
sd(tr.taj.miss06.z$TajimaD,na.rm=T)
sd(ty.taj.miss06.z$TajimaD,na.rm=T)
sd(rugu.taj.miss06.z$TajimaD,na.rm=T)
sd(ruty.taj.miss06.z$TajimaD,na.rm=T)
sd(guty.taj.miss06.z$TajimaD,na.rm=T)

### Plot density distributions------------------------------------------------

## Male-only
d.er.taj.male.a <- density(er.taj.male.a$TajimaD,na.rm=T)
d.gu.taj.male.a <- density(gu.taj.male.a$TajimaD,na.rm=T)
d.ru.taj.male.a <- density(ru.taj.male.a$TajimaD,na.rm=T)
d.sa.taj.male.a <- density(sa.taj.male.a$TajimaD,na.rm=T)
d.tr.taj.male.a <- density(tr.taj.male.a$TajimaD,na.rm=T)
d.ty.taj.male.a <- density(ty.taj.male.a$TajimaD,na.rm=T)
d.rugu.taj.male.a <- density(rugu.taj.male.a$TajimaD,na.rm=T)
d.guty.taj.male.a <- density(guty.taj.male.a$TajimaD,na.rm=T)
d.ruty.taj.male.a <- density(ruty.taj.male.a$TajimaD,na.rm=T)

d.er.taj.male.z <- density(er.taj.male.z$TajimaD,na.rm=T)
d.gu.taj.male.z <- density(gu.taj.male.z$TajimaD,na.rm=T)
d.ru.taj.male.z <- density(ru.taj.male.z$TajimaD,na.rm=T)
d.sa.taj.male.z <- density(sa.taj.male.z$TajimaD,na.rm=T)
d.tr.taj.male.z <- density(tr.taj.male.z$TajimaD,na.rm=T)
d.ty.taj.male.z <- density(ty.taj.male.z$TajimaD,na.rm=T)
d.rugu.taj.male.z <- density(rugu.taj.male.z$TajimaD,na.rm=T)
d.guty.taj.male.z <- density(guty.taj.male.z$TajimaD,na.rm=T)
d.ruty.taj.male.z <- density(ruty.taj.male.z$TajimaD,na.rm=T)

par(mfrow=c(3,3))
plot(density(er.taj.male.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D erythrogaster")
polygon(density(er.taj.male.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(er.taj.male.a$TajimaD,na.rm=T),y0=0,x1=median(er.taj.male.a$TajimaD,na.rm=T),y1=max(d.er.taj.male.a$y),col='grey')
lines(density(er.taj.male.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(er.taj.male.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(er.taj.male.z$TajimaD,na.rm=T),y0=0,x1=median(er.taj.male.z$TajimaD,na.rm=T),y1=max(d.er.taj.male.z$y),col='seagreen4')

plot(density(gu.taj.male.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis")
polygon(density(gu.taj.male.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(gu.taj.male.a$TajimaD,na.rm=T),y0=0,x1=median(gu.taj.male.a$TajimaD,na.rm=T),y1=max(d.gu.taj.male.a$y),col='grey')
lines(density(gu.taj.male.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(gu.taj.male.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(gu.taj.male.z$TajimaD,na.rm=T),y0=0,x1=median(gu.taj.male.z$TajimaD,na.rm=T),y1=max(d.gu.taj.male.z$y),col='seagreen4')

plot(density(ru.taj.male.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D rustica")
polygon(density(ru.taj.male.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(ru.taj.male.a$TajimaD,na.rm=T),y0=0,x1=median(ru.taj.male.a$TajimaD,na.rm=T),y1=max(d.ru.taj.male.a$y),col='grey')
lines(density(ru.taj.male.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ru.taj.male.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(ru.taj.male.z$TajimaD,na.rm=T),y0=0,x1=median(ru.taj.male.z$TajimaD,na.rm=T),y1=max(d.ru.taj.male.z$y),col='seagreen4')

plot(density(sa.taj.male.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D savignii")
polygon(density(sa.taj.male.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(sa.taj.male.a$TajimaD,na.rm=T),y0=0,x1=median(sa.taj.male.a$TajimaD,na.rm=T),y1=max(d.sa.taj.male.a$y),col='grey')
lines(density(sa.taj.male.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(sa.taj.male.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(sa.taj.male.z$TajimaD,na.rm=T),y0=0,x1=median(sa.taj.male.z$TajimaD,na.rm=T),y1=max(d.sa.taj.male.z$y),col='seagreen4')

plot(density(tr.taj.male.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D transitiva")
polygon(density(tr.taj.male.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(tr.taj.male.a$TajimaD,na.rm=T),y0=0,x1=median(tr.taj.male.a$TajimaD,na.rm=T),y1=max(d.tr.taj.male.a$y),col='grey')
lines(density(tr.taj.male.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(tr.taj.male.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(tr.taj.male.z$TajimaD,na.rm=T),y0=0,x1=median(tr.taj.male.z$TajimaD,na.rm=T),y1=max(d.tr.taj.male.z$y),col='seagreen4')

plot(density(ty.taj.male.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D tytleri")
polygon(density(ty.taj.male.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(ty.taj.male.a$TajimaD,na.rm=T),y0=0,x1=median(ty.taj.male.a$TajimaD,na.rm=T),y1=max(d.ty.taj.male.a$y),col='grey')
lines(density(ty.taj.male.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ty.taj.male.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(ty.taj.male.z$TajimaD,na.rm=T),y0=0,x1=median(ty.taj.male.z$TajimaD,na.rm=T),y1=max(d.ty.taj.male.z$y),col='seagreen4')

plot(density(guty.taj.male.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis x rustica")
polygon(density(guty.taj.male.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(guty.taj.male.a$TajimaD,na.rm=T),y0=0,x1=median(guty.taj.male.a$TajimaD,na.rm=T),y1=max(d.guty.taj.male.a$y),col='grey')
lines(density(guty.taj.male.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(guty.taj.male.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(guty.taj.male.z$TajimaD,na.rm=T),y0=0,x1=median(guty.taj.male.z$TajimaD,na.rm=T),y1=max(d.guty.taj.male.z$y),col='seagreen4')

plot(density(rugu.taj.male.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis x tytleri")
polygon(density(rugu.taj.male.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(rugu.taj.male.a$TajimaD,na.rm=T),y0=0,x1=median(rugu.taj.male.a$TajimaD,na.rm=T),y1=max(d.rugu.taj.male.a$y),col='grey')
lines(density(rugu.taj.male.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(rugu.taj.male.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(rugu.taj.male.z$TajimaD,na.rm=T),y0=0,x1=median(rugu.taj.male.z$TajimaD,na.rm=T),y1=max(d.rugu.taj.male.z$y),col='seagreen4')

plot(density(ruty.taj.male.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D rustica x tytleri")
polygon(density(ruty.taj.male.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(ruty.taj.male.a$TajimaD,na.rm=T),y0=0,x1=median(ruty.taj.male.a$TajimaD,na.rm=T),y1=max(d.ruty.taj.male.a$y),col='grey')
lines(density(ruty.taj.male.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ruty.taj.male.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(ruty.taj.male.z$TajimaD,na.rm=T),y0=0,x1=median(ruty.taj.male.z$TajimaD,na.rm=T),y1=max(d.ruty.taj.male.z$y),col='seagreen4')

## Max-missing = 0.2
d.er.taj.miss02.a <- density(er.taj.miss02.a$TajimaD,na.rm=T)
d.gu.taj.miss02.a <- density(gu.taj.miss02.a$TajimaD,na.rm=T)
d.ru.taj.miss02.a <- density(ru.taj.miss02.a$TajimaD,na.rm=T)
d.sa.taj.miss02.a <- density(sa.taj.miss02.a$TajimaD,na.rm=T)
d.tr.taj.miss02.a <- density(tr.taj.miss02.a$TajimaD,na.rm=T)
d.ty.taj.miss02.a <- density(ty.taj.miss02.a$TajimaD,na.rm=T)
d.rugu.taj.miss02.a <- density(rugu.taj.miss02.a$TajimaD,na.rm=T)
d.guty.taj.miss02.a <- density(guty.taj.miss02.a$TajimaD,na.rm=T)
d.ruty.taj.miss02.a <- density(ruty.taj.miss02.a$TajimaD,na.rm=T)

d.er.taj.miss02.z <- density(er.taj.miss02.z$TajimaD,na.rm=T)
d.gu.taj.miss02.z <- density(gu.taj.miss02.z$TajimaD,na.rm=T)
d.ru.taj.miss02.z <- density(ru.taj.miss02.z$TajimaD,na.rm=T)
d.sa.taj.miss02.z <- density(sa.taj.miss02.z$TajimaD,na.rm=T)
d.tr.taj.miss02.z <- density(tr.taj.miss02.z$TajimaD,na.rm=T)
d.ty.taj.miss02.z <- density(ty.taj.miss02.z$TajimaD,na.rm=T)
d.rugu.taj.miss02.z <- density(rugu.taj.miss02.z$TajimaD,na.rm=T)
d.guty.taj.miss02.z <- density(guty.taj.miss02.z$TajimaD,na.rm=T)
d.ruty.taj.miss02.z <- density(ruty.taj.miss02.z$TajimaD,na.rm=T)

par(mfrow=c(3,3))
plot(density(er.taj.miss02.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D erythrogaster")
polygon(density(er.taj.miss02.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(er.taj.miss02.a$TajimaD,na.rm=T),y0=0,x1=median(er.taj.miss02.a$TajimaD,na.rm=T),y1=max(d.er.taj.miss02.a$y),col='grey')
lines(density(er.taj.miss02.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(er.taj.miss02.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(er.taj.miss02.z$TajimaD,na.rm=T),y0=0,x1=median(er.taj.miss02.z$TajimaD,na.rm=T),y1=max(d.er.taj.miss02.z$y),col='seagreen4')

plot(density(gu.taj.miss02.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis")
polygon(density(gu.taj.miss02.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(gu.taj.miss02.a$TajimaD,na.rm=T),y0=0,x1=median(gu.taj.miss02.a$TajimaD,na.rm=T),y1=max(d.gu.taj.miss02.a$y),col='grey')
lines(density(gu.taj.miss02.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(gu.taj.miss02.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(gu.taj.miss02.z$TajimaD,na.rm=T),y0=0,x1=median(gu.taj.miss02.z$TajimaD,na.rm=T),y1=max(d.gu.taj.miss02.z$y),col='seagreen4')

plot(density(ru.taj.miss02.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D rustica")
polygon(density(ru.taj.miss02.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(ru.taj.miss02.a$TajimaD,na.rm=T),y0=0,x1=median(ru.taj.miss02.a$TajimaD,na.rm=T),y1=max(d.ru.taj.miss02.a$y),col='grey')
lines(density(ru.taj.miss02.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ru.taj.miss02.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(ru.taj.miss02.z$TajimaD,na.rm=T),y0=0,x1=median(ru.taj.miss02.z$TajimaD,na.rm=T),y1=max(d.ru.taj.miss02.z$y),col='seagreen4')

plot(density(sa.taj.miss02.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D savignii")
polygon(density(sa.taj.miss02.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(sa.taj.miss02.a$TajimaD,na.rm=T),y0=0,x1=median(sa.taj.miss02.a$TajimaD,na.rm=T),y1=max(d.sa.taj.miss02.a$y),col='grey')
lines(density(sa.taj.miss02.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(sa.taj.miss02.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(sa.taj.miss02.z$TajimaD,na.rm=T),y0=0,x1=median(sa.taj.miss02.z$TajimaD,na.rm=T),y1=max(d.sa.taj.miss02.z$y),col='seagreen4')

plot(density(tr.taj.miss02.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D transitiva")
polygon(density(tr.taj.miss02.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(tr.taj.miss02.a$TajimaD,na.rm=T),y0=0,x1=median(tr.taj.miss02.a$TajimaD,na.rm=T),y1=max(d.tr.taj.miss02.a$y),col='grey')
lines(density(tr.taj.miss02.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(tr.taj.miss02.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(tr.taj.miss02.z$TajimaD,na.rm=T),y0=0,x1=median(tr.taj.miss02.z$TajimaD,na.rm=T),y1=max(d.tr.taj.miss02.z$y),col='seagreen4')

plot(density(ty.taj.miss02.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D tytleri")
polygon(density(ty.taj.miss02.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(ty.taj.miss02.a$TajimaD,na.rm=T),y0=0,x1=median(ty.taj.miss02.a$TajimaD,na.rm=T),y1=max(d.ty.taj.miss02.a$y),col='grey')
lines(density(ty.taj.miss02.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ty.taj.miss02.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(ty.taj.miss02.z$TajimaD,na.rm=T),y0=0,x1=median(ty.taj.miss02.z$TajimaD,na.rm=T),y1=max(d.ty.taj.miss02.z$y),col='seagreen4')

plot(density(guty.taj.miss02.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis x rustica")
polygon(density(guty.taj.miss02.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(guty.taj.miss02.a$TajimaD,na.rm=T),y0=0,x1=median(guty.taj.miss02.a$TajimaD,na.rm=T),y1=max(d.guty.taj.miss02.a$y),col='grey')
lines(density(guty.taj.miss02.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(guty.taj.miss02.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(guty.taj.miss02.z$TajimaD,na.rm=T),y0=0,x1=median(guty.taj.miss02.z$TajimaD,na.rm=T),y1=max(d.guty.taj.miss02.z$y),col='seagreen4')

plot(density(rugu.taj.miss02.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis x tytleri")
polygon(density(rugu.taj.miss02.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(rugu.taj.miss02.a$TajimaD,na.rm=T),y0=0,x1=median(rugu.taj.miss02.a$TajimaD,na.rm=T),y1=max(d.rugu.taj.miss02.a$y),col='grey')
lines(density(rugu.taj.miss02.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(rugu.taj.miss02.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(rugu.taj.miss02.z$TajimaD,na.rm=T),y0=0,x1=median(rugu.taj.miss02.z$TajimaD,na.rm=T),y1=max(d.rugu.taj.miss02.z$y),col='seagreen4')

plot(density(ruty.taj.miss02.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D rustica x tytleri")
polygon(density(ruty.taj.miss02.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(ruty.taj.miss02.a$TajimaD,na.rm=T),y0=0,x1=median(ruty.taj.miss02.a$TajimaD,na.rm=T),y1=max(d.ruty.taj.miss02.a$y),col='grey')
lines(density(ruty.taj.miss02.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ruty.taj.miss02.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(ruty.taj.miss02.z$TajimaD,na.rm=T),y0=0,x1=median(ruty.taj.miss02.z$TajimaD,na.rm=T),y1=max(d.ruty.taj.miss02.z$y),col='seagreen4')

## Max-missing = 0.4
d.er.taj.miss04.a <- density(er.taj.miss04.a$TajimaD,na.rm=T)
d.gu.taj.miss04.a <- density(gu.taj.miss04.a$TajimaD,na.rm=T)
d.ru.taj.miss04.a <- density(ru.taj.miss04.a$TajimaD,na.rm=T)
d.sa.taj.miss04.a <- density(sa.taj.miss04.a$TajimaD,na.rm=T)
d.tr.taj.miss04.a <- density(tr.taj.miss04.a$TajimaD,na.rm=T)
d.ty.taj.miss04.a <- density(ty.taj.miss04.a$TajimaD,na.rm=T)
d.rugu.taj.miss04.a <- density(rugu.taj.miss04.a$TajimaD,na.rm=T)
d.guty.taj.miss04.a <- density(guty.taj.miss04.a$TajimaD,na.rm=T)
d.ruty.taj.miss04.a <- density(ruty.taj.miss04.a$TajimaD,na.rm=T)

d.er.taj.miss04.z <- density(er.taj.miss04.z$TajimaD,na.rm=T)
d.gu.taj.miss04.z <- density(gu.taj.miss04.z$TajimaD,na.rm=T)
d.ru.taj.miss04.z <- density(ru.taj.miss04.z$TajimaD,na.rm=T)
d.sa.taj.miss04.z <- density(sa.taj.miss04.z$TajimaD,na.rm=T)
d.tr.taj.miss04.z <- density(tr.taj.miss04.z$TajimaD,na.rm=T)
d.ty.taj.miss04.z <- density(ty.taj.miss04.z$TajimaD,na.rm=T)
d.rugu.taj.miss04.z <- density(rugu.taj.miss04.z$TajimaD,na.rm=T)
d.guty.taj.miss04.z <- density(guty.taj.miss04.z$TajimaD,na.rm=T)
d.ruty.taj.miss04.z <- density(ruty.taj.miss04.z$TajimaD,na.rm=T)

par(mfrow=c(3,3))
plot(density(er.taj.miss04.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D erythrogaster")
polygon(density(er.taj.miss04.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(er.taj.miss04.a$TajimaD,na.rm=T),y0=0,x1=median(er.taj.miss04.a$TajimaD,na.rm=T),y1=max(d.er.taj.miss04.a$y),col='grey')
lines(density(er.taj.miss04.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(er.taj.miss04.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(er.taj.miss04.z$TajimaD,na.rm=T),y0=0,x1=median(er.taj.miss04.z$TajimaD,na.rm=T),y1=max(d.er.taj.miss04.z$y),col='seagreen4')

plot(density(gu.taj.miss04.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis")
polygon(density(gu.taj.miss04.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(gu.taj.miss04.a$TajimaD,na.rm=T),y0=0,x1=median(gu.taj.miss04.a$TajimaD,na.rm=T),y1=max(d.gu.taj.miss04.a$y),col='grey')
lines(density(gu.taj.miss04.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(gu.taj.miss04.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(gu.taj.miss04.z$TajimaD,na.rm=T),y0=0,x1=median(gu.taj.miss04.z$TajimaD,na.rm=T),y1=max(d.gu.taj.miss04.z$y),col='seagreen4')

plot(density(ru.taj.miss04.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D rustica")
polygon(density(ru.taj.miss04.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(ru.taj.miss04.a$TajimaD,na.rm=T),y0=0,x1=median(ru.taj.miss04.a$TajimaD,na.rm=T),y1=max(d.ru.taj.miss04.a$y),col='grey')
lines(density(ru.taj.miss04.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ru.taj.miss04.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(ru.taj.miss04.z$TajimaD,na.rm=T),y0=0,x1=median(ru.taj.miss04.z$TajimaD,na.rm=T),y1=max(d.ru.taj.miss04.z$y),col='seagreen4')

plot(density(sa.taj.miss04.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D savignii")
polygon(density(sa.taj.miss04.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(sa.taj.miss04.a$TajimaD,na.rm=T),y0=0,x1=median(sa.taj.miss04.a$TajimaD,na.rm=T),y1=max(d.sa.taj.miss04.a$y),col='grey')
lines(density(sa.taj.miss04.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(sa.taj.miss04.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(sa.taj.miss04.z$TajimaD,na.rm=T),y0=0,x1=median(sa.taj.miss04.z$TajimaD,na.rm=T),y1=max(d.sa.taj.miss04.z$y),col='seagreen4')

plot(density(tr.taj.miss04.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D transitiva")
polygon(density(tr.taj.miss04.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(tr.taj.miss04.a$TajimaD,na.rm=T),y0=0,x1=median(tr.taj.miss04.a$TajimaD,na.rm=T),y1=max(d.tr.taj.miss04.a$y),col='grey')
lines(density(tr.taj.miss04.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(tr.taj.miss04.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(tr.taj.miss04.z$TajimaD,na.rm=T),y0=0,x1=median(tr.taj.miss04.z$TajimaD,na.rm=T),y1=max(d.tr.taj.miss04.z$y),col='seagreen4')

plot(density(ty.taj.miss04.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D tytleri")
polygon(density(ty.taj.miss04.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(ty.taj.miss04.a$TajimaD,na.rm=T),y0=0,x1=median(ty.taj.miss04.a$TajimaD,na.rm=T),y1=max(d.ty.taj.miss04.a$y),col='grey')
lines(density(ty.taj.miss04.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ty.taj.miss04.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(ty.taj.miss04.z$TajimaD,na.rm=T),y0=0,x1=median(ty.taj.miss04.z$TajimaD,na.rm=T),y1=max(d.ty.taj.miss04.z$y),col='seagreen4')

plot(density(guty.taj.miss04.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis x rustica")
polygon(density(guty.taj.miss04.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(guty.taj.miss04.a$TajimaD,na.rm=T),y0=0,x1=median(guty.taj.miss04.a$TajimaD,na.rm=T),y1=max(d.guty.taj.miss04.a$y),col='grey')
lines(density(guty.taj.miss04.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(guty.taj.miss04.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(guty.taj.miss04.z$TajimaD,na.rm=T),y0=0,x1=median(guty.taj.miss04.z$TajimaD,na.rm=T),y1=max(d.guty.taj.miss04.z$y),col='seagreen4')

plot(density(rugu.taj.miss04.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis x tytleri")
polygon(density(rugu.taj.miss04.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(rugu.taj.miss04.a$TajimaD,na.rm=T),y0=0,x1=median(rugu.taj.miss04.a$TajimaD,na.rm=T),y1=max(d.rugu.taj.miss04.a$y),col='grey')
lines(density(rugu.taj.miss04.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(rugu.taj.miss04.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(rugu.taj.miss04.z$TajimaD,na.rm=T),y0=0,x1=median(rugu.taj.miss04.z$TajimaD,na.rm=T),y1=max(d.rugu.taj.miss04.z$y),col='seagreen4')

plot(density(ruty.taj.miss04.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D rustica x tytleri")
polygon(density(ruty.taj.miss04.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(ruty.taj.miss04.a$TajimaD,na.rm=T),y0=0,x1=median(ruty.taj.miss04.a$TajimaD,na.rm=T),y1=max(d.ruty.taj.miss04.a$y),col='grey')
lines(density(ruty.taj.miss04.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ruty.taj.miss04.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(ruty.taj.miss04.z$TajimaD,na.rm=T),y0=0,x1=median(ruty.taj.miss04.z$TajimaD,na.rm=T),y1=max(d.ruty.taj.miss04.z$y),col='seagreen4')

## Max-missing = 0.6
d.er.taj.miss06.a <- density(er.taj.miss06.a$TajimaD,na.rm=T)
d.gu.taj.miss06.a <- density(gu.taj.miss06.a$TajimaD,na.rm=T)
d.ru.taj.miss06.a <- density(ru.taj.miss06.a$TajimaD,na.rm=T)
d.sa.taj.miss06.a <- density(sa.taj.miss06.a$TajimaD,na.rm=T)
d.tr.taj.miss06.a <- density(tr.taj.miss06.a$TajimaD,na.rm=T)
d.ty.taj.miss06.a <- density(ty.taj.miss06.a$TajimaD,na.rm=T)
d.rugu.taj.miss06.a <- density(rugu.taj.miss06.a$TajimaD,na.rm=T)
d.guty.taj.miss06.a <- density(guty.taj.miss06.a$TajimaD,na.rm=T)
d.ruty.taj.miss06.a <- density(ruty.taj.miss06.a$TajimaD,na.rm=T)

d.er.taj.miss06.z <- density(er.taj.miss06.z$TajimaD,na.rm=T)
d.gu.taj.miss06.z <- density(gu.taj.miss06.z$TajimaD,na.rm=T)
d.ru.taj.miss06.z <- density(ru.taj.miss06.z$TajimaD,na.rm=T)
d.sa.taj.miss06.z <- density(sa.taj.miss06.z$TajimaD,na.rm=T)
d.tr.taj.miss06.z <- density(tr.taj.miss06.z$TajimaD,na.rm=T)
d.ty.taj.miss06.z <- density(ty.taj.miss06.z$TajimaD,na.rm=T)
d.rugu.taj.miss06.z <- density(rugu.taj.miss06.z$TajimaD,na.rm=T)
d.guty.taj.miss06.z <- density(guty.taj.miss06.z$TajimaD,na.rm=T)
d.ruty.taj.miss06.z <- density(ruty.taj.miss06.z$TajimaD,na.rm=T)

par(mfrow=c(3,3))
plot(density(er.taj.miss06.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D erythrogaster")
polygon(density(er.taj.miss06.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(er.taj.miss06.a$TajimaD,na.rm=T),y0=0,x1=median(er.taj.miss06.a$TajimaD,na.rm=T),y1=max(d.er.taj.miss06.a$y),col='grey')
lines(density(er.taj.miss06.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(er.taj.miss06.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(er.taj.miss06.z$TajimaD,na.rm=T),y0=0,x1=median(er.taj.miss06.z$TajimaD,na.rm=T),y1=max(d.er.taj.miss06.z$y),col='seagreen4')

plot(density(gu.taj.miss06.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis")
polygon(density(gu.taj.miss06.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(gu.taj.miss06.a$TajimaD,na.rm=T),y0=0,x1=median(gu.taj.miss06.a$TajimaD,na.rm=T),y1=max(d.gu.taj.miss06.a$y),col='grey')
lines(density(gu.taj.miss06.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(gu.taj.miss06.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(gu.taj.miss06.z$TajimaD,na.rm=T),y0=0,x1=median(gu.taj.miss06.z$TajimaD,na.rm=T),y1=max(d.gu.taj.miss06.z$y),col='seagreen4')

plot(density(ru.taj.miss06.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D rustica")
polygon(density(ru.taj.miss06.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(ru.taj.miss06.a$TajimaD,na.rm=T),y0=0,x1=median(ru.taj.miss06.a$TajimaD,na.rm=T),y1=max(d.ru.taj.miss06.a$y),col='grey')
lines(density(ru.taj.miss06.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ru.taj.miss06.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(ru.taj.miss06.z$TajimaD,na.rm=T),y0=0,x1=median(ru.taj.miss06.z$TajimaD,na.rm=T),y1=max(d.ru.taj.miss06.z$y),col='seagreen4')

plot(density(sa.taj.miss06.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D savignii")
polygon(density(sa.taj.miss06.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(sa.taj.miss06.a$TajimaD,na.rm=T),y0=0,x1=median(sa.taj.miss06.a$TajimaD,na.rm=T),y1=max(d.sa.taj.miss06.a$y),col='grey')
lines(density(sa.taj.miss06.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(sa.taj.miss06.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(sa.taj.miss06.z$TajimaD,na.rm=T),y0=0,x1=median(sa.taj.miss06.z$TajimaD,na.rm=T),y1=max(d.sa.taj.miss06.z$y),col='seagreen4')

plot(density(tr.taj.miss06.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D transitiva")
polygon(density(tr.taj.miss06.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(tr.taj.miss06.a$TajimaD,na.rm=T),y0=0,x1=median(tr.taj.miss06.a$TajimaD,na.rm=T),y1=max(d.tr.taj.miss06.a$y),col='grey')
lines(density(tr.taj.miss06.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(tr.taj.miss06.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(tr.taj.miss06.z$TajimaD,na.rm=T),y0=0,x1=median(tr.taj.miss06.z$TajimaD,na.rm=T),y1=max(d.tr.taj.miss06.z$y),col='seagreen4')

plot(density(ty.taj.miss06.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D tytleri")
polygon(density(ty.taj.miss06.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(ty.taj.miss06.a$TajimaD,na.rm=T),y0=0,x1=median(ty.taj.miss06.a$TajimaD,na.rm=T),y1=max(d.ty.taj.miss06.a$y),col='grey')
lines(density(ty.taj.miss06.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ty.taj.miss06.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(ty.taj.miss06.z$TajimaD,na.rm=T),y0=0,x1=median(ty.taj.miss06.z$TajimaD,na.rm=T),y1=max(d.ty.taj.miss06.z$y),col='seagreen4')

plot(density(guty.taj.miss06.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis x rustica")
polygon(density(guty.taj.miss06.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(guty.taj.miss06.a$TajimaD,na.rm=T),y0=0,x1=median(guty.taj.miss06.a$TajimaD,na.rm=T),y1=max(d.guty.taj.miss06.a$y),col='grey')
lines(density(guty.taj.miss06.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(guty.taj.miss06.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(guty.taj.miss06.z$TajimaD,na.rm=T),y0=0,x1=median(guty.taj.miss06.z$TajimaD,na.rm=T),y1=max(d.guty.taj.miss06.z$y),col='seagreen4')

plot(density(rugu.taj.miss06.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis x tytleri")
polygon(density(rugu.taj.miss06.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(rugu.taj.miss06.a$TajimaD,na.rm=T),y0=0,x1=median(rugu.taj.miss06.a$TajimaD,na.rm=T),y1=max(d.rugu.taj.miss06.a$y),col='grey')
lines(density(rugu.taj.miss06.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(rugu.taj.miss06.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(rugu.taj.miss06.z$TajimaD,na.rm=T),y0=0,x1=median(rugu.taj.miss06.z$TajimaD,na.rm=T),y1=max(d.rugu.taj.miss06.z$y),col='seagreen4')

plot(density(ruty.taj.miss06.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D rustica x tytleri")
polygon(density(ruty.taj.miss06.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
segments(x0=median(ruty.taj.miss06.a$TajimaD,na.rm=T),y0=0,x1=median(ruty.taj.miss06.a$TajimaD,na.rm=T),y1=max(d.ruty.taj.miss06.a$y),col='grey')
lines(density(ruty.taj.miss06.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ruty.taj.miss06.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))
segments(x0=median(ruty.taj.miss06.z$TajimaD,na.rm=T),y0=0,x1=median(ruty.taj.miss06.z$TajimaD,na.rm=T),y1=max(d.ruty.taj.miss06.z$y),col='seagreen4')



#####################################OLDER BELOW-DER
# Without median lines

plot(density(gu.taj.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis")
polygon(density(gu.taj.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(gu.taj.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(gu.taj.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(ru.taj.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D rustica")
polygon(density(ru.taj.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(ru.taj.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ru.taj.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(sa.taj.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D savignii")
polygon(density(sa.taj.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(sa.taj.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(sa.taj.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(tr.taj.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D transitiva")
polygon(density(tr.taj.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(tr.taj.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(tr.taj.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(ty.taj.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D tytleri")
polygon(density(ty.taj.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(ty.taj.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ty.taj.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(rugu.taj.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis x rustica")
polygon(density(rugu.taj.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(rugu.taj.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(rugu.taj.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(guty.taj.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D gutturalis x tytleri")
polygon(density(guty.taj.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(guty.taj.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(guty.taj.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))

plot(density(ruty.taj.a$TajimaD,na.rm=T),col='lightgrey',xlim=c(-2,2.5),main=NA,xlab="Tajima's D rustica x tytleri")
polygon(density(ruty.taj.a$TajimaD,na.rm=T),col=alpha('lightgrey',0.5))
lines(density(ruty.taj.z$TajimaD,na.rm=T),col='seagreen4')
polygon(density(ruty.taj.z$TajimaD,na.rm=T),col=alpha('seagreen4',0.5))


### Basic plotting----------------------------------------------------------

par(mfrow=c(6,1))
plot(er.taj$TajimaD,type='l')
plot(gu.taj$TajimaD,type='l')
plot(ru.taj$TajimaD,type='l')
plot(sa.taj$TajimaD,type='l')
plot(tr.taj$TajimaD,type='l')
plot(ty.taj$TajimaD,type='l')

plot(rugu.taj$TajimaD,type='l')
plot(ruty.taj$TajimaD,type='l')
plot(guty.taj$TajimaD,type='l')

par(mfrow=c(1,1))
boxplot(er.taj.a$TajimaD,er.taj.z$TajimaD,
        gu.taj.a$TajimaD,gu.taj.z$TajimaD,
        ru.taj.a$TajimaD,ru.taj.z$TajimaD,
        sa.taj.a$TajimaD,sa.taj.z$TajimaD,
        tr.taj.a$TajimaD,tr.taj.z$TajimaD,
        ty.taj.a$TajimaD,ty.taj.z$TajimaD,
        guty.taj.a$TajimaD,guty.taj.z$TajimaD,
        rugu.taj.a$TajimaD,rugu.taj.z$TajimaD,
        ruty.taj.a$TajimaD,ruty.taj.z$TajimaD,
        pch=20,ylab="Tajima's D",
        col=c('seagreen4','orangered3'),outline=F,
        names=c('erythrogaster','','gutturalis','','rustica','','savignii','','transitiva','','tytleri','','gutturalis-tytleri','','rustica-gutturalis','','rustica-tytleri',''))

legend("topleft", legend=c("Autosomes", "Z Chromosome"),
       col=c("seagreen4", "orangered3"), lty=1,lwd=6, cex=1.2,bty = 'n')

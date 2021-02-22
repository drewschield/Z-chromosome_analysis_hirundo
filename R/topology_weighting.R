############################################################################
# Barn Swallow topology weighting
############################################################################

### Goal: quantify relative support for alternative topologies across the
### genome. As a demonstration of different relative support, we'll examine
### the triplet topology for rustica, savignii, and transtiva, and whether
### there is a greater amount of support for one topology over the others
### generally, and on the Z Chromosome, specifically.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('../twisst_results/')

install.packages('zoo')

library(data.table)
library(dplyr)
library(zoo)

### Read in Z-linked and autosomal scaffold lists---------------------------

list.z <- read.table('../processing_files/hirundo_rustica_scaffold_list.Z.txt',header=F)
list.a <- read.table('../processing_files/hirundo_rustica_scaffold_list.auto.txt',header=F)

### ------------------------------------------------------------------------
### rustica, savignii, transitiva triplet (triplet 1)
### ------------------------------------------------------------------------
### Topologies we are considering:
### 
###   topo 1:         topo 2:         topo 3:
###
###    /-sm             /-sm             /-sm      
###   |                |              --|          
### --|      /-ru    --|      /-ru      |   /-ru   
###   |   /-|          |   /-|           \-|       
###    \-|   \-sa       \-|   \-tr         |   /-sa
###      |                |                 \-|    
###       \-tr             \-sa                \-tr


### Read in data and split by chromosome-------------------------------------

w1kb <- read.table('output.ru-sa-tr.phyml_bionj.w1000.data.weights.chrom.txt',header=T)

w1kb.a <- setDT(w1kb)[scaffold %chin% list.a$V1]
w1kb.z <- setDT(w1kb)[scaffold %chin% list.z$V1]

### Calculate total topology weight and proportions of subtrees--------------

tot <- w1kb$topo1[1]+w1kb$topo2[1]+w1kb$topo3[1]

t1.a.prop <- w1kb.a$topo1/tot
t2.a.prop <- w1kb.a$topo2/tot
t3.a.prop <- w1kb.a$topo3/tot

t1.z.prop <- w1kb.z$topo1/tot
t2.z.prop <- w1kb.z$topo2/tot
t3.z.prop <- w1kb.z$topo3/tot

### Plot boxplots-------------------------------------------------------------

par(mfrow=c(1,2))

boxplot(t1.a.prop,t2.a.prop,t3.a.prop,pch=20,col=c('goldenrod','aquamarine4','skyblue'),names=c('t1','t2','t3'),ylab='Proportion of Subtrees',ylim=c(0,1),outline=F)
boxplot(t1.z.prop,t2.z.prop,t3.z.prop,pch=20,col=c('goldenrod','aquamarine4','skyblue'),names=c('t1','t2','t3'),ylab='Proportion of Subtrees',ylim=c(0,1),outline=F)

### Kruskal-Wallis tests (i.e., non-parametric ANOVA) for autosomes and Z-----

tr1.a <- tibble(t1.a.prop,t2.a.prop,t3.a.prop)
tr1.z <- tibble(t1.z.prop,t2.z.prop,t3.z.prop)

kruskal.test(t1.a.prop~t2.a.prop,data = tr1.a)
kruskal.test(t1.a.prop~t3.a.prop,data = tr1.a)
kruskal.test(t2.a.prop~t3.a.prop,data = tr1.a)

kruskal.test(t1.z.prop~t2.z.prop,data = tr1.z)
kruskal.test(t1.z.prop~t3.z.prop,data = tr1.z)
kruskal.test(t2.z.prop~t3.z.prop,data = tr1.z)

wilcox.test(t1.a.prop,t3.a.prop)
wilcox.test(t2.a.prop,t3.a.prop)

wilcox.test(t1.z.prop,t3.z.prop)
wilcox.test(t2.z.prop,t3.z.prop)

### Mann-Whitney U tests comparing topology support between autos and Z-------

wilcox.test(t3.a.prop,t3.z.prop)

wilcox.test(tr1.a$t1.a.prop,tr1.z$t1.z.prop)
wilcox.test(tr1.a$t2.a.prop,tr1.z$t2.z.prop)
wilcox.test(tr1.a$t3.a.prop,tr1.z$t3.z.prop)

# ----

t1.a.approx <- na.approx(t1.a.prop)
t2.a.approx <- na.approx(t2.a.prop)
t3.a.approx <- na.approx(t3.a.prop)

t1.z.approx <- na.approx(t1.z.prop)
t2.z.approx <- na.approx(t2.z.prop)
t3.z.approx <- na.approx(t3.z.prop)

t1s.a <- smooth.spline(t1.a.approx,spar=0.2)
t2s.a <- smooth.spline(t2.a.approx,spar=0.2)
t3s.a <- smooth.spline(t3.a.approx,spar=0.2)

t1s.z <- smooth.spline(t1.z.approx,spar=0.2)
t2s.z <- smooth.spline(t2.z.approx,spar=0.2)
t3s.z <- smooth.spline(t3.z.approx,spar=0.2)

par(mfrow=c(4,1))
plot(t1.z.prop,pch=20,col=alpha('red',0.25),ylim=c(0,1))
points(t2.z.prop,pch=20,col=alpha('lightblue',0.25))
points(t3.z.prop,pch=20,col=alpha('darkgreen',0.25))

lines(t1s.z,col='red',lwd=2)
lines(t2s.z,col='blue',lwd=2)
lines(t3s.z,col='darkgreen',lwd=2)


### ------------------------------------------------------------------------
### erythrogaster, gutturalis, tytleri triplet (triplet 2)
### ------------------------------------------------------------------------
### Topologies we are considering:
###                                                
###    /-sa             /-sa             /-sa      
###   |                |              --|          
### --|      /-er    --|      /-er      |   /-er   
###   |   /-|          |   /-|           \-|       
###    \-|   \-gu       \-|   \-ty         |   /-gu
###      |                |                 \-|    
###       \-ty             \-gu                \-ty
### 

### Read in data and split by chromosome-------------------------------------

w1kb <- read.table('output.er-gu-ty.phyml_bionj.w1000.data.weights.chrom.txt',header=T)

w1kb.a <- setDT(w1kb)[scaffold %chin% list.a$V1]
w1kb.z <- setDT(w1kb)[scaffold %chin% list.z$V1]

### Calculate total topology weight and proportions of subtrees--------------

tot <- w1kb$topo1[1]+w1kb$topo2[1]+w1kb$topo3[1]

t1.a.prop <- w1kb.a$topo1/tot
t2.a.prop <- w1kb.a$topo2/tot
t3.a.prop <- w1kb.a$topo3/tot

t1.z.prop <- w1kb.z$topo1/tot
t2.z.prop <- w1kb.z$topo2/tot
t3.z.prop <- w1kb.z$topo3/tot

### Plot boxplots-------------------------------------------------------------

par(mfrow=c(1,2))

boxplot(t1.a.prop,t2.a.prop,t3.a.prop,pch=20,col=c('goldenrod','aquamarine4','skyblue'),names=c('t1','t2','t3'),ylab='Proportion of Subtrees',ylim=c(0,1),outline=F)
boxplot(t1.z.prop,t2.z.prop,t3.z.prop,pch=20,col=c('goldenrod','aquamarine4','skyblue'),names=c('t1','t2','t3'),ylab='Proportion of Subtrees',ylim=c(0,1),outline=F)

### Kruskal-Wallis tests (i.e., non-parametric ANOVA) for autosomes and Z-----

tr1.a <- tibble(t1.a.prop,t2.a.prop,t3.a.prop)
tr1.z <- tibble(t1.z.prop,t2.z.prop,t3.z.prop)

kruskal.test(t1.a.prop~t2.a.prop,data = tr1.a)
kruskal.test(t1.a.prop~t3.a.prop,data = tr1.a)
kruskal.test(t2.a.prop~t3.a.prop,data = tr1.a)

kruskal.test(t1.z.prop~t2.z.prop,data = tr1.z)
kruskal.test(t1.z.prop~t3.z.prop,data = tr1.z)
kruskal.test(t2.z.prop~t3.z.prop,data = tr1.z)

### Mann-Whitney U tests comparing topology support between autos and Z-------

wilcox.test(tr1.a$t1.a.prop,tr1.z$t1.z.prop)
wilcox.test(tr1.a$t2.a.prop,tr1.z$t2.z.prop)
wilcox.test(tr1.a$t3.a.prop,tr1.z$t3.z.prop)

# test Z support for t2 over t1 and t3
wilcox.test(tr1.z$t2.z.prop,tr1.z$t1.z.prop)
wilcox.test(tr1.z$t2.z.prop,tr1.z$t3.z.prop)

### ------------------------------------------------------------------------
### savignii, transitiva, erythrogaster triplet (triplet 3)
### ------------------------------------------------------------------------
### Topologies we are considering:
###                                                
###    /-sm             /-sm             /-sm      
###   |                |              --|          
### --|      /-sa    --|      /-sa      |   /-sa   
###   |   /-|          |   /-|           \-|       
###    \-|   \-tr       \-|   \-er         |   /-tr
###      |                |                 \-|    
###       \-er             \-tr                \-er
### 

### Read in data and split by chromosome-------------------------------------

w1kb <- read.table('output.sa-tr-er.phyml_bionj.w1000.data.weights.chrom.txt',header=T)

w1kb.a <- setDT(w1kb)[scaffold %chin% list.a$V1]
w1kb.z <- setDT(w1kb)[scaffold %chin% list.z$V1]

### Calculate total topology weight and proportions of subtrees--------------

tot <- w1kb$topo1[1]+w1kb$topo2[1]+w1kb$topo3[1]

t1.a.prop <- w1kb.a$topo1/tot
t2.a.prop <- w1kb.a$topo2/tot
t3.a.prop <- w1kb.a$topo3/tot

t1.z.prop <- w1kb.z$topo1/tot
t2.z.prop <- w1kb.z$topo2/tot
t3.z.prop <- w1kb.z$topo3/tot

### Plot boxplots-------------------------------------------------------------

par(mfrow=c(1,2))

boxplot(t1.a.prop,t2.a.prop,t3.a.prop,pch=20,col=c('goldenrod','aquamarine4','skyblue'),names=c('t1','t2','t3'),ylab='Proportion of Subtrees',ylim=c(0,1),outline=F)
boxplot(t1.z.prop,t2.z.prop,t3.z.prop,pch=20,col=c('goldenrod','aquamarine4','skyblue'),names=c('t1','t2','t3'),ylab='Proportion of Subtrees',ylim=c(0,1),outline=F)

### Kruskal-Wallis tests (i.e., non-parametric ANOVA) for autosomes and Z-----

tr1.a <- tibble(t1.a.prop,t2.a.prop,t3.a.prop)
tr1.z <- tibble(t1.z.prop,t2.z.prop,t3.z.prop)

kruskal.test(t1.a.prop~t2.a.prop,data = tr1.a)
kruskal.test(t1.a.prop~t3.a.prop,data = tr1.a)
kruskal.test(t2.a.prop~t3.a.prop,data = tr1.a)

kruskal.test(t1.z.prop~t2.z.prop,data = tr1.z)
kruskal.test(t1.z.prop~t3.z.prop,data = tr1.z)
kruskal.test(t2.z.prop~t3.z.prop,data = tr1.z)

### Mann-Whitney U tests comparing topology support between autos and Z-------

wilcox.test(tr1.a$t1.a.prop,tr1.z$t1.z.prop)
wilcox.test(tr1.a$t2.a.prop,tr1.z$t2.z.prop)
wilcox.test(tr1.a$t3.a.prop,tr1.z$t3.z.prop)

# t1 versus t2 and t3 for autosomes
wilcox.test(tr1.a$t1.a.prop,tr1.a$t2.a.prop)
wilcox.test(tr1.a$t1.a.prop,tr1.a$t3.a.prop)

# t1 versus t2 and t3 for Z
wilcox.test(tr1.z$t1.z.prop,tr1.z$t2.z.prop)
wilcox.test(tr1.z$t1.z.prop,tr1.z$t3.z.prop)

# Z t1 versus auto t1
wilcox.test(tr1.a$t1.a.prop,tr1.z$t1.z.prop)

### ------------------------------------------------------------------------
### erythrogaster, tytleri, savignii triplet (triplet 4)
### ------------------------------------------------------------------------
### Topologies we are considering:
###                                                
###    /-sm             /-sm             /-sm      
###   |                |              --|          
### --|      /-er    --|      /-er      |   /-er   
###   |   /-|          |   /-|           \-|       
###    \-|   \-ty       \-|   \-sa         |   /-ty
###      |                |                 \-|    
###       \-sa             \-ty                \-sa
### 

### Read in data and split by chromosome-------------------------------------

w1kb <- read.table('output.er-ty-sa.phyml_bionj.w1000.data.weights.chrom.txt',header=T)

w1kb.a <- setDT(w1kb)[scaffold %chin% list.a$V1]
w1kb.z <- setDT(w1kb)[scaffold %chin% list.z$V1]

### Calculate total topology weight and proportions of subtrees--------------

tot <- w1kb$topo1[1]+w1kb$topo2[1]+w1kb$topo3[1]

t1.a.prop <- w1kb.a$topo1/tot
t2.a.prop <- w1kb.a$topo2/tot
t3.a.prop <- w1kb.a$topo3/tot

t1.z.prop <- w1kb.z$topo1/tot
t2.z.prop <- w1kb.z$topo2/tot
t3.z.prop <- w1kb.z$topo3/tot

### Plot boxplots-------------------------------------------------------------

par(mfrow=c(1,2))

boxplot(t1.a.prop,t2.a.prop,t3.a.prop,pch=20,col=c('goldenrod','aquamarine4','skyblue'),names=c('t1','t2','t3'),ylab='Proportion of Subtrees',ylim=c(0,1),outline=F)
boxplot(t1.z.prop,t2.z.prop,t3.z.prop,pch=20,col=c('goldenrod','aquamarine4','skyblue'),names=c('t1','t2','t3'),ylab='Proportion of Subtrees',ylim=c(0,1),outline=F)

### Kruskal-Wallis tests (i.e., non-parametric ANOVA) for autosomes and Z-----

tr1.a <- tibble(t1.a.prop,t2.a.prop,t3.a.prop)
tr1.z <- tibble(t1.z.prop,t2.z.prop,t3.z.prop)

kruskal.test(t1.a.prop~t2.a.prop,data = tr1.a)
kruskal.test(t1.a.prop~t3.a.prop,data = tr1.a)
kruskal.test(t2.a.prop~t3.a.prop,data = tr1.a)

kruskal.test(t1.z.prop~t2.z.prop,data = tr1.z)
kruskal.test(t1.z.prop~t3.z.prop,data = tr1.z)
kruskal.test(t2.z.prop~t3.z.prop,data = tr1.z)

### Mann-Whitney U tests comparing topology support between autos and Z-------

wilcox.test(tr1.a$t1.a.prop,tr1.z$t1.z.prop)
wilcox.test(tr1.a$t2.a.prop,tr1.z$t2.z.prop)
wilcox.test(tr1.a$t3.a.prop,tr1.z$t3.z.prop)

# t1 versus t2 and t3 for autosomes
wilcox.test(tr1.a$t1.a.prop,tr1.a$t2.a.prop)
wilcox.test(tr1.a$t1.a.prop,tr1.a$t3.a.prop)

# t1 versus t2 and t3 for Z
wilcox.test(tr1.z$t1.z.prop,tr1.z$t2.z.prop)
wilcox.test(tr1.z$t1.z.prop,tr1.z$t3.z.prop)

# Z t1 versus auto t1
wilcox.test(tr1.a$t1.a.prop,tr1.z$t1.z.prop)

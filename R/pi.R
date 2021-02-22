############################################################################
# Barn Swallow Z-linked and autosomal nucleotide diversity
############################################################################

### Goal: compare distributions of nucleotide diversity among Barn Swallow
### subspecies, populations, and hybrid zones.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('../pixy_results')

library(data.table)
library(dplyr)
library(scales)

### Read in Z-linked and autosomal scaffold lists---------------------------

list.z <- read.table('../processing_files/hirundo_rustica_scaffold_list.Z.txt',header=F)
list.a <- read.table('../processing_files/hirundo_rustica_scaffold_list.auto.txt',header=F)

### ------------------------------------------------------------------------
### Nucleotide diversity statistics
### ------------------------------------------------------------------------

### Read in data------------------------------------------------------------

## Concatenated pixy output for subspecies/hybrids and specific localities
pixy.pi <- read.table('pixy_100kb_all_pi.txt',header=T)
pixy.pi.loc <- read.table('pixy_localities_100kb_all_pi.txt',header=T)

# Subspecies
er <- pixy.pi[which(pixy.pi$pop=='er'),]
gu <- pixy.pi[which(pixy.pi$pop=='gu'),]
ru <- pixy.pi[which(pixy.pi$pop=='ru'),]
sa <- pixy.pi[which(pixy.pi$pop=='sa'),]
tr <- pixy.pi[which(pixy.pi$pop=='tr'),]
ty <- pixy.pi[which(pixy.pi$pop=='ty'),]
# Hybrids
rugu <- pixy.pi[which(pixy.pi$pop=='rg'),]
ruty <- pixy.pi[which(pixy.pi$pop=='rt'),]
guty <- pixy.pi[which(pixy.pi$pop=='tg'),]
# Within subspecies
gu_boatu <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_boatu'),]
gu_changchun <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_changchun'),]
gu_changsha <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_changsha'),]
gu_hainan <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_hainan'),]
gu_harbin <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_harbin'),]
gu_hokkaido <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_hokkaido'),]
gu_nanning <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_nanning'),]
gu_qinhuangdao <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_qinhuangdao'),]
gu_qiqihar <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_qiqihar'),]
gu_shuang <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_shuang'),]
gu_tokyo <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_tokyo'),]
gu_xian <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_xian'),]
gu_zhengzhou <- pixy.pi.loc[which(pixy.pi.loc$pop=='gu_zhengzhou'),]
ru_agadir <- pixy.pi.loc[which(pixy.pi.loc$pop=='ru_agadir'),]
ru_beni <- pixy.pi.loc[which(pixy.pi.loc$pop=='ru_beni'),]
ru_durgun <- pixy.pi.loc[which(pixy.pi.loc$pop=='ru_durgun'),]
ru_khovd <- pixy.pi.loc[which(pixy.pi.loc$pop=='ru_khovd'),]
ru_krasnoyarsk <- pixy.pi.loc[which(pixy.pi.loc$pop=='ru_krasnoyarsk'),]
ru_marrakech <- pixy.pi.loc[which(pixy.pi.loc$pop=='ru_marrakech'),]
ru_moscow <- pixy.pi.loc[which(pixy.pi.loc$pop=='ru_moscow'),]
ru_urumqi <- pixy.pi.loc[which(pixy.pi.loc$pop=='ru_urumqi'),]
ru_yekaterinburg <- pixy.pi.loc[which(pixy.pi.loc$pop=='ru_yekaterinburg'),]
ty_kytyleek <- pixy.pi.loc[which(pixy.pi.loc$pop=='ty_kytyleek'),]
ty_malamolevo <- pixy.pi.loc[which(pixy.pi.loc$pop=='ty_malamolevo'),]
ty_zakaltoose <- pixy.pi.loc[which(pixy.pi.loc$pop=='ty_zakaltoose'),]

## Subset data by Z and autosome lists--------------------------------------

# Subspecies
er.z <- setDT(er)[chromosome %chin% list.z$V1]
er.a <- setDT(er)[chromosome %chin% list.a$V1]
gu.z <- setDT(gu)[chromosome %chin% list.z$V1]
gu.a <- setDT(gu)[chromosome %chin% list.a$V1]
ru.z <- setDT(ru)[chromosome %chin% list.z$V1]
ru.a <- setDT(ru)[chromosome %chin% list.a$V1]
sa.z <- setDT(sa)[chromosome %chin% list.z$V1]
sa.a <- setDT(sa)[chromosome %chin% list.a$V1]
tr.z <- setDT(tr)[chromosome %chin% list.z$V1]
tr.a <- setDT(tr)[chromosome %chin% list.a$V1]
ty.z <- setDT(ty)[chromosome %chin% list.z$V1]
ty.a <- setDT(ty)[chromosome %chin% list.a$V1]
# Hybrids
rugu.z <- setDT(rugu)[chromosome %chin% list.z$V1]
rugu.a <- setDT(rugu)[chromosome %chin% list.a$V1]
ruty.z <- setDT(ruty)[chromosome %chin% list.z$V1]
ruty.a <- setDT(ruty)[chromosome %chin% list.a$V1]
guty.z <- setDT(guty)[chromosome %chin% list.z$V1]
guty.a <- setDT(guty)[chromosome %chin% list.a$V1]
# Within subspecies
gu_boatu.a <- setDT(gu_boatu)[chromosome %chin% list.a$V1]
gu_changchun.a <- setDT(gu_changchun)[chromosome %chin% list.a$V1]
gu_changsha.a <- setDT(gu_changsha)[chromosome %chin% list.a$V1]
gu_hainan.a <- setDT(gu_hainan)[chromosome %chin% list.a$V1]
gu_harbin.a <- setDT(gu_harbin)[chromosome %chin% list.a$V1]
gu_hokkaido.a <- setDT(gu_hokkaido)[chromosome %chin% list.a$V1]
gu_nanning.a <- setDT(gu_nanning)[chromosome %chin% list.a$V1]
gu_qinhuangdao.a <- setDT(gu_qinhuangdao)[chromosome %chin% list.a$V1]
gu_qiqihar.a <- setDT(gu_qiqihar)[chromosome %chin% list.a$V1]
gu_shuang.a <- setDT(gu_shuang)[chromosome %chin% list.a$V1]
gu_tokyo.a <- setDT(gu_tokyo)[chromosome %chin% list.a$V1]
gu_xian.a <- setDT(gu_xian)[chromosome %chin% list.a$V1]
gu_zhengzhou.a <- setDT(gu_zhengzhou)[chromosome %chin% list.a$V1]
ru_agadir.a <- setDT(ru_agadir)[chromosome %chin% list.a$V1]
ru_beni.a <- setDT(ru_beni)[chromosome %chin% list.a$V1]
ru_durgun.a <- setDT(ru_durgun)[chromosome %chin% list.a$V1]
ru_khovd.a <- setDT(ru_khovd)[chromosome %chin% list.a$V1]
ru_krasnoyarsk.a <- setDT(ru_krasnoyarsk)[chromosome %chin% list.a$V1]
ru_marrakech.a <- setDT(ru_marrakech)[chromosome %chin% list.a$V1]
ru_moscow.a <- setDT(ru_moscow)[chromosome %chin% list.a$V1]
ru_urumqi.a <- setDT(ru_urumqi)[chromosome %chin% list.a$V1]
ru_yekaterinburg.a <- setDT(ru_yekaterinburg)[chromosome %chin% list.a$V1]
ty_kytyleek.a <- setDT(ty_kytyleek)[chromosome %chin% list.a$V1]
ty_malamolevo.a <- setDT(ty_malamolevo)[chromosome %chin% list.a$V1]
ty_zakaltoose.a <- setDT(ty_zakaltoose)[chromosome %chin% list.a$V1]

gu_boatu.z <- setDT(gu_boatu)[chromosome %chin% list.z$V1]
gu_changchun.z <- setDT(gu_changchun)[chromosome %chin% list.z$V1]
gu_changsha.z <- setDT(gu_changsha)[chromosome %chin% list.z$V1]
gu_hainan.z <- setDT(gu_hainan)[chromosome %chin% list.z$V1]
gu_harbin.z <- setDT(gu_harbin)[chromosome %chin% list.z$V1]
gu_hokkaido.z <- setDT(gu_hokkaido)[chromosome %chin% list.z$V1]
gu_nanning.z <- setDT(gu_nanning)[chromosome %chin% list.z$V1]
gu_qinhuangdao.z <- setDT(gu_qinhuangdao)[chromosome %chin% list.z$V1]
gu_qiqihar.z <- setDT(gu_qiqihar)[chromosome %chin% list.z$V1]
gu_shuang.z <- setDT(gu_shuang)[chromosome %chin% list.z$V1]
gu_tokyo.z <- setDT(gu_tokyo)[chromosome %chin% list.z$V1]
gu_xian.z <- setDT(gu_xian)[chromosome %chin% list.z$V1]
gu_zhengzhou.z <- setDT(gu_zhengzhou)[chromosome %chin% list.z$V1]
ru_agadir.z <- setDT(ru_agadir)[chromosome %chin% list.z$V1]
ru_beni.z <- setDT(ru_beni)[chromosome %chin% list.z$V1]
ru_durgun.z <- setDT(ru_durgun)[chromosome %chin% list.z$V1]
ru_khovd.z <- setDT(ru_khovd)[chromosome %chin% list.z$V1]
ru_krasnoyarsk.z <- setDT(ru_krasnoyarsk)[chromosome %chin% list.z$V1]
ru_marrakech.z <- setDT(ru_marrakech)[chromosome %chin% list.z$V1]
ru_moscow.z <- setDT(ru_moscow)[chromosome %chin% list.z$V1]
ru_urumqi.z <- setDT(ru_urumqi)[chromosome %chin% list.z$V1]
ru_yekaterinburg.z <- setDT(ru_yekaterinburg)[chromosome %chin% list.z$V1]
ty_kytyleek.z <- setDT(ty_kytyleek)[chromosome %chin% list.z$V1]
ty_malamolevo.z <- setDT(ty_malamolevo)[chromosome %chin% list.z$V1]
ty_zakaltoose.z <- setDT(ty_zakaltoose)[chromosome %chin% list.z$V1]

### Basic plotting----------------------------------------------------------

# Main barplot
par(mfrow=c(1,1))
boxplot(er.a$avg_pi,er.z$avg_pi,
        gu.a$avg_pi,gu.z$avg_pi,
        ru.a$avg_pi,ru.z$avg_pi,
        sa.a$avg_pi,sa.z$avg_pi,
        tr.a$avg_pi,tr.z$avg_pi,
        ty.a$avg_pi,ty.z$avg_pi,
        guty.a$avg_pi,guty.z$avg_pi,
        rugu.a$avg_pi,rugu.z$avg_pi,
        ruty.a$avg_pi,ruty.z$avg_pi,
        pch=20,ylab='Nucleotide Diversity',ylim=c(0,0.012),
        col=c('seagreen4','orangered3'),outline=F,
        names=c('erythrogaster','','gutturalis','','rustica','','savignii','','transitiva','','tytleri','','gutturalis-tytleri','','rustica-gutturalis','','rustica-tytleri',''))

legend("topleft", legend=c("Autosomes", "Z Chromosome"),
       col=c("seagreen4", "orangered3"), lty=1,lwd=6, cex=1.2,bty = 'n')

# Hybrids
boxplot(rugu.a$PI,rugu.z$PI,ruty.a$PI,ruty.z$PI,guty.a$PI,guty.z$PI,col=c('seagreen4','orangered3'),ylim=c(0,0.006))

# Make a more attractive GGplot version with jittering

df.er.a <- data.frame(group = 'erA', value = er.a$PI)
df.er.z <- data.frame(group = 'erZ', value = er.z$PI)
df.gu.a <- data.frame(group = 'guA', value = gu.a$PI)
df.gu.z <- data.frame(group = 'guZ', value = gu.z$PI)
df.ru.a <- data.frame(group = 'ruA', value = ru.a$PI)
df.ru.z <- data.frame(group = 'ruZ', value = ru.z$PI)
df.sa.a <- data.frame(group = 'saA', value = sa.a$PI)
df.sa.z <- data.frame(group = 'saZ', value = sa.z$PI)
df.tr.a <- data.frame(group = 'trA', value = tr.a$PI)
df.tr.z <- data.frame(group = 'trZ', value = tr.z$PI)
df.ty.a <- data.frame(group = 'tyA', value = ty.a$PI)
df.ty.z <- data.frame(group = 'tyZ', value = ty.z$PI)
df.guty.a <- data.frame(group = 'gutyA', value = guty.a$PI)
df.guty.z <- data.frame(group = 'gutyZ', value = guty.z$PI)
df.rugu.a <- data.frame(group = 'ruguA', value = rugu.a$PI)
df.rugu.z <- data.frame(group = 'ruguZ', value = rugu.z$PI)
df.ruty.a <- data.frame(group = 'rutyA', value = ruty.a$PI)
df.ruty.z <- data.frame(group = 'rutyZ', value = ruty.z$PI)

plot.data <- rbind(df.er.a,
                   df.er.z,
                   df.gu.a,
                   df.gu.z,
                   df.ru.a,
                   df.ru.z,
                   df.sa.a,
                   df.sa.z,
                   df.tr.a,
                   df.tr.z,
                   df.ty.a,
                   df.ty.z,
                   df.guty.a,
                   df.guty.z,
                   df.rugu.a,
                   df.rugu.z,
                   df.ruty.a,
                   df.ruty.z)

### Summary statistics------------------------------------------------------

# Mean Pi for autosomes and the Z Chromosome
# Subspecies
mean(er.a$avg_pi)
mean(er.z$avg_pi)
mean(gu.a$avg_pi)
mean(gu.z$avg_pi)
mean(ru.a$avg_pi)
mean(ru.z$avg_pi)
mean(sa.a$avg_pi)
mean(sa.z$avg_pi)
mean(tr.a$avg_pi)
mean(tr.z$avg_pi)
mean(ty.a$avg_pi)
mean(ty.z$avg_pi)
# Hybrids
mean(rugu.a$avg_pi)
mean(rugu.z$avg_pi)
mean(ruty.a$avg_pi)
mean(ruty.z$avg_pi)
mean(guty.a$avg_pi)
mean(guty.z$avg_pi)
# Within subspecies
mean.within.a <- c()
mean.within.a <- c(mean.within.a, mean(gu_boatu.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(gu_changchun.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(gu_changsha.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(gu_hainan.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(gu_harbin.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(gu_hokkaido.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(gu_nanning.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(gu_qinhuangdao.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(gu_qiqihar.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(gu_shuang.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(gu_tokyo.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(gu_xian.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(gu_zhengzhou.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(ru_agadir.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(ru_beni.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(ru_durgun.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(ru_khovd.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(ru_krasnoyarsk.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(ru_marrakech.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(ru_moscow.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(ru_urumqi.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(ru_yekaterinburg.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(ty_kytyleek.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(ty_malamolevo.a$avg_pi));
mean.within.a <- c(mean.within.a, mean(ty_zakaltoose.a$avg_pi));

mean.within.z <- c()
mean.within.z <- c(mean.within.z, mean(gu_boatu.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(gu_changchun.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(gu_changsha.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(gu_hainan.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(gu_harbin.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(gu_hokkaido.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(gu_nanning.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(gu_qinhuangdao.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(gu_qiqihar.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(gu_shuang.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(gu_tokyo.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(gu_xian.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(gu_zhengzhou.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(ru_agadir.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(ru_beni.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(ru_durgun.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(ru_khovd.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(ru_krasnoyarsk.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(ru_marrakech.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(ru_moscow.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(ru_urumqi.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(ru_yekaterinburg.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(ty_kytyleek.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(ty_malamolevo.z$avg_pi));
mean.within.z <- c(mean.within.z, mean(ty_zakaltoose.z$avg_pi));

mean.within <- data.frame(mean.within.a,mean.within.z)

# Ratio of means

pi.mean <- read.table('../pixy_results/mean_pi_summary.txt',header=T)

par(mfrow=c(1,2))
# Plot at different x-axis scales
scatt <- palette(c('red','purple','black'))
palette()
color = rep(NA, length=length(pi.mean$category))
color[which(pi.mean$category=="subspecies")] = "red"
color[which(pi.mean$category=="hybrid")] = "purple"
color[which(pi.mean$category=="within")] = "black"
plot(pi.mean$mean.pi.auto,pi.mean$mean.pi.z.corr,pch=20,xlim=c(0,0.006),ylim=c(0,0.006),xlab='Autosomal Pi',ylab='Z Chromosome Pi',col=color)
theo.pi.a <- c(0,0.0015,0.003,0.0045,0.006)
theo.pi.z <- theo.pi.a*0.75
theo.pi.var <- theo.pi.a*0.5625
theo.pi.var2 <- theo.pi.a*1.125
abline(lm(theo.pi.z~theo.pi.a),col='black',lty=1)
abline(lm(theo.pi.var~theo.pi.a),col='orange',lty=1)
abline(lm(theo.pi.var2~theo.pi.a),col='skyblue',lty=1)

cor.test(pi.mean$mean.pi.auto,pi.mean$mean.pi.z.corr,method = "spearman")

pi.mean.mean <- (pi.mean$mean.pi.z.corr + pi.mean$mean.pi.auto)/2
pi.za <- pi.mean$mean.pi.z.corr/pi.mean$mean.pi.auto

plot(pi.mean.mean,pi.za,pch=20,ylim=c(0,1.2),col=color,xlab='Mean Pi',ylab='PiZA')
abline(h=0.75,col='black')
abline(h=0.5625,col='orange')
abline(h=1.125,col='skyblue')

### Two-sample t-tests between autosomes and Z------------------------------

t.test(er.a$avg_pi,er.z$avg_pi)
t.test(gu.a$avg_pi,gu.z$avg_pi)
t.test(ru.a$avg_pi,ru.z$avg_pi)
t.test(sa.a$avg_pi,sa.z$avg_pi)
t.test(tr.a$avg_pi,tr.z$avg_pi)
t.test(ty.a$avg_pi,ty.z$avg_pi)

t.test(guty.a$avg_pi,guty.z$avg_pi)
t.test(rugu.a$avg_pi,rugu.z$avg_pi)
t.test(ruty.a$avg_pi,ruty.z$avg_pi)

t.test(gu_boatu.a$avg_pi,gu_boatu.z$avg_pi)
t.test(gu_changchun.a$avg_pi,gu_changchun.z$avg_pi)
t.test(gu_changsha.a$avg_pi,gu_changsha.z$avg_pi)
t.test(gu_hainan.a$avg_pi,gu_hainan.z$avg_pi)
t.test(gu_harbin.a$avg_pi,gu_harbin.z$avg_pi)
t.test(gu_hokkaido.a$avg_pi,gu_hokkaido.z$avg_pi)
t.test(gu_nanning.a$avg_pi,gu_nanning.z$avg_pi)
t.test(gu_qinhuangdao.a$avg_pi,gu_qinhuangdao.z$avg_pi)
t.test(gu_qiqihar.a$avg_pi,gu_qiqihar.z$avg_pi)
t.test(gu_shuang.a$avg_pi,gu_shuang.z$avg_pi)
t.test(gu_tokyo.a$avg_pi,gu_tokyo.z$avg_pi)
t.test(gu_xian.a$avg_pi,gu_xian.z$avg_pi)
t.test(gu_zhengzhou.a$avg_pi,gu_zhengzhou.z$avg_pi)
t.test(ru_agadir.a$avg_pi,ru_agadir.z$avg_pi)
t.test(ru_beni.a$avg_pi,ru_beni.z$avg_pi)
t.test(ru_durgun.a$avg_pi,ru_durgun.z$avg_pi)
t.test(ru_khovd.a$avg_pi,ru_khovd.z$avg_pi)
t.test(ru_krasnoyarsk.a$avg_pi,ru_krasnoyarsk.z$avg_pi)
t.test(ru_marrakech.a$avg_pi,ru_marrakech.z$avg_pi)
t.test(ru_moscow.a$avg_pi,ru_moscow.z$avg_pi)
t.test(ru_urumqi.a$avg_pi,ru_urumqi.z$avg_pi)
t.test(ru_yekaterinburg.a$avg_pi,ru_yekaterinburg.z$avg_pi)
t.test(ty_kytyleek.a$avg_pi,ty_kytyleek.z$avg_pi)
t.test(ty_malamolevo.a$avg_pi,ty_malamolevo.z$avg_pi)
t.test(ty_zakaltoose.a$avg_pi,ty_zakaltoose.z$avg_pi)

### Correlations between subspecies-----------------------------------------

pi.comp.erty <-merge(er,ty,by=c("chromosome","window_pos_1")) 
pi.comp.ergu <-merge(er,gu,by=c("chromosome","window_pos_1")) 
pi.comp.erru <-merge(er,ru,by=c("chromosome","window_pos_1")) 
pi.comp.ersa <-merge(er,sa,by=c("chromosome","window_pos_1")) 
pi.comp.ertr <-merge(er,tr,by=c("chromosome","window_pos_1"))
pi.comp.guru <-merge(gu,ru,by=c("chromosome","window_pos_1"))
pi.comp.gusa <-merge(gu,sa,by=c("chromosome","window_pos_1"))
pi.comp.gutr <-merge(gu,tr,by=c("chromosome","window_pos_1"))
pi.comp.guty <-merge(gu,ty,by=c("chromosome","window_pos_1"))
pi.comp.rusa <-merge(ru,sa,by=c("chromosome","window_pos_1"))
pi.comp.rutr <-merge(ru,tr,by=c("chromosome","window_pos_1"))
pi.comp.ruty <-merge(ru,ty,by=c("chromosome","window_pos_1"))
pi.comp.satr <-merge(sa,tr,by=c("chromosome","window_pos_1"))
pi.comp.saty <-merge(sa,ty,by=c("chromosome","window_pos_1"))
pi.comp.trty <-merge(tr,ty,by=c("chromosome","window_pos_1"))

par(mfrow=c(5,3))
plot(pi.comp.ergu$avg_pi.x,pi.comp.ergu$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='erythrogaster',ylab='gutturalis')
plot(pi.comp.erru$avg_pi.x,pi.comp.erru$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='erythrogaster',ylab='rustica')
plot(pi.comp.ersa$avg_pi.x,pi.comp.ersa$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='erythrogaster',ylab='savignii')
plot(pi.comp.ertr$avg_pi.x,pi.comp.ertr$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='erythrogaster',ylab='transitiva')
plot(pi.comp.erty$avg_pi.x,pi.comp.erty$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='erythrogaster',ylab='tytleri')
plot(pi.comp.guru$avg_pi.x,pi.comp.guru$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='gutturalis',ylab='rustica')
plot(pi.comp.gusa$avg_pi.x,pi.comp.gusa$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='gutturalis',ylab='savignii')
plot(pi.comp.gutr$avg_pi.x,pi.comp.gutr$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='gutturalis',ylab='transitiva')
plot(pi.comp.guty$avg_pi.x,pi.comp.guty$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='gutturalis',ylab='tytleri')
plot(pi.comp.rusa$avg_pi.x,pi.comp.rusa$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='rustica',ylab='savignii')
plot(pi.comp.rutr$avg_pi.x,pi.comp.rutr$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='rustica',ylab='transitiva')
plot(pi.comp.ruty$avg_pi.x,pi.comp.ruty$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='rustica',ylab='transitiva')
plot(pi.comp.satr$avg_pi.x,pi.comp.satr$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='savignii',ylab='transitiva')
plot(pi.comp.saty$avg_pi.x,pi.comp.saty$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='savginii',ylab='tytleri')
plot(pi.comp.trty$avg_pi.x,pi.comp.trty$avg_pi.y,pch=20,col=alpha('black',0.25),xlab='transitiva',ylab='tytleri')

cor.test(pi.comp.erty$avg_pi.x,pi.comp.erty$avg_pi.y)
cor.test(pi.comp.ergu$avg_pi.x,pi.comp.ergu$avg_pi.y)
cor.test(pi.comp.erru$avg_pi.x,pi.comp.erru$avg_pi.y)
cor.test(pi.comp.ersa$avg_pi.x,pi.comp.ersa$avg_pi.y)
cor.test(pi.comp.ertr$avg_pi.x,pi.comp.ertr$avg_pi.y)
cor.test(pi.comp.guru$avg_pi.x,pi.comp.guru$avg_pi.y)
cor.test(pi.comp.gusa$avg_pi.x,pi.comp.gusa$avg_pi.y)
cor.test(pi.comp.gutr$avg_pi.x,pi.comp.gutr$avg_pi.y)
cor.test(pi.comp.guty$avg_pi.x,pi.comp.guty$avg_pi.y)
cor.test(pi.comp.rusa$avg_pi.x,pi.comp.rusa$avg_pi.y)
cor.test(pi.comp.rutr$avg_pi.x,pi.comp.rutr$avg_pi.y)
cor.test(pi.comp.ruty$avg_pi.x,pi.comp.ruty$avg_pi.y)
cor.test(pi.comp.satr$avg_pi.x,pi.comp.satr$avg_pi.y)
cor.test(pi.comp.saty$avg_pi.x,pi.comp.saty$avg_pi.y)
cor.test(pi.comp.trty$avg_pi.x,pi.comp.trty$avg_pi.y)

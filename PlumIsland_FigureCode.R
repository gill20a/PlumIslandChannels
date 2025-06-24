setwd("~/Dropbox/Williams/Students/Armstrong, Michael/")
depth<-read.csv("PIE_Excel_Simple_Field_MA.csv")
head(depth)
setwd("~/Dropbox/Williams/Students/Armstrong, Michael/Figure Set/")
sat<-read.csv("20230428_SatalliteData.csv")
slumps<-read.csv("20230601_Slumps.csv")
library(readr)
library(stats4)
library(bbmle)
library(mltools)
library(Metrics)
library(nlme)

LLS = function(y,k){
  Mhat=1*exp(-k*xNA$t) #creates a vector of=length to obs of preds
  ssq = sum((Mhat - xNA$Mt)^2)
  sigma = sqrt(ssq/length(xNA$Mt))
  return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
}




head(sat)
unique(sat$Stream)
# install.packages("RColorBrewer")
library(RColorBrewer)
colnames(sat)
max(sat$W_2019)
dev.off()
par(mfcol=c(2,2), mar=c(4,4,0.5,0.5), mgp=c(2,0.8,0))
# tiff("20230621_2_Figure2.tiff", width = 6, height = 6, units = 'in', res = 300)

pick1<-sat$Stream=="SC"
plot(sat[pick1,"Dist"], sat[pick1,"W_2019"], 
     xlab="Distance Downstream (m)", ylab="Width (m)", col="#4b9645", xlim=c(0,5000))
pick1<-sat$Stream=="WC"
points(sat[pick1,"Dist"], sat[pick1,"W_2019"], 
       xlab="Distance Downstream (m)", ylab="Depth (m)", col="#67B9b4")
pick1<-sat$Stream=="CHC"
points(sat[pick1,"Dist"], sat[pick1,"W_2019"], 
     xlab="Distance Downstream (m)", ylab="Depth (m)", col="#904596")
pick1<-sat$Stream=="SWC"
points(sat[pick1,"Dist"], sat[pick1,"W_2019"], 
       xlab="Distance Downstream (m)", ylab="Depth (m)", col="grey50")
pick1<-sat$Stream=="SHC"
points(sat[pick1,"Dist"], sat[pick1,"W_2019"], 
       xlab="Distance Downstream (m)", ylab="Depth (m)", col="#D95F02")
mtext("a", side=1, adj=0.95,line=-1.5, font=2)
legend("topleft", c("SC", "WC", "CHC", "SWC", "SHC"), bty="n",
       col=c("#4b9645", "#67B9b4", "#904596","grey50",  "#D95F02"), cex=0.9, pch=16)
fit.lme<-lme(W_2019~Dist, random = ~1|Stream, data=sat, na.action = na.omit)
summary(fit.lme)

library(MuMIn)
r.squaredGLMM(fit.lme)

summary(fit.lme)
anova(fit.lme)

library(nlme)
##


pick1<-depth$Stream=="CHC"
plot(depth[pick1,"Prop_Downstream"], depth[pick1,"Depth"], 
       xlab="Proportion Downstream", ylab="Depth (m)", col="#904596", ylim=c(1.5,3.5), xlim=c(0,1))
fit.CHC<-lm(depth[pick1,"Depth"]~depth[pick1,"Prop_Downstream"])
summary(fit.CHC)
pick1<-depth$Stream=="SC"
points(depth[pick1,"Prop_Downstream"], depth[pick1,"Depth"], 
     xlab="Proportion Downstream", ylab="Depth (m)", col="#4b9645")
fit.SC<-lm(depth[pick1,"Depth"]~depth[pick1,"Prop_Downstream"], na.action = na.omit)
summary(fit.SC)

fit<-lm(depth[,"Depth"]~depth[,"Prop_Downstream"], na.action = na.omit)
summary(fit)
head(depth)
unique(depth$Stream)
fit.lme<-lme(Depth~Prop_Downstream, random = ~1|Stream, data=depth, na.action = na.omit)
summary(fit.lme)
anova(fit.lme)
r.squaredGLMM(fit.lme)
coef<-c(1.907764,1.398446 )
abline(coef, lwd=2, col='grey80')
mtext("b", side=1, adj=0.95,line=-1.5, font=2)
legend("topleft", c("SC","CHC"), bty="n",
       col=c("#4b9645",  "#904596"), cex=0.9, pch=16)

##


jose<-read.csv("Jose Calcs.csv")
hist(log(sat[,"Curvepos"]), breaks=30)
sat$Log_Curvepos<-log(sat[,"Curvepos"])
plot(sat[,"Proportion_Downstream"], sat[,"Curvepos"], ylim=c(0,0.05),
     xlab="Proportion Downstream", ylab=expression("Curvature (rad m"^-1*")"), col=rgb(0.33, 0.33, 0.33, 0.3), xlim=c(0,1))
fit.CHC<-lm(depth[pick1,"Depth"]~depth[pick1,"Prop_Downstream"])
summary(fit.CHC)
pick1<-depth$Stream=="SC"
jose$prop_downstream<-jose$X..Dist/100
points(jose$X..Dist/100, jose$Mean, col='black', pch=16)
class(jose)
fit<-lm(jose$Mean~(jose$prop_downstream), na.action = na.omit)
summary(fit)
# abline(fit)

fit.lme<-lme(Curvepos~Proportion_Downstream, random = ~1|Stream, data=sat, na.action = na.omit)
fit.lme<-lme(Log_Curvepos~Proportion_Downstream, random = ~1|Stream, data=sat, na.action = na.omit)
summary(fit.lme)
anova(fit.lme)
r.squaredGLMM(fit.lme)
coef<-c(0.012764533,-0.008072929 )
abline(coef, lwd=1)
mtext("c", side=3, adj=0.95,line=-1.5, font=2)
# legend("topleft", 'p = 0.004', lty=1, bty='n')

pick1<-depth$Stream=="CHC"
plot(depth[pick1,"Prop_Downstream"], depth[pick1,"MaxBS"], 
     xlab="Proportion Downstream", ylab="Outer Bank Slope", col="#904596", ylim=c(0,0.7), xlim=c(0,1))
fit.CHC<-lm(depth[pick1,"Depth"]~depth[pick1,"Prop_Downstream"])
summary(fit.CHC)
pick1<-depth$Stream=="SC"
points(depth[pick1,"Prop_Downstream"], depth[pick1,"MaxBS"], 
       xlab="Proportion Downstream", ylab="Outer Bank Slope", col="#4b9645")
fit.SC<-lm(depth[pick1,"Depth"]~depth[pick1,"Prop_Downstream"], na.action = na.omit)
summary(fit.SC)

fit.lme<-lme(MaxBS~Prop_Downstream, random = ~1|Stream, data=depth, na.action = na.omit)
summary(fit.lme)
anova(fit.lme)
coef<-c(0.4646718,-0.2438880)
r.squaredGLMM(fit.lme)
abline(coef, lwd=2, col='grey80')

fit.lme$Intercept
fit.lme

mtext("d", side=1, adj=0.95,line=-1.5, font=2)
legend("topleft", c("SC","CHC"), bty="n",
       col=c("#4b9645",  "#904596"), cex=0.9, pch=16)

#########################################

head(sat)
colnames(sat)
tiff("20230622_Figure3.tiff", width = 6, height = 6, units = 'in', res = 300)
par(mfcol=c(2,2), mar=c(4,4,0.75,0.5))
sat$MR_year<-sat[,"MR"]/67
plot(sat[,"Proportion_Downstream"], sat[,"MR"]/67, 
     xlab="Proportion Downstream",  ylab=expression("Migration Rate (m yr"^-1*")"), col=rgb(0.33, 0.33, 0.33, 0.5), xlim=c(0,1))
mtext("a", side=3, adj=0.95,line=-1.5, font=2)
fit.lme<-lme(MR_year~Proportion_Downstream, random = ~1|Stream, data=sat, na.action = na.omit)
summary(fit.lme)
anova(fit.lme)
r.squaredGLMM(fit.lme)
coef<-c(0.0001934917,0.0010301166)
abline(coef, lwd=2)

plot(sat[,"Proportion_Downstream"], sat[,"MRnorm"], 
     xlab="Proportion Downstream", ylab=expression("Migration Rate (ch-w yr"^-1*")"), 
     col=rgb(0.33, 0.33, 0.33, 0.5), xlim=c(0,1))
mtext("b", side=3, adj=0.95,line=-1.5, font=2)
fit.lme<-lme(MRnorm~Proportion_Downstream, random = ~1|Stream, data=sat, na.action = na.omit)
summary(fit.lme)
anova(fit.lme)
r.squaredGLMM(fit.lme)
coef<-c(0.004360273,-0.004297205)
abline(coef, lwd=2)
# 
# t=sat[,"Proportion_Downstream"]
# t1<-seq(0,1,length.out=10000)
# Mt=sat[,"MRnorm"]
# x <- data.frame(t, Mt)
# xNA<- na.exclude(x)
# Mt<-xNA$Mt
# t<-xNA$t
# singleLL = mle2(minuslogl = LLS, start = list(k = 0.5), data = list(y=Mt),
#                 method="L-BFGS-B", lower=c(0),upper=c(100000))
# 
# #add 1 to K (in AIC caculation) for estimating sigma
# attr(singleLL ,"df") = attr(singleLL,"df")+1
# summary(singleLL)
# AIC(singleLL)
# ssq = sum(((1*exp(-single.k[i,1]*xNA$t)) - xNA$Mt)^2)
# #Extract k value
# single.k=coef(singleLL)[1]
# single.k[i,2]<-sum(((1*exp(-single.k[i,1]*xNA$t)) - xNA$Mt)^2)/length(xNA$t[!is.na(xNA$t)])
# single.exp.fit<-exp(-single.k*t1)

# linear<-lm(sat[,"MRnorm"]~sat[,"Proportion_Downstream"])
# AIC(linear)

plot(slumps[,"Dist."], slumps[,"Shape_Area"], 
     xlab="Proportion Downstream", ylab=expression("Bank Slump Area (m"^2*")"), 
     col=rgb(0.33, 0.33, 0.33, 0.5), xlim=c(0,1))
mtext("c", side=3, adj=0.95,line=-1.5, font=2)
fit.lme<-lme(Shape_Area~Dist., random = ~1|Stream, data=slumps, na.action = na.omit)
summary(fit.lme)
anova(fit.lme)
r.squaredGLMM(fit.lme)
coef<-c(-0.90925,81.73552)
abline(coef, lwd=2)

hist(slumps$Dist., breaks=100, border=F, col='grey50', main="",
     ylab="Bank slump frequency", xlab="Proportion Downstream")
box()
mtext("d", side=3, adj=0.95,line=-1.5, font=2)





dev.off()
pick1<-depth$Stream=="CHC"
plot(depth[pick1,"Dist"], depth[pick1,"Depth"], 
     xlab="Distance Downstream", ylab="Depth (m)", col="#904596", ylim=c(1.5,3.5), xlim=c(0,5000))
fit.CHC<-lm(depth[pick1,"Depth"]~depth[pick1,"Dist"])
summary(fit.CHC)
pick1<-depth$Stream=="SC"
points(depth[pick1,"Dist"], depth[pick1,"Depth"], 
       xlab="Distance Downstream", ylab="Depth (m)", col="#4b9645")
fit.SC<-lm(depth[pick1,"Depth"]~depth[pick1,"Dist"], na.action=na.omit)
summary(fit.SC)
mtext("b", side=1, adj=0.95,line=-1.5, font=2)

##
plot(sat[,"Distance.downstream"], sat[,"Curvature"], 
     xlab="Distance Downstream (m)", ylab="Curvature (rad/min)",col="grey50", 
     xlim=c(0,5000), ylim=c(0,0.06))

##
plot(sat[,"Distance.downstream"], sat[,"Curvature"], 
     xlab="Distance Downstream (m)", ylab="Curvature (rad/min)",col="grey50", 
     xlim=c(0,5000), ylim=c(0,0.025))

##
hist(slumps$Dist., breaks=100, border=F, col='grey50', main="",
     ylab="Bank slump frequency", xlab="Proportion Downstream")
box()

head(sat)
head(sat)
install.packages('OneR')
library(OneR)
sat$bins<-bin(sat$Proportion_Downstream, nbins = 100, labels = c(1:100), method = c("length"), na.omit = TRUE)
bins
head(sat)
barplot(as.numeric(bins))
hist(as.numeric(bins))
barplot(tapply(sat$Slump, sat$bins, sum, na.rm=T))

sum(sat$Slump)
head(sat)

barplot(tapply(sat$Slump, sat$bins, sum, na.rm=T)/sum(sat$Slump))

head(slumps)
hist(slumps$Dist., breaks=100, border=F, col='grey50', main="",
     ylab="Bank slump frequency", xlab="Proportion Downstream")
box()
?hist

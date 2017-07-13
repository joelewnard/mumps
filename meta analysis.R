setwd('Google drive/mumps/data')

out = fitdist(data=c(0.01994302,0.054273504),dist='beta',method='qme',probs=c(0.025,0.975))
mat = cbind(c(0.1994302,0.74074074,0.62678063,0.82621083,2.93447293,3.13390313,2.25071225,2.70655271,7.20797721,8.4045584,0.68376068),
            c(5.04273504,5.75498575,4.98575499,6.01139601,9.54415954,10.02849,9.059820906,9.37321937,16.5811966,20.02849,17.6638177))
mat = mat/100

library('fitdistrplus')
pars = matrix(NA,11,2)
for (i in 1:11){
  pars[i,] = fitdist(data=mat[i,],dist='beta',method='qme',probs=c(0.025,0.975))$estimate
}

pars
sum(pars)
sum(pars[,1])/sum(pars)

file = read.csv('studies.csv')


par(mfrow=c(1,1))
plot(file$lnRR~log(file$dt))

#install.packages('metafor')
library(metafor)
?rma.uni

logRR = file$lnRR
logV = file$var
study = file$study
logdt = log(file$dt)
doses = file$doses
#doses = c(1,2,NA,1,1,1,1,2)

a = study%in%c('a1','a2'); a1 = study=='a1'; a2 = study=='a2'
b = study=='b'
c = study%in%c('c1','c2'); c1 = study=='c1'; c2 = study=='c2'
d = study=='d'
e = study=='e'
f = study=='f'
doses[study=='b'] = NA

#install.packages('gee')
#library(gee)
#mod = glm(logRR~doses+study,weights=(1/logV)); summary(mod)

#install.packages('lme4')
library(lme4)
summary(mod)

mod = lmer(logRR~logdt+(1|study),weights=(1/logV))
BIC(mod)
##### favor FE by AIC
summary(glm(logRR~doses+study,weights=(1/logV)))

mod1 = glm(logRR~logdt+study,weights=(1/logV)); mod0 = glm(logRR~study,weights=(1/logV)); #,subset=(is.na(doses)==F) for model comparison with dose model
BIC(mod0)
names(mod1)
1 - (mod1$deviance/mod1$null.deviance)
1 - (sum((mod1$residuals^2)*mod1$weights))/(sum((mod0$residuals^2)*mod0$weights))
mod = mod1
pars = mvrnorm(1e5,mod$coefficients,vcov(mod))

modDoses = glm(logRR~logdt*doses+study,weights=(1/logV))

BIC(modDoses)
parsDoses = mvrnorm(1e5,modDoses$coefficients[is.na(modDoses$coefficients)==F],vcov(modDoses))

quantile(exp(parsDoses[,9]),c(0.5,0.025,0.975))

rr2dose10y = exp(parsDoses[,1] + parsDoses[,2]*log(10) + parsDoses[,3]*2 + parsDoses[,9]*2*log(10))/exp(parsDoses[,1] + parsDoses[,2]*log(1) + parsDoses[,3]*2 + parsDoses[,9]*2*log(1))
rr1dose10y = exp(parsDoses[,1] + parsDoses[,2]*log(10) + parsDoses[,3]*1 + parsDoses[,9]*1*log(10))/exp(parsDoses[,1] + parsDoses[,2]*log(1) + parsDoses[,3]*1 + parsDoses[,9]*1*log(1))

quantile((1-rr1dose10y)/(1-rr2dose10y),c(0.5,0.025,0.975))
mean(rr2dose10y>rr1dose10y)
rr1dose10y         

x = seq(-1,5,0.01)
summary(modDoses)
est = matrix(NA,1e5,length(x))
for (i in 1:1e5){
  est[i,] = pars[i,1] + pars[i,2]*x + pars[i,3]*mean(study=='a2') + pars[i,4]*mean(study=='b') +
    pars[i,5]*mean(study=='c1') + pars[i,6]*mean(study=='c2') + pars[i,7]*mean(study=='d') + pars[i,8]*mean(study=='e') + pars[i,9]*mean(study=='f')
  print(i)
}
est0 = c()
for (i in 1:1e5){
  est0[i] = pars[i,1] + pars[i,2]*log(0.5) + pars[i,3]*mean(study=='a2') + pars[i,4]*mean(study=='b') +
    pars[i,5]*mean(study=='c1') + pars[i,6]*mean(study=='c2') + pars[i,7]*mean(study=='d') + pars[i,8]*mean(study=='e') + pars[i,9]*mean(study=='f')
}

quantile(1-exp(est0),c(0.5,0.025,0.975))

modMean = apply(est,2,mean)
q95fn = function(x){return(quantile(x,c(0.025,0.975)))}
modCI = apply(est,2,q95fn)

opt.fn = function(lambda,dat,x){
  sqdiff = (dat[1:501] - exp(-lambda*exp(x[1:501])))^2
  return(sum(sqdiff))
}
lambda = c()
for (i in 1:2e3){
  dat = (1-exp(est[i,]))/(1-exp(est0[i]))
  lambda[i] = optim(0.5,fn=opt.fn,dat=dat,x=x)$par
  print(i)
}
quantile(1/lambda,c(0.5,0.025,0.975))

veMat1 = veMat0 = matrix(NA,2e3,length(x))
for (i in 1:2e3){
  veMat1[i,] = (1-exp(est0[i]))*exp(-lambda[i]*exp(x))
  veMat0[i,] = (1-exp(est0[i]))*exp(-lambda[i]*c(exp(x[1:371]),rep(exp(x[371]),230)))
}
veMat1 = apply(veMat1,2,q95fn)
veMat0 = apply(veMat0,2,q95fn)

nc = read.csv('nocomp.csv',header=T)
under5 =rma.uni(yi=nc$logRR[c(1,4,7)],vi=nc$var[c(1,4,7)],method='FE')
over10 = rma.uni(yi=nc$logRR[c(3,6,9)],vi=nc$var[c(3,6,9)],method='FE')

1-exp(quantile(est[,which.min((exp(x)-0)^2)],c(0.5,0.025,0.975)))
1-exp(quantile(est[,which.min((exp(x)-5)^2)],c(0.5,0.025,0.975)))
1-exp(quantile(est[,which.min((exp(x)-10)^2)],c(0.5,0.025,0.975)))
1-exp(quantile(est[,which.min((exp(x)-15)^2)],c(0.5,0.025,0.975)))
1-exp(quantile(est[,which.min((exp(x)-20)^2)],c(0.5,0.025,0.975)))






#layout(matrix(1:4,nrow=2,ncol=2,byrow=T),heights=c(1,0.5),width=c(1,0.75))
par(mgp=c(3,0.35,0))
layout(matrix(c(1,2,1,3,1,4,5,5),nrow=4,ncol=2,byrow=T),heights=c(1.2,0.3,0.9,0.3),widths=c(1,0.5))
par(mar=c(0,2.5,1.5,0.5))
set.seed(2)
par(lwd=0.5)
plot(1,type='n',ylim=c(-5,0.5),xlim=c(-1,3),axes=F,ann=F)
polygon(y=c(modCI[1,],rev(modCI[2,])),x=c(x,rev(x)),lty=0,col='light grey')
lines(y=modMean,x=x)
box(bty='l')
xoffset = runif(length(logdt),-0.075,0.075)
bgs = cols = shapes = c()
bgs[study%in%c('a1','a2')] = 'light blue'; cols[study%in%c('a1','a2')] = 'blue'; shapes[study=='a1'] = 24; shapes[study=='a2'] = 23
bgs[study=='b'] = 'yellow'; cols[study=='b'] = 'dark orange'; shapes[study=='b'] = 22
bgs[study%in%c('c1','c2')] = 'pink'; cols[study%in%c('c1','c2')] = 'red'; shapes[study=='c1'] = 24; shapes[study=='c2'] = 25
bgs[study=='d'] = 'violet'; cols[study=='d'] = 'purple'; shapes[study=='d'] = 24
bgs[study=='e'] = 'green'; cols[study=='e'] = 'dark green'; shapes[study=='e'] = 24
bgs[study=='f'] = 'white'; cols[study=='f'] = 'black'; shapes[study=='f'] = 23
abline(h=0,lty='dashed',col='light grey')
for (i in 1:length(logRR)){
  lines(y=logRR[i]+c(-1.96,1.96)*sqrt(logV[i]),x=rep(logdt[i]+xoffset[i],2),col=cols[i],lwd=0.75)
}

points(y=logRR,x=logdt+xoffset,pch=shapes,bg=bgs,col=cols,cex=0.5 + 0.75*(log(mod1$weights) - min(log(mod1$weights))),lwd=0.75)
axis(2,lwd=0,lwd.ticks=0.5,cex.axis=0.65,at=log(c(0.01,0.1,1,10)),labels=c(0.01,0.1,1,10),las=1,tck=-0.01)
axis(2,lwd=0,lwd.ticks=0.5,cex.axis=0.65,at=log(c(seq(0.02,0.09,0.01),seq(0.2,0.9,0.1),seq(2,9,1))),tck=-0.005,lab=NA)
axis(1,lwd=0,lwd.ticks=0.5,cex.axis=0.65,at=log(c(0.1,0.5,1,2,5,10,15,20)),labels=c(0.1,0.5,1,2,5,10,15,20),tck=-0.01)
mtext('Relative risk of illness given vaccination',cex=0.6,side=2,line=1.5,las=0)
mtext('Time since last dose (y)',cex=0.6,side=1,line=1.5); mtext('A',side=3,adj=0,cex=0.65,font=2)
studies = c('Cohen et al., Emerg Infect Dis 2007\n(1 dose)',
            'Cohen et al., Emerg Infect Dis 2007 \n(2 doses)',
            'Vandermeule et al., Vaccine 2004\n(mixed dosing)',
            'Hilleman et al., NEJM 1967\n(1 dose, school contacts)',
            'Hilleman et al., NEJM 1967\n(1 dose, household contacts)',
            'Tizes et al., MMWR 1973\n(1 dose)',
            'Schlegel et al., BMJ 1999\n(1 dose)',
            'Marin et al., Vaccine 2008\n(2 doses)')
ind = c(1,6,10,21,22,23,24,25)
par(mar=c(0,0,0,0))
plot(1,type='n',ylim=c(0,100),xlim=c(0,100),axes=F,ann=F)
for (i in 1:length(studies)){
  lines(x=c(1,25),col=cols[ind[i]],y=rep(100-10*i,2))
  points(x=mean(c(1,25)),col=cols[ind[i]],bg=bgs[ind[i]],y=100-10*i,pch=shapes[ind[i]],cex=1,lwd=0.75)
  text(adj=0,x=30,y=100-10*i,studies[i],cex=0.65)
}
polygon(x=c(1,25,25,1),y=c(5,5,10,10),lty=0,col='light grey'); lines(x=c(1,25),y=rep(7.5,2)); text(adj=0,x=30,y=7.5,'Pooled est. (95% CI)',cex=0.65)


par(mar=c(0.5,1.75,0.5,0.5))
hist((1-exp(est0))*100,axes=F,ann=F,breaks=200,lty=0,col='light grey',xlim=c(88,100),freq=F)
box(bty='l')
axis(1,lwd.ticks=0.5,lwd=0,cex.axis=0.65,at=seq(88,100,4),tck=-0.04)
axis(2,lwd.ticks=0.5,lwd=0,cex.axis=0.65,las=1,tck=-0.04)
mtext('B',side=3,adj=0,cex=0.65,font=2)
mtext('Effectiveness (%) at 6 months',cex=0.6,side=1,line=1)
mtext('Density',cex=0.6,line=1.25,side=2,las=0)

par(mar=c(0,1.75,2.5,0.5))
#mtext('Time-varying',cex=0.4,side=3,line=1,adj=0);
plot(veMat1[1,],ylim=log(c(0.1,1)),xlim=exp(c(0,4.38)),type='n',axes=F,ann=F)
mtext('C',side=3,adj=0,cex=0.65,font=2)
polygon(y=log(c(veMat1[1,],rev(veMat1[2,]))),x=c(exp(x),rev(exp(x))),lty=0,col='light grey')#rgb(0,0,1,alpha=0.25))
lines(y=log(apply(veMat1,2,mean)),x=exp(x))#col='dark blue')
#polygon(y=c(veMat0[1,],rev(veMat0[2,])),x=c(exp(x),rev(exp(x))),lty=0,col=rgb(1,0,0,alpha=0.25))
#lines(y=apply(veMat0,2,mean),x=exp(x),col='dark red')
#text(x=40,y=0.9,adj=0,cex=0.6,'Continuous waning',col='dark blue',xpd=T)
#text(x=40,y=0.8,adj=0,cex=0.6,'Limited waning',col='dark red',xpd=T)
points(y=log(1-exp(logRR)),x=exp(logdt),pch=shapes,bg=bgs,col=cols,cex=0.25+0.25*(log(mod1$weights) - min(log(mod1$weights))),lwd=0.5)
box(bty='l')
axis(1,lwd=0,lwd.ticks=0.5,cex.axis=0.65,tck=-0.04)
axis(2,lwd=0,lwd.ticks=0.5,cex.axis=0.65,las=1,at=log(seq(0,1,0.1)),labels=c(NA,10,rep(NA,8),100),tck=-0.04)
mtext('Vaccine effectiveness (%)',cex=0.6,side=2,line=1.25,las=0)
mtext('Time since last dose (y)',cex=0.6,side=1,line=1.25); #mtext('C',side=3,adj=0,cex=0.65,font=2)











nc$logRR[1] = log(0.05)
par(mar=c(3,4,1.5,0.5))
cols = rep(c('red','blue','blue'),2); bgs = rep(c('pink','light blue','light blue'),2)
plot(1,type='n',xlim=c(-3,2.3),ylim=c(-1.8,0.2),axes=F,ann=F); abline(v=0,lty='dashed',col='light grey')
ind = c(1,4,7,3,6,9); xs = c(0,0.2,0.4,1,1.2,1.4)
for (i in 1:length(ind)){
  lines(x=nc$logRR[ind[i]]+c(ifelse(i==1,0,-1.96),1.96)*sqrt(nc$var[ind[i]]),y=-rep(xs[i],2),col=cols[i])
}
points(x=nc$logRR[ind],y=-xs,cex=0.5,pch=rep(c(21,24,25),2),col=cols,bg=bgs,lwd=0.75)
lines(x=c(under5$ci.lb,under5$ci.ub),y=-rep(0.6,2),lwd=1)
points(x=under5$b,y=-0.6,pch=22,bg='light grey',lwd=1,cex=0.75)
lines(x=c(over10$ci.lb,over10$ci.ub),y=-rep(1.6,2),lwd=0.75)
points(x=over10$b,y=-1.6,pch=22,bg='light grey',lwd=1,cex=0.75)
text(y=-c(0.4,1.4),adj=0,x=c(-3.75,-3.75),c('0-5\nyears','11+\nyears'),cex=0.65,xpd=T)
axis(1,lwd=0.5,lwd.ticks=0.5,cex.axis=0.65,at=log(c(0.05,0.1,1,10)),labels=c(0,0.1,1,10),tck=-0.04,las=1); axis.break(axis=1,breakpos=log(0.075),style='zigzag',brw=0.04)
axis(1,lwd=0,lwd.ticks=0.5,cex.axis=0.65,at=log(c(seq(0.2,0.9,0.1),seq(2,9,1))),tck=-0.02,lab=NA)
axis(2,at=-c(0,0.8),lwd=0.5,labels=NA,tck=-0.04); axis(2,at=-c(1,1.8),lwd=0.5,labels=NA,tck=-.04)
mtext('Relative risk of illness (ref. 6-10y since last dose)',cex=0.6,side=1,line=1.5);mtext('D',side=3,adj=0,cex=0.65,font=2)
mtext('Time since last dose',line=3,side=2,cex=0.6,las=0)
studies = c('Schaffzin et al., Pediatrics 2007\n(2 doses)','Preeta et al., PIDJ 2014\n(mixed dosing, boys)',
            'Preeta et al., PIDJ 2014\n(mixed dosing, girls)')
par(mar=c(0,0,1.5,0))
plot(1,type='n',ylim=c(0,100),xlim=c(0,100),axes=F,ann=F)
for (i in 1:length(studies)){
  lines(x=c(1,25),col=cols[i],y=rep(100-20*i,2))
  points(x=mean(c(1,25)),col=cols[i],bg=bgs[i],y=100-20*i,pch=c(21,24,25)[i],cex=0.5,lwd=0.75)
  text(adj=0,x=30,y=100-20*i,studies[i],cex=0.65)
}
lines(x=c(1,25),y=rep(20,2),lwd=1); points(x=mean(c(1,25)),y=20,pch=22,bg='light grey',lwd=1,cex=0.75); text(adj=0,x=30,y=20,'Pooled est. (95% CI)',cex=0.65)




pars = tryopt$par
prop = N/sum(N)

beta0 = betaE = matrix(NA,10,10)
for (i in 1:10){
  beta0[i,] = contact[i,]*exp(pars[i])*exp(pars[11])*prop[i]/prop
  betaE[i,] = contact[i,]*exp(pars[i])*exp(pars[11])*(startS[i]/sum(N))/prop
}

r0 = beta0/(1/5); rE = betaE/(1/5)

#### max eigenvalues
max(eigen(r0)$values)
max(eigen(rE)$values)

#### age specific values
plot(apply(r0,2,sum),type='l',ylim=c(0,6))
lines(apply(rE,2,sum),col='red')


##### time varying component
beta0t = betaEt = array(NA,dim=c(50,10,10))
for (t in 1:50) for (i in 1:10){
  beta0t[t,i,] = alpha[t,i]*contact[i,]*exp(pars[i])*exp(pars[11])*prop[i]/prop
  betaEt[t,i,] = alpha[t,i]*contact[i,]*exp(pars[i])*exp(pars[11])*((S[t,i]+P[t,i])/sum(N))/prop
}
r0t = rEt = c()
for (t in 1:50){
  if (is.complex(eigen(beta0t[t,,]*5)$values)==F){
    r0t[t] = max(eigen(beta0t[t,,]*5)$values)
  }
  if (is.complex(eigen(betaEt[t,,]*5)$values)==F){
    rEt[t] = max(eigen(betaEt[t,,]*5)$values) 
  }
}
plot(na.approx(r0t),type='l',ylim=c(0,6),ylab='R0',x=1968:(1968+49))


r0age = rEage = matrix(NA,50,10)
for (t in 1:50){
  r0age[t,] = apply(beta0t[t,,]*5,2,sum)
  rEage[t,] = apply(betaEt[t,,]*5,2,sum)  
}
plot(r0age[,10])




r0age = r0age[2:50,]; rEage = rEage[2:50,]



ages = c('<1y','1-4y','5-9y','10-14y','15-19y','20-24y','25-29y','30-39y','40-64y','65+y')
layout(matrix(c(11,12,1:10),byrow=T,nrow=6,ncol=2),heights=c(1.35,rep(1,4),1.35));  par(mgp=c(3,0.25,0)); par(tck=-0.03)
par(mar=c(0.5,2.5,1,2.5)); par(lwd=0.5)
for (a in 1:10){
  if (a==9){
    par(mar=c(3,2.5,1,2.5)); par(lwd=0.5)
  }
  plot(r0age[,a],x=1968:2016,type='l',col='red',lwd=1,ylim=c(0,7),axes=F,ann=F)
  if (a==1){
    text(expression(R[0]),cex=0.6,col='red',x=1970,y=6,adj=0)
    text(expression(R[E]),cex=0.6,col='blue',x=1970,y=4.5,adj=0)
  }
  abline(h=1,lty='dotted',col='grey')
  lines(rEage[,a],x=1968:2016,lwd=1,col='blue')
  box(bty='l')
  axis(2,las=1,cex.axis=0.5,lwd=0,lwd.ticks=0.5)
  if (a<9){
    axis(1,cex.axis=0.5,lwd=0,lwd.ticks=0.5,labels=NA)
  } else{
    axis(1,cex.axis=0.5,lwd=0,lwd.ticks=0.5,at=seq(1968,2018,4))
    mtext('Year',cex=0.45,line=1.25,side=1)
  }
  mtext(ages[a],side=3,adj=0,cex=0.5,font=2,line=0)
  mtext(expression(paste(italic(R)[0],' (infector age)',sep='')),cex=0.45,side=2,line=1)
}
par(mar=c(2,2.5,1.5,2.5));
plot(na.approx(r0t[2:50]),x=1968:2016,type='l',col='red',lwd=1.5,ylim=c(0,6),axes=F,ann=F)
abline(h=1,lty='dotted',col='grey')
mtext(expression(italic(R)[0]),cex=0.45,side=2,line=1)
mtext('Year',cex=0.45,line=1.25,side=1)
axis(2,las=1,cex.axis=0.5,lwd=0,lwd.ticks=0.5)
axis(1,cex.axis=0.5,lwd=0,lwd.ticks=0.5,at=seq(1968,2018,4))
#lines(na.approx(rEt[2:50]),x=1968:2016,lwd=1,col='blue')
box(bty='l')
mtext(expression(paste('Decline in ',italic(R)[0],' over time',sep='')),side=3,adj=0,cex=0.5,font=2,line=0)

plot(y=apply(r0,2,sum),type='l',ylim=c(0,6),x=c(0.5,3,seq(7.5,27.5,5),35,52.5,72.5),axes=F,ann=F,col='red',lwd=1)
text(expression(R[0]),cex=0.6,col='red',x=40,y=5.5,adj=0)
text(expression(R[E]),cex=0.6,col='blue',x=40,y=4.5,adj=0)
lines(apply(rE,2,sum),col='blue',x=c(0.5,3,seq(7.5,27.5,5),35,52.5,72.5),lwd=1)
abline(h=1,lty='dotted',col='grey')
mtext(expression(italic(R)[0]),cex=0.45,side=2,line=1)
mtext('Age',cex=0.45,line=1.25,side=1)
axis(2,las=1,cex.axis=0.5,lwd=0,lwd.ticks=0.5)
axis(1,cex.axis=0.5,lwd=0,lwd.ticks=0.5)
box(bty='l')
mtext(expression(paste(italic(R)[0],' by infector age, 1967 equilibrium',sep='')),side=3,adj=0,cex=0.5,font=2,line=0)

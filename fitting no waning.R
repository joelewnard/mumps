load('d1ts.Rdata'); d1ts = c(0,d1ts)
load('d2ts.Rdata'); d2ts = c(0,d2ts)

v1 = rbind(d1ts,matrix(0,9,length(d1ts)))
v2 = rbind(rep(0,length(d2ts)),d2ts,matrix(0,8,length(d2ts)))

srVax.fn = function(t,y,parms,v1,v2){
  alpha = 1/(365*c(1,4,5,5,5,5,5,10,25,15))
  
  rho = parms[1]
  w = 0#parms[2]/365
  p = 0.7
  pv = parms[3]
  nu = 1-parms[4]
  
  lambda = parms[5:14]
  
  S = y[1:10]
  R = y[11:20]
  V = y[21:30]
  P = y[31:40]
  N = sum(y[1:40])
  
  lambda = lambda*(S+R+V+P)/(p*S+pv*P)
  
  mu = c(1/80,rep(0,9))/365
  
  dS = mu*N - (lambda + alpha)*S
  dS[2:10] = dS[2:10] + alpha[1:9]*S[1:9]*(1-v1[1:9])
  
  dR = lambda*(S + P) - alpha*R
  dR[2:10] = dR[2:10] + alpha[1:9]*R[1:9]
  
  dV = -(alpha + w)*V
  dV[2:10] = dV[2:10] + alpha[1:9]*(V[1:9] + nu*v1[1:9]*S[1:9] + rho*nu*v2[1:9]*P[1:9])
  
  dP = -(alpha + lambda)*P + w*V
  dP[2:10] = dP[2:10] + alpha[1:9]*(P[1:9]*(1 - rho*nu*v2[1:9]) + (1-nu)*v1[1:9]*S[1:9])
  
  return(list(c(dS,dR,dV,dP)))  
}

#setwd('~/Google drive/mumps/data')
#save(init2,file='init2.Rdata') ### this is the initial values needed to RUN EVERYTHING
#init2 = init1[[1]][[1]]

curState = array(NA,dim=c(1000,49,40))
for (j in 1:1000){
  curState[j,1,] = init2[1:40]
  for (t in 2:49){
    curState[j,t,] = ode(y=curState[j,t-1,],parms=c(1,omega[j],pv[j],rr[j],(1/repCas)*exp(vals.obs[,t+6])/(1e5*365)),times=(0:1)*365,
                         v1=v1[,t],v2=v2[,t],func=srVax.fn)[2,2:41]
  }
  print(j)
}

multNW = prevNW = array(NA,dim=c(1000,49,10))
for (j in 1:1000){
  for (t in 2:49){
    Lambda = (1/repCas)*exp(vals.obs[,t+6])/(1e5*365)
    
    lambda = Lambda*ns/(0.7*curState[j,t,1:10] + pv[j]*curState[j,t,31:40])
    prevNW[j,t,] = lambda*(curState[j,t,1:10] + curState[j,t,31:40])*5
    for (i in 1:10){
      multNW[j,t,i] = lambda[i]/sum(beta*exp(opt1967$par[i])*contact[i,]*prevNW[j,t,]/ns)
    }
    lambda.save[j,t,] = lambda
  }
}
#curStateNW = curState
setwd('~/Google drive/mumps/code')
outputNoWan = list(lambda.save,init1,curState[,,1:10],curState[,,11:20],curState[,,21:30],curState[,,31:40])
#save(outputNoWan,file='outputNoWan.Rdata')
#save(multNW,file='multNW.Rdata')
#curStateNW = curState
#save(curStateNW,file='curStateNW.Rdata')

#r0sNoWan = r0t.fn(outputNoWan)
#save(r0sNoWan,file='r0sNoWan.Rdata')
###########################################################
#### base situation (with waning)
###########################################################

srVax.fn = function(t,y,parms,v1,v2){
  alpha = 1/(365*c(1,4,5,5,5,5,5,10,25,15))
  
  rho = parms[1]
  w = parms[2]/365
  p = 0.7
  pv = parms[3]
  nu = 1-parms[4]
  
  lambda = parms[5:14]
  
  S = y[1:10]
  R = y[11:20]
  V = y[21:30]
  P = y[31:40]
  N = sum(y[1:40])
  
  lambda = lambda*(S+R+V+P)/(p*S+pv*P)
  
  mu = c(1/80,rep(0,9))/365
  
  dS = mu*N - (lambda + alpha)*S
  dS[2:10] = dS[2:10] + alpha[1:9]*S[1:9]*(1-v1[1:9])
  
  dR = lambda*(S + P) - alpha*R
  dR[2:10] = dR[2:10] + alpha[1:9]*R[1:9]
  
  dV = -(alpha + w)*V
  dV[2:10] = dV[2:10] + alpha[1:9]*(V[1:9] + nu*v1[1:9]*S[1:9] + rho*nu*v2[1:9]*P[1:9])
  
  dP = -(alpha + lambda)*P + w*V
  dP[2:10] = dP[2:10] + alpha[1:9]*(P[1:9]*(1 - rho*nu*v2[1:9]) + (1-nu)*v1[1:9]*S[1:9])
  
  return(list(c(dS,dR,dV,dP)))  
}


curState = array(NA,dim=c(1000,49,40))
for (j in 1:1000){
  curState[j,1,] = init2[1:40]
  for (t in 2:49){
    curState[j,t,] = ode(y=curState[j,t-1,],parms=c(1,omega[j],pv[j],rr[j],(1/repCas)*exp(vals.obs[,t+6])/(1e5*365)),times=(0:1)*365,
                         v1=v1[,t],v2=v2[,t],func=srVax.fn)[2,2:41]
  }
  print(j)
}
#load('outputSolved.Rdata')
#load('alphas.Rdata')
#lambda = outputSolved[[1]]; init1 = outputSolved[[2]]; S = outputSolved[[3]]; R = outputSolved[[4]]; V = outputSolved[[5]]; P = outputSolved[[6]]; #preds = outputSolved[[7]]


mult = lambda.save = curI = array(NA,dim=c(1000,49,10))
for (j in 1:1000){
  for (t in 2:49){
    Lambda = (1/repCas)*exp(vals.obs[,t+6])/(1e5*365)
    
    lambda = Lambda*ns/(0.7*curState[j,t,1:10] + pv[j]*curState[j,t,31:40])
    prev = lambda*(curState[j,t,1:10] + curState[j,t,31:40])*5
    curI[j,t,] = prev
    for (i in 1:10){
      mult[j,t,i] = lambda[i]/sum(beta*exp(opt1967$par[i])*contact[i,]*prev/ns)
    }
    lambda.save[j,t,] = lambda
  }
}

#save(curI,file='curI.Rdata')
#alphas = mult; save(alphas,file='alphas.Rdata')
#outputSolved = list(lambda.save,init1,curState[,,1:10],curState[,,11:20],curState[,,21:30],curState[,,31:40])
#save(outputSolved,file='outputSolved.Rdata')


r0t.fn = function(input){
#init = outputSolved
  lambda = input[[1]]; init1 = input[[2]]; S = input[[3]]; R = input[[4]]; V = input[[5]]; P = input[[6]]; #preds = input[[7]]
  
  for (i in 1:10) for (j in 1:10){
    contact[i,j] = contact[j,i] = mean(c(contact[i,j],contact[j,i]))
  }
  
  I = array(NA,dim=c(1000,50,10))
  for (j in 1:1000) for (i in 1:10) for (t in 2:49){
    I[j,t,i] = lambda[j,t,i]*(S[j,t,i] + P[j,t,i])*5
  }
  
  pars = opt1967$par; b = exp(pars[1:10]); beta = exp(pars[11])
  alphas = array(NA,dim=c(1000,50,10))
  for (j in 1:1000) for (i in 1:10) for (t in 2:49){
    alphas[j,t,i] = lambda[j,t,i]/(b[i]*beta*sum(contact[i,]*I[j,t,]/N))
  }
  
  beta0t = array(NA,dim=c(1000,50,10,10))
  for (j in 1:1000) for (t in 1:49) for (i in 1:10){
    beta0t[j,t,i,] = alphas[j,t,i]*contact[i,]*exp(pars[i])*exp(pars[11])*prop[i]/prop
  }

  r0t = matrix(NA,1000,49)
  for (j in 1:1000) for (t in 2:49){
    r0t[j,t] = max(eigen(beta0t[j,t,,]*5)$values,na.rm=T)
  }
  #setwd('~/Google drive/mumps/papers/figures')
  #save(r0t,file='r0t.Rdata')
  t1 = 1:49; t2 = t1^2; t3 = t1^3; t4 = t1^4; t5 = t1^5; t6 = t1^6; t7 = t1^7; t8 = t1^8; t9 = t1^9; t10 = t1^10; t11 = t1^11; t12 = t1^12; t13 = t1^13; t14 = t1^14; t15 = t1^15
  r0s = matrix(NA,1000,49)
  for (j in 1:1000){
    r0modA = glm(log(r0t[j,])~t1);
    r0modB = glm(log(r0t[j,])~t1+t2);
    if (BIC(r0modB)<BIC(r0modA)){
      r0modA = r0modB
      r0modB = glm(log(r0t[j,])~t1+t2+t3)
      if (BIC(r0modB)<BIC(r0modA)){
        r0modA = r0modB
        r0modB = glm(log(r0t[j,])~t1+t2+t3+t4)
        if (BIC(r0modB)<BIC(r0modA)){
          r0modA = r0modB
          r0modB = glm(log(r0t[j,])~t1+t2+t3+t4+t5)
          if (BIC(r0modB)<BIC(r0modA)){
            r0modA = r0modB
            r0modB = glm(log(r0t[j,])~t1+t2+t3+t4+t5+t6)
            if (BIC(r0modB)<BIC(r0modA)){
              r0modA = r0modB
              r0modB = glm(log(r0t[j,])~t1+t2+t3+t4+t5+t6+t7)
              if (BIC(r0modB)<BIC(r0modA)){
                r0modA = r0modB
                r0modB = glm(log(r0t[j,])~t1+t2+t3+t4+t5+t6+t7+t8)
                if (BIC(r0modB)<BIC(r0modA)){
                  r0modA = r0modB
                  r0modB = glm(log(r0t[j,])~t1+t2+t3+t4+t5+t6+t7+t8+t9)
                  if (BIC(r0modB)<BIC(r0modA)){
                    r0modA = r0modB
                    r0modB = glm(log(r0t[j,])~t1+t2+t3+t4+t5+t6+t7+t8+t9+t10)
                  }
                }
              }
            }
          }
        }
      }
    }
    x = mvrnorm(1,coef(r0modB),vcov(r0modB))
    for (t in 2:49){
      r0s[j,t] = exp(x%*%(t^(0:(length(coef(r0modB))-1))))
    }
  }
  q95fn = function(x){return(quantile(x,c(0.5,0.025,0.975),na.rm=T))}
  r0s = apply(r0s,2,q95fn)
  return(r0s)
}

r0s = r0t.fn(outputSolved)
plot(r0s[2,],ylim=c(0,6),type='l')
lines(r0s[3,])
lines(r0s[1,],lwd=1)
#save(r0s,file='r0s.Rdata')

###########################################################
###########################################################
#### waning with 1% decline in reporting per year
###########################################################
###########################################################

init2 = c(startS,startR,rep(0,20))
curState = array(NA,dim=c(1000,49,40))
for (j in 1:1000){
  curState[j,1,] = init2
  for (t in 2:49){
    curState[j,t,] = ode(y=curState[j,t-1,],parms=c(1,omega[j],pv[j],rr[j],(1/(repCas*0.99^(t-1)))*exp(vals.obs[,t+6])/(1e5*365)),times=(0:1)*365,
                         v1=v1[,t],v2=v2[,t],func=srVax.fn)[2,2:41]
  }
  print(j)
}

mult = lambda.save = array(NA,dim=c(1000,49,10))
for (j in 1:1000){
  for (t in 2:49){
    Lambda = (1/repCas)*exp(vals.obs[,t+6])/(1e5*365)
    
    lambda = Lambda*ns/(0.7*curState[j,t,1:10] + pv[j]*curState[j,t,31:40])
    prev = lambda*(curState[j,t,1:10] + curState[j,t,31:40])*5
    for (i in 1:10){
      mult[j,t,i] = lambda[i]/sum(beta*exp(opt1967$par[i])*contact[i,]*prev/ns)
    }
    lambda.save[j,t,] = lambda
  }
}

outputSolvedRep1 = list(lambda.save,init1,curState[,,1:10],curState[,,11:20],curState[,,21:30],curState[,,31:40])
#save(outputSolvedRep1,file='outputSolvedRep1.Rdata')

r0sRep1 = r0t.fn(outputSolvedRep1)
plot(r0s[2,],ylim=c(0,6),type='l')
lines(r0s[3,])
lines(r0s[1,],lwd=1)

lines(r0sRep1[3,],col='blue')
lines(r0sRep1[2,],col='blue')


###########################################################
###########################################################
#### waning with 1% decline in reporting per year
###########################################################
###########################################################

init2 = c(startS,startR,rep(0,20))
curState = array(NA,dim=c(1000,49,40))
for (j in 1:1000){
  curState[j,1,] = init2
  for (t in 2:49){
    curState[j,t,] = ode(y=curState[j,t-1,],parms=c(1,omega[j],pv[j],rr[j],(1/(repCas*0.98^(t-1)))*exp(vals.obs[,t+6])/(1e5*365)),times=(0:1)*365,
                         v1=v1[,t],v2=v2[,t],func=srVax.fn)[2,2:41]
  }
  print(j)
}

mult = lambda.save = array(NA,dim=c(1000,49,10))
for (j in 1:1000){
  for (t in 2:49){
    Lambda = (1/repCas)*exp(vals.obs[,t+6])/(1e5*365)
    
    lambda = Lambda*ns/(0.7*curState[j,t,1:10] + pv[j]*curState[j,t,31:40])
    prev = lambda*(curState[j,t,1:10] + curState[j,t,31:40])*5
    for (i in 1:10){
      mult[j,t,i] = lambda[i]/sum(beta*exp(opt1967$par[i])*contact[i,]*prev/ns)
    }
    lambda.save[j,t,] = lambda
  }
}

outputSolvedRep2 = list(lambda.save,init1,curState[,,1:10],curState[,,11:20],curState[,,21:30],curState[,,31:40])
save(outputSolvedRep2,file='outputSolvedRep2.Rdata')

r0sRep2 = r0t.fn(outputSolvedRep2)
plot(r0s[2,],ylim=c(0,6),type='l')
lines(r0s[3,])
lines(r0s[1,],lwd=1)

lines(r0sRep2[3,],col='blue')
lines(r0sRep2[2,],col='blue')

######################################################################
######################################################################
##### defining susx for analyses #####################################
######################################################################
######################################################################

setwd('~/Google drive/mumps/code')
susx = susx1 = susx2 = susxNoWan = array(NA,dim=c(10,48,1000))

for (i in 1:10) for (t in 1:48) for (j in 1:1000){
  
  S = outputSolved[[3]]; R = outputSolved[[4]]; V = outputSolved[[5]]; P = outputSolved[[6]]; #preds = outputSolved[[7]]
  susx[i,t,j] = (S[j,t+1,i] + P[j,t+1,i])/N[i]
  
  S = outputSolvedRep1[[3]]; R = outputSolvedRep1[[4]]; V = outputSolvedRep1[[5]]; P = outputSolvedRep1[[6]]; #preds = outputSolvedRep1[[7]]
  susx1[i,t,j] = (S[j,t+1,i] + P[j,t+1,i])/N[i]
  
  S = outputSolvedRep2[[3]]; R = outputSolvedRep2[[4]]; V = outputSolvedRep2[[5]]; P = outputSolvedRep2[[6]]; #preds = outputSolvedRep1[[7]]
  susx2[i,t,j] = (S[j,t+1,i] + P[j,t+1,i])/N[i]
  
  S = outputNoWan[[3]]; R = outputNoWan[[4]]; V = outputNoWan[[5]]; P = outputNoWan[[6]]; #preds = outputSolvedRep1[[7]]
  susxNoWan[i,t,j] = (S[j,t+1,i] + P[j,t+1,i])/N[i]
}
#save(susxNoWan,file='susxNoWan.Rdata')
#save(susx,file='susx.Rdata'); save(susx1,file='susx1.Rdata'); save(susx2,file='susx2.Rdata')

#################################################################
#################################################################
#################### plotting R0s ###############################
#################################################################
#################################################################

setwd('~/Google drive/mumps/papers/figures')
#save(r0s,file='r0s.Rdata'); save(r0sRep1,file='r0sRep1.Rdata'); save(r0sRep2,file='r0sRep2.Rdata')
load('r0s.Rdata'); load('r0sRep1.Rdata'); load('r0sRep2.Rdata')



par(lwd=0.5); par(tck=-0.03); par(mgp=c(0,0.25,0)); par(mar=c(3,3,1,1))
plot(r0s[1,],type='n',axes=F,ann=F,x=1968:2016,ylim=c(0,7))
box(bty='l')
axis(1,cex.axis=0.5,lwd=0,lwd.ticks=0.5,at=seq(0,75,15),las=1)
axis(1,cex.axis=0.5,at=seq(1968,2016,4),lwd.ticks=0.5,lwd=0)
axis(2,cex.axis=0.5,at=seq(0,6,1),lwd.ticks=0.5,lwd=0,las=1)
mtext(expression(italic(R)[0](italic(t))),cex=0.65,line=1,xpd=T,side=2,las=1)
mtext('Year',cex=0.65,line=2,xpd=T,side=1)

polygon(y=c(r0s[2,],rev(r0s[3,])),x=c(1968:2016,2016:1968),lty=0,col='light grey')
lines(y=r0s[2,],x=1968:2016); lines(y=r0s[3,],x=1968:2016)

polygon(y=c(r0sRep1[2,],rev(r0sRep1[3,])),x=c(1968:2016,2016:1968),lty=0,col=rgb(0,0.125,1,0.25))
lines(r0sRep1[2,],col='dark blue',lwd=1,x=1968:2016); lines(r0sRep1[3,],col='dark blue',lwd=1,x=1968:2016)

polygon(y=c(r0sRep2[2,],rev(r0sRep2[3,])),x=c(1968:2016,2016:1968),lty=0,col=rgb(1,0.125,0,0.25))
lines(r0sRep2[2,],col='dark red',lwd=1,x=1968:2016); lines(r0sRep2[3,],col='dark red',lwd=1,x=1968:2016)

polygon(y=c(7,7,6.5,6.5),x=c(1984,1990,1990,1984),col='light grey',lty=0)
lines(y=c(7,7),x=c(1984,1990)); lines(y=c(6.5,6.5),x=c(1984,1990))
text(y=6.75,x=1993,adj=0,'No change in reporting',cex=0.5)

polygon(y=c(6,6,5.5,5.5),x=c(1984,1990,1990,1984),col=rgb(0,0.125,1,0.25),lty=0)
lines(y=c(6,6),x=c(1984,1990),lwd=1,col='dark blue'); lines(y=c(5.5,5.5),x=c(1984,1990),lwd=1,col='dark blue')
text(y=5.75,x=1993,adj=0,'1% annual decline in reporting',cex=0.5)

polygon(y=c(5,5,4.5,4.5),x=c(1984,1990,1990,1984),col=rgb(1,0.125,0,0.25),lty=0)
lines(y=c(5,5),x=c(1984,1990),lwd=1,col='dark red'); lines(y=c(4.5,4.5),x=c(1984,1990),lwd=1,col='dark red')
text(y=4.75,x=1993,adj=0,'2% annual decline in reporting',cex=0.5)


########################################################################################
########################################################################################
#################### plotting ests of suscept population ###############################
########################################################################################
########################################################################################



layout(matrix(c(1:10,11,11),nrow=4,ncol=3,byrow=T))#,heights=c(rep(1,4),0.25))
par(mar=c(2.5,3.5,1.5,2.5))
for (i in 1:10){
  
  par(tck=-0.035)
  
  plot(y=c(susx[i,,2],NA),x=1968:2016,type='n',axes=F,ann=F,ylim=c(-0.12,0.02))
  mtext(LETTERS[i],cex=0.5,font=2,side=3,adj=0)
  mtext(paste('Ages ',ages[i],sep=''),font=3,side=3,adj=1,cex=0.5)
  abline(h=0,col='dark grey',lwd=0.5)
  
  polyS = apply((susx1[i,,]-susx[i,,])/susx[i,,],1,q95fn)
  #abline(v=seq(1,49,4),lwd=0.25,col='light grey')
  polygon(y=c(polyS[2,],rev(polyS[3,])),x=c(1968:2015,2015:1968),col=rgb(0,0.125,1,alpha=0.25),lty=0)
  lines(y=polyS[2,],x=1968:2015,col='dark blue',lwd=0.25);   lines(y=polyS[3,],x=1968:2015,col='dark blue',lwd=0.25)
  
  polyS = apply((susx2[i,,]-susx[i,,])/susx[i,,],1,q95fn)
  #abline(v=seq(1,49,4),lwd=0.25,col='light grey')
  polygon(y=c(polyS[2,],rev(polyS[3,])),x=c(1968:2015,2015:1968),col=rgb(1,0.125,0,alpha=0.25),lty=0)
  lines(y=polyS[2,],x=1968:2015,col='dark red',lwd=0.25);   lines(y=polyS[3,],x=1968:2015,col='dark red',lwd=0.25)
  
  
  
  #  polyS = apply(susx1[i,,],1,q95fn)
  #  polygon(y=c(polyS[2,],rev(polyS[3,])),x=c(1968:2015,2015:1968),col=rgb(1,0.125,0,alpha=0.25),lty=0)
  #  lines(y=polyS[2,],x=1968:2015,col='dark red',lwd=0.25);   lines(y=polyS[3,],x=1968:2015,col='dark red',lwd=0.25)
  
  axis(2,at=seq(-0.1,0.1,0.05),labels=seq(-10,10,5),las=1,cex.axis=0.6,lwd.ticks=0.5,lwd=0)
  axis(1,cex.axis=0.6,at=seq(1968,2016,4),lwd.ticks=0.5,lwd=0)
  mtext(side=2,line=1.5,'Change (%) in\nsusceptible population',cex=0.5)
  mtext(side=1,line=1.5,'Year',cex=0.5)
  box(bty='u')
  if (i%in%c(9,10)){
    axis(1,at=seq(1,49,4),cex.axis=0.6,lwd.ticks=0.5,lwd=0,labels=seq(1968,2016,4))  
  } else{
    axis(1,at=seq(1,49,4),cex.axis=0.6,lwd.ticks=0.5,lwd=0,labels=NA)
  }
}

par(mar=c(0,3,0,3))
plot(1,ylim=c(0,1),xlim=c(0,1),axes=F,ann=F,type='n')
polygon(y=c(0.5,0.5,0.6,0.6),x=c(0,0.25,0.25,0),col=rgb(0,0.125,1,alpha=0.25),lty=0)
lines(y=c(0.5,0.5),x=c(0,0.25),col='dark blue')
lines(y=c(0.6,0.6),x=c(0,0.25),col='dark blue')
text(y=0.55,x=0.3,adj=0,'1% reduction per year (40% overall) in reporting probability',xpd=T,cex=0.65)

polygon(y=c(0.3,0.3,0.4,0.4),x=c(0,0.25,0.25,0),col=rgb(1,0.125,0,alpha=0.25),lty=0)
lines(y=c(0.3,0.3),x=c(0,0.25),col='dark red')
lines(y=c(0.4,0.4),x=c(0,0.25),col='dark red')
text(y=0.35,x=0.3,adj=0,'2% reduction per year (63% overall) in reporting probability',xpd=T,cex=0.65)


########################################################################################
########################################################################################
#################### plotting ests of suscept population and R0 without waning #########
########################################################################################
########################################################################################


layout(matrix(c(1:11,11),nrow=4,ncol=3,byrow=T))#,heights=c(rep(1,4),0.25))
par(mar=c(2.5,3.5,1.5,2.5) )
for (i in 1:10){
  
  par(tck=-0.035)
  
  plot(y=c(susx[i,,2],NA),x=1968:2016,type='n',axes=F,ann=F,ylim=c(0,1))
  mtext(LETTERS[i],cex=0.5,font=2,side=3,adj=0)
  mtext(paste('Ages ',ages[i],sep=''),font=3,side=3,adj=1,cex=0.5)
 # abline(h=0,col='dark grey',lwd=0.5)
  
  polyS = apply(susx[i,,],1,q95fn)
  #abline(v=seq(1,49,4),lwd=0.25,col='light grey')
  polygon(y=c(polyS[2,],rev(polyS[3,])),x=c(1968:2015,2015:1968),col=rgb(0.5,0.49,0.51,0.125),lty=0)
  lines(y=polyS[2,],x=1968:2015,lwd=0.25);   lines(y=polyS[3,],x=1968:2015,lwd=0.25)
  
  polyS = apply(susxNoWan[i,,],1,q95fn)
  #abline(v=seq(1,49,4),lwd=0.25,col='light grey')
  polygon(y=c(polyS[2,],rev(polyS[3,])),x=c(1968:2015,2015:1968),col=rgb(34/255,1,34/255,alpha=0.5),lty=0)
  lines(y=polyS[2,],x=1968:2015,col='dark green',lwd=1);   lines(y=polyS[3,],x=1968:2015,col='dark green',lwd=1)
  
  axis(2,at=seq(0,1,0.2),labels=seq(0,100,20),las=1,cex.axis=0.6,lwd.ticks=0.5,lwd=0)
  axis(1,cex.axis=0.6,at=seq(1968,2016,4),lwd.ticks=0.5,lwd=0)
  mtext(side=2,line=1.5,'Susceptible population (%)',cex=0.5)
  mtext(side=1,line=1.5,'Year',cex=0.5)
  box(bty='l')
}
r0sNoWan[,49]

plot(r0s[1,],type='n',axes=F,ann=F,x=1968:2016,ylim=c(0,10))
polygon(y=c(r0s[2,],rev(r0s[3,])),x=c(1968:2016,2016:1968),lty=0,col=rgb(0.5,0.49,0.51,0.125))
lines(y=r0s[2,],x=1968:2016,lwd=0.25); lines(y=r0s[3,],x=1968:2016,lwd=0.25)

polygon(y=c(r0sNoWan[2,],rev(r0sNoWan[3,])),x=c(1968:2016,2016:1968),lty=0,col=rgb(34/255,1,34/255,0.5))
lines(y=r0sNoWan[2,],x=1968:2016,lwd=1,col='dark green'); lines(y=r0sNoWan[3,],x=1968:2016,lwd=1,col='dark green')

lines(y=c(10.5,10.5),x=c(1968,1976),lwd=0.25,xpd=T); lines(y=c(9.5,9.5),x=c(1968,1976),lwd=0.25,xpd=T)
polygon(x=c(1968,1976,1976,1968),y=c(10.5,10.5,9.5,9.5),col=rgb(0.5,0.49,0.51,0.125),lty=0,xpd=T)
text(y=10,'Waning efficacy',cex=0.65,adj=0,x=1979)

lines(y=c(8.5,8.5),x=c(1968,1976),lwd=1,col='dark green'); lines(y=c(7.5,7.5),x=c(1968,1976),lwd=1,col='dark green')
polygon(x=c(1968,1976,1976,1968),y=c(8.5,8.5,7.5,7.5),col=rgb(34/255,1,34/255,0.5),lty=0,xpd=T)
text(y=8,'No waning',cex=0.65,adj=0,x=1979)

box(bty='l')
axis(1,cex.axis=0.5,lwd=0,lwd.ticks=0.5,at=seq(0,75,15),las=1)
axis(1,cex.axis=0.5,at=seq(1968,2016,4),lwd.ticks=0.5,lwd=0)
axis(2,cex.axis=0.5,at=seq(0,10,1),lwd.ticks=0.5,lwd=0,las=1)
mtext(expression(italic(R)[0](italic(t))),cex=0.5,line=1,xpd=T,side=2,las=1)
mtext('Year',cex=0.5,line=1.5,xpd=T,side=1)
mtext(side=3,cex=0.5,'K',adj=0,font=2)


susxNW0 = susxNW10 = susxNW25 = susxNW50 = susxNW75 = susxNW100 = array(NA,dim=c(10,48,1000))
susxNW0a = susxNW10a = susxNW25a = susxNW50a = susxNW75a = susxNW100a = array(NA,dim=c(10,48,1000))
for (i in 1:10) for (t in 1:48) for (j in 1:1000){
  susxNW0[i,t,j] = (curStateNW[j,t+1,i] + curStateNW[j,t+1,i+30] + curStateNW[j,t+1,i+20])/sum(curStateNW[j,t+1,seq(i,40,10)])
  susxNW10[i,t,j] = (curStateNW[j,t+1,i] + curStateNW[j,t+1,i+30] + 0.9*curStateNW[j,t+1,i+20])/sum(curStateNW[j,t+1,seq(i,40,10)])
  susxNW25[i,t,j] = (curStateNW[j,t+1,i] + curStateNW[j,t+1,i+30] + 0.75*curStateNW[j,t+1,i+20])/sum(curStateNW[j,t+1,seq(i,40,10)])
  susxNW50[i,t,j] = (curStateNW[j,t+1,i] + curStateNW[j,t+1,i+30] + 0.50*curStateNW[j,t+1,i+20])/sum(curStateNW[j,t+1,seq(i,40,10)])
  susxNW75[i,t,j] = (curStateNW[j,t+1,i] + curStateNW[j,t+1,i+30] + 0.25*curStateNW[j,t+1,i+20])/sum(curStateNW[j,t+1,seq(i,40,10)])
  susxNW100[i,t,j] = (curStateNW[j,t+1,i] + curStateNW[j,t+1,i+30] + 0*curStateNW[j,t+1,i+20])/sum(curStateNW[j,t+1,seq(i,40,10)])
  
  susxNW0a[i,t,j] = (curStateNW[j,t+1,i] + curStateNW[j,t+1,i+30] + curStateNW[j,t+1,i+20] + curStateNW[j,t+1,i+10])/sum(curStateNW[j,t+1,seq(i,40,10)])
  susxNW10a[i,t,j] = (curStateNW[j,t+1,i] + curStateNW[j,t+1,i+30] + 0.9*(curStateNW[j,t+1,i+20] + curState[j,t+1,i+10]))/sum(curStateNW[j,t+1,seq(i,40,10)])
  susxNW25a[i,t,j] = (curStateNW[j,t+1,i] + curStateNW[j,t+1,i+30] + 0.75*(curStateNW[j,t+1,i+20] + curState[j,t+1,i+10]))/sum(curStateNW[j,t+1,seq(i,40,10)])
  susxNW50a[i,t,j] = (curStateNW[j,t+1,i] + curStateNW[j,t+1,i+30] + 0.50*(curStateNW[j,t+1,i+20] + curState[j,t+1,i+10]))/sum(curStateNW[j,t+1,seq(i,40,10)])
  susxNW75a[i,t,j] = (curStateNW[j,t+1,i] + curStateNW[j,t+1,i+30] + 0.25*(curStateNW[j,t+1,i+20] + curState[j,t+1,i+10]))/sum(curStateNW[j,t+1,seq(i,40,10)])
  susxNW100a[i,t,j] = (curStateNW[j,t+1,i] + curStateNW[j,t+1,i+30] + 0*curStateNW[j,t+1,i+20])/sum(curStateNW[j,t+1,seq(i,40,10)])
}

xs = c(2.5,7.5,12.5,17.5,22.5,27.5,35,52.5,72.5)


layout(matrix(c(1,4,7,
                2,5,8,
                3,6,9),byrow=T,nrow=3),heights=c(1,1,1))
par(lwd=0.5)
par(mar=c(4,2.5,2,1))
ynum = c(27,48,55)
yearlab = c(1986,2006,2013)
par(mgp=c(0,0.35,0))
par(tck=-0.025)
for (k in 1:3){

  plot(x=xs,y=rep(0,9),ylim=c(0,max(exp(vals.obs[,ynum[k]]))),type='n',axes=F,ann=F)
  #plot(rep(0,9),ylim=c(-15,max(vals.obs[,ynum[k]])),type='n',axes=F,ann=F)
  lines(x=xs,y=exp(vals.obs[2:10,ynum[k]])) ## 31  
  box(bty='l')
  axis(side=1,at=seq(0,80,20),cex.axis=0.6,lwd.ticks=0.5,lwd=0)
  axis(side=2,cex.axis=0.6,lwd.ticks=0.5,lwd=0,las=1)
  mtext(side=3,adj=0.5,yearlab[k],font=2,cex=0.75,line=1)
  mtext(expression(paste('Notifications per ',10^5,sep='')),cex=0.5,line=1.5,xpd=T,side=2)
  mtext(side=1,'Age (y)',cex=0.5,line=1.5)
  
  plot(x=xs,y=xs,type='n',axes=F,ann=F,ylim=c(0,1))
  polygon(x=c(xs,rev(xs)),y=c(apply(susx[2:10,ynum[k]-7,],1,q95fn)[2,],
                           rev(apply(susx[2:10,ynum[k]-7,],1,q95fn)[3,])),lty=0,col='light grey')
  lines(x=c(xs),y=apply(susx[2:10,ynum[k]-7,],1,q95fn)[2,],col='black',lwd=0.5)
  lines(x=c(xs),y=apply(susx[2:10,ynum[k]-7,],1,q95fn)[3,],col='black',lwd=0.5)
  box(bty='l')
  axis(side=1,at=seq(0,80,20),cex.axis=0.6,lwd.ticks=0.5,lwd=0)
  axis(side=2,cex.axis=0.6,lwd.ticks=0.5,lwd=0,las=1,at=seq(0,1,0.25),labels=seq(0,100,25))
  mtext('Proportion susceptible (%)',cex=0.5,line=1.5,xpd=T,side=2)
  mtext(side=1,'Age (y)',cex=0.5,line=1.5)
  if (k==1){
    mtext(side=3,'Vaccine waning scenario',font=4,line=1,cex=0.75,xpd=T,adj=0)
  }

  plot(x=xs,y=xs,type='n',axes=F,ann=F,ylim=c(0,1))
  polygon(x=c(xs,rev(xs)),y=c(apply(susxNW0[2:10,ynum[k]-7,],1,q95fn)[2,],
                           rev(apply(susxNW0[2:10,ynum[k]-7,],1,q95fn)[3,])),lty=0,col=rgb(1,0,0,0.15))
  lines(x=c(xs),y=apply(susxNW0[2:10,ynum[k]-7,],1,q95fn)[2,],col=rgb(1,0,0.5,1),lwd=0.5)
  lines(x=c(xs),y=apply(susxNW0[2:10,ynum[k]-7,],1,q95fn)[3,],col=rgb(1,0,0.5,1),lwd=0.5)
#  lines(x=c(xs),y=apply(susxNW0a[2:10,ynum[k]-7,],1,q95fn)[2,],col=rgb(1,0,0.5,1),lwd=0.5,lty='dotted')
#  lines(x=c(xs),y=apply(susxNW0a[2:10,ynum[k]-7,],1,q95fn)[3,],col=rgb(1,0,0.5,1),lwd=0.5,lty='dotted')
  
  polygon(x=c(xs,rev(xs)),y=c(apply(susxNW10[2:10,ynum[k]-7,],1,q95fn)[2,],
                           rev(apply(susxNW10[2:10,ynum[k]-7,],1,q95fn)[3,])),lty=0,col=rgb(0.8,0,0.2,0.15))
  lines(x=c(xs),y=apply(susxNW10[2:10,ynum[k]-7,],1,q95fn)[2,],col=rgb(0.8,0,0.6,1),lwd=0.5)
  lines(x=c(xs),y=apply(susxNW10[2:10,ynum[k]-7,],1,q95fn)[3,],col=rgb(0.8,0,0.6,1),lwd=0.5)
 # lines(x=c(xs),y=apply(susxNW10a[2:10,ynum[k]-7,],1,q95fn)[2,],col=rgb(0.8,0,0.6,1),lwd=0.5,lty='dotted')
 # lines(x=c(xs),y=apply(susxNW10a[2:10,ynum[k]-7,],1,q95fn)[3,],col=rgb(0.8,0,0.6,1),lwd=0.5,lty='dotted')
  
  polygon(x=c(xs,rev(xs)),y=c(apply(susxNW25[2:10,ynum[k]-7,],1,q95fn)[2,],
                           rev(apply(susxNW25[2:10,ynum[k]-7,],1,q95fn)[3,])),lty=0,col=rgb(0.6,0,0.4,0.15))
  lines(x=c(xs),y=apply(susxNW25[2:10,ynum[k]-7,],1,q95fn)[2,],col=rgb(0.6,0,0.7,1),lwd=0.5)
  lines(x=c(xs),y=apply(susxNW25[2:10,ynum[k]-7,],1,q95fn)[3,],col=rgb(0.6,0,0.7,1),lwd=0.5)
  #lines(x=c(xs),y=apply(susxNW25a[2:10,ynum[k]-7,],1,q95fn)[2,],col=rgb(0.6,0,0.7,1),lwd=0.5,lty='dotted')
  #lines(x=c(xs),y=apply(susxNW25a[2:10,ynum[k]-7,],1,q95fn)[3,],col=rgb(0.6,0,0.7,1),lwd=0.5,lty='dotted')
  
  polygon(x=c(xs,rev(xs)),y=c(apply(susxNW50[2:10,ynum[k]-7,],1,q95fn)[2,],
                           rev(apply(susxNW50[2:10,ynum[k]-7,],1,q95fn)[3,])),lty=0,col=rgb(0.4,0,0.6,0.15))
  lines(x=c(xs),y=apply(susxNW50[2:10,ynum[k]-7,],1,q95fn)[2,],col=rgb(0.4,0,0.8,1),lwd=0.5)
  lines(x=c(xs),y=apply(susxNW50[2:10,ynum[k]-7,],1,q95fn)[3,],col=rgb(0.4,0,0.8,1),lwd=0.5)
 # lines(x=c(xs),y=apply(susxNW50a[2:10,ynum[k]-7,],1,q95fn)[2,],col=rgb(0.4,0,0.8,1),lwd=0.5,lty='dotted')
 # lines(x=c(xs),y=apply(susxNW50a[2:10,ynum[k]-7,],1,q95fn)[3,],col=rgb(0.4,0,0.8,1),lwd=0.5,lty='dotted')
  
  polygon(x=c(xs,rev(xs)),y=c(apply(susxNW75[2:10,ynum[k]-7,],1,q95fn)[2,],
                           rev(apply(susxNW75[2:10,ynum[k]-7,],1,q95fn)[3,])),lty=0,col=rgb(0.2,0,0.8,0.15))
  lines(x=c(xs),y=apply(susxNW75[2:10,ynum[k]-7,],1,q95fn)[2,],col=rgb(0.2,0,0.9,1),lwd=0.5)
  lines(x=c(xs),y=apply(susxNW75[2:10,ynum[k]-7,],1,q95fn)[3,],col=rgb(0.2,0,0.9,1),lwd=0.5)
  #lines(x=c(xs),y=apply(susxNW75a[2:10,ynum[k]-7,],1,q95fn)[2,],col=rgb(0.3,0,0.9,1),lwd=0.5,lty='dotted')
  #lines(x=c(xs),y=apply(susxNW75a[2:10,ynum[k]-7,],1,q95fn)[3,],col=rgb(0.3,0,0.9,1),lwd=0.5,lty='dotted')
  
  box(bty='l')
  axis(side=1,at=seq(0,80,20),cex.axis=0.6,lwd.ticks=0.5,lwd=0)
  axis(side=2,cex.axis=0.6,lwd.ticks=0.5,lwd=0,las=1,at=seq(0,1,0.25),labels=seq(0,100,25))
  mtext('Proportion susceptible (%)',cex=0.5,line=1.5,xpd=T,side=2)
  mtext(side=1,'Age (y)',cex=0.5,line=1.5)
  if (k==1){
    mtext(side=3,'Vaccine escape scenario',font=4,line=1,cex=0.75,xpd=T,adj=0)    
  }
  
  if (k==1){
    polygon(y=c(1,1,0.98,0.98),x=c(30,45,45,30),col=rgb(1,0,0,0.15),lty=0)
    lines(y=c(1,1),x=c(30,45),col=rgb(1,0,0.5,1),lwd=0.5)
    lines(y=c(0.98,0.98),x=c(30,45),col=rgb(1,0,0.5,1),lwd=0.5)
    text(y=0.99,x=48,adj=0,expression(paste(VE[WT]==0,'%',sep='')),cex=0.75)
    
    polygon(y=c(0.85,0.85,0.83,0.83),x=c(30,45,45,30),col=rgb(0.8,0,0.2,0.15),lty=0)
    lines(y=c(0.85,0.85),x=c(30,45),col=rgb(0.8,0,0.6,1),lwd=0.5)
    lines(y=c(0.83,0.83),x=c(30,45),col=rgb(0.8,0,0.6,1),lwd=0.5)
    text(y=0.84,x=48,adj=0,expression(paste(VE[WT]==10,'%',sep='')),cex=0.75)
    
    polygon(y=c(0.7,0.7,0.68,0.68),x=c(30,45,45,30),col=rgb(0.6,0,0.4,0.15),lty=0)
    lines(y=c(0.7,0.7),x=c(30,45),col=rgb(0.6,0,0.7,1),lwd=0.5)
    lines(y=c(0.68,0.68),x=c(30,45),col=rgb(0.6,0,0.7,1),lwd=0.5)
    text(y=0.69,x=48,adj=0,expression(paste(VE[WT]==25,'%',sep='')),cex=0.75)
    
    polygon(y=c(0.55,0.55,0.53,0.53),x=c(30,45,45,30),col=rgb(0.4,0,0.6,0.15),lty=0)
    lines(y=c(0.55,0.55),x=c(30,45),col=rgb(0.4,0,0.8,1),lwd=0.5)
    lines(y=c(0.53,0.53),x=c(30,45),col=rgb(0.4,0,0.8,1),lwd=0.5)
    text(y=0.54,x=48,adj=0,expression(paste(VE[WT]==50,'%',sep='')),cex=0.75)
    
    polygon(y=c(0.4,0.4,0.38,0.38),x=c(30,45,45,30),col=rgb(0.2,0,0.8,0.15),lty=0)
    lines(y=c(0.4,0.4),x=c(30,45),col=rgb(0.2,0,0.9,1),lwd=0.5)
    lines(y=c(0.38,0.38),x=c(30,45),col=rgb(0.2,0,0.9,1),lwd=0.5)
    text(y=0.39,x=48,adj=0,expression(paste(VE[WT]==75,'%',sep='')),cex=0.75)
  }
}


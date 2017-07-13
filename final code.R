
setwd('~/Google drive/mumps/data')
template = read.csv('template.csv',header=T)
template = template[1:80,2:40]
years = 1960:2060

inc = read.csv('data.csv',header=T)
inc = inc$cases/inc$pop
inc = c(rep(NA,6),inc) ########## gives nationwide cases/pop from 1960 on consistent with mat timing


vals = read.csv('vals.csv',header=T)[,2]
vals[vals==0] = 0.01

mat = matrix(NA,80,39)
for (i in 1:80) for (j in 1:39){
  if (is.na(template[i,j])==F){
    mat[i,j] = vals[template[i,j]]    
  } else{
    mat[i,j] = NA
  }
}

props = c(58.48533,rep(337.74533,4))
props = props/sum(props)
for (i in c(8:17,22:23)){
  mat[1:5,i] = props*mat[1:5,i]
}


vals = matrix(NA,80,39)

###### FILL IN BLANKS/AVERAGES ACCORDING TO INCIDENCE DATA
for (i in 1:80) for (j in 8:12){
  vals[i,j] = (inc[j]/inc[10])*mat[i,j]
}
for (i in 1:80) for (j in 13:17){
  vals[i,j] = (inc[j]/inc[15])*mat[i,j]
}
for (i in 1:80) for (j in 22:23){
  vals[i,j] = (inc[j]/mean(inc[22:23]))*mat[i,j]
}
for (i in 1:80){
  vals[i,17:22] = na.approx(vals[i,17:22])
}
vals = log(vals)
for (i in 1:80){
  vals[i,23:39] = na.approx(log(mat[i,23:39]))
}

####### FILL IN AGGREGATED AGE GROUPS
old = exp(vals[,35])
meanold = mean(old[16:80])

prop = old[16:80]/meanold
vals.interp = vals
for (t in 6:24){
  vals.interp[16:80,t] = log(prop) + vals[16:80,t]
}
vals.interp[16:20,24] = vals[16,24]

##### do same for 20 and up for 1968:1983 and 1987:1995
meanold = mean(old[31:80])
prop = old[31:80]/meanold
for (t in 25:28){
  vals.interp[31:80,t] = log(prop) + vals[31:80,t]
}
meanold = mean(old[21:80])
prop = old[21:80]/meanold
for (t in 29:34){
  vals.interp[21:80,t] = log(prop) + vals[21:80,t]
}


#save(vals.interp,file='vals.interp.Rdata')
###### Add in years after 1998

vals.in = exp(vals.interp[c(1,2,6,11,16,21,26,31,41,66),8]); plot(vals.in) # values as of 1967
vals.obs = vals.interp[c(1,2,6,11,16,21,26,31,41,66),]

interyears = vals.interp[c(1,2,6,11,16,21,26,31,41,66),39] - log(2)
interyears = matrix(rep(interyears,8),nrow=10,ncol=8,byrow=F)
#exp(interyears)/exp(pars[12])

y2006 = c(0.7977208,        ##### enter data by age; <1
          rep(8.13390313,4),  # 1-4
          rep(9.7008547,5),   # 5-9
          rep(12.4786325,4),  # 10-13
          rep(17.3931624,4),  # 14-17
          rep(30.9259259,6),  # 18-24
          rep(11.8376068,5),  # 25-29
          rep(11.1965812,10), # 30-39
          rep(8.49002849,10), # 40-49
          rep(5.71225071,10), # 50-59
          rep(2.50712251,21)) # 60+ 
y2006 = y2006*6500/4017 ### scaling to final size (reports included 4017 of 6500 patients)

pop = 298.4; popstates = 12.64 + 2.983 + 2.763 + 5.164 + 5.843 + 1.773 + .783 + 5.578
y2006 = y2006*(popstates/pop)
y2006 = c(y2006[1],   ################ <1
          mean(y2006[2:5]),   ################ 1-4
          mean(y2006[6:10]),   ################ 5-9
          mean(y2006[11:15]),  ################ 10-14
          mean(y2006[16:20]),  ################ 15-19
          mean(y2006[21:25]),  ################ 20-24
          mean(y2006[26:30]),  ################ 25-29
          mean(y2006[31:40]),  ################ 30-39
          mean(y2006[41:65]),  ################ 40-64
          mean(y2006[66:80]))  ################ 65-79
later = rbind(c(0.10,0.60,0.53,0.53,0.41,0.41,0.20,0.20,0.16,0.07),
              c(0.07,1.36,1.47,1.47,0.35,0.35,0.12,0.12,0.09,0.02),
              c(0.56,0.70,1.66,1.66,1.69,1.69,0.45,0.45,0.15,0.08),
              c(0.68,1.57,2.18,2.18,1.78,1.78,0.71,0.71,0.19,0.06),
              c(0.07,0.24,0.17,0.17,0.28,0.28,0.11,0.11,0.08,0.04),
              c(0.10,0.23,0.13,0.13,0.05,0.05,0.06,0.06,0.06,0.04),
              c(0.05,0.21,0.12,0.12,0.65,0.65,0.16,0.16,0.07,0.03),
              c(0.15,0.39,0.29,0.29,1.09,1.09,0.38,0.38,0.26,0.10))
later[,3] = later[,3]*(y2006[3]/mean(y2006[3:4]))
later[,4] = later[,4]*(y2006[4]/mean(y2006[3:4]))
later[,5] = later[,5]*(y2006[5]/mean(y2006[5:6]))
later[,6] = later[,6]*(y2006[6]/mean(y2006[5:6]))
later[,7] = later[,7]*(y2006[7]/mean(y2006[7:8]))
later[,8] = later[,8]*(y2006[8]/mean(y2006[7:8]))

append = matrix(NA,10,17)


append[,1:8] = interyears[,1]
append[,9] = log(y2006)
for (i in 10:17){
  append[,i] = log(later[i-9,])
}

vals.obs = cbind(vals.obs,append)
#save(vals.obs,file='vals.obs.Rdata') #### THIS IS THE INCIDENCE DATA
load('vals.obs.Rdata')

##################################################
#### Function for steady-state dynamics ##########
##################################################
setwd('~/Google drive/mumps/data')
contact = read.csv('contact10x10.csv',header=F)
contact = as.matrix(contact)

### divide by pop in each age group to normalize contacat matrix from Britain 2000s (mossong) to age structure of simulated pop
ageweights = c(6.2/5, 4*6.2/5, 5.6, 5.8, 6.3, 6.8, 6.8, 6.5+6.6, 7.3+7.3+6.5+5.7+6.0, 4.8+3.9+3.2)
ageweights = (ageweights/sum(ageweights))
for (i in 1:10){
  contact[i,] = contact[i,]/ageweights
  contact[i,] = contact[i,]*c(1,4,rep(5,5),10,25,15)/80
}

for (i in 1:10) for (j in 1:10){
  contact[i,j] = mean(c(contact[i,j],contact[j,i]))
}


sir.fn = function(t,y,parms,contact){
  
  alpha = 1/(365*c(1,4,5,5,5,5,5,10,25,15))
  
  b = parms[1:10]
  beta = parms[11]
  
  S = y[1:10]
  E = y[11:20]
  I = y[21:30]
  R = y[31:40]
  
  N = c()
  for (i in 1:10){
    N[i] = sum(y[seq(i,40,10)])
  }
  
  lambda = rep(0,10)
  for (i in 1:10){
    lambda[i] = beta*sum(contact[i,]*I/N)
  }
  
  mu = c(1/80,rep(0,9))/365
  #mu = c(1,rep(0,9))/365
  sigma = 1/17
  gamma = 1/5
  
  dS = mu*sum(N) - (lambda*b + alpha)*S
  dS[2:10] = dS[2:10] + alpha[1:9]*S[1:9]
  
  dE = lambda*b*S - (alpha + sigma)*E
  dE[2:10] = dE[2:10] + alpha[1:9]*E[1:9]
  
  dI = sigma*E - (alpha + gamma)*I
  dI[2:10] = dI[2:10] + alpha[1:9]*I[1:9]
  
  dR = gamma*I - alpha*R
  dR[2:10] = dR[2:10] + alpha[1:9]*R[1:9]
  
  return(list(c(dS,dE,dI,dR)))#,dY))) 
}

init = c(1e4*c(1,4,rep(5,5),10,25,15),rep(0,10),rep(100*c(1,4,rep(5,5),10,25,15)),rep(0,10))
vals.1967 = exp(vals.obs[,which(years==1967)]) ### want to fit against this as if it represents steady state pre-vaccine

######### objective function (weighted least squares); use for fitting

sse.fn = function(pars,contact,init,vals.in){
  
  out = ode(y=init,times=(0:70)*365,func=sir.fn,parms=exp(pars[1:11]),contact=contact,method='ode45')[,2:41]
  ns = c()
  for (i in 1:10){
    ns[i] = sum(out[70,seq(i,40,10)])
  }
  prev = out[70,21:30]/ns
  inc = 1e5*prev/(5/365)
  plot(inc*exp(pars[12]),type='l')
  
  error = vals.in - inc*exp(pars[12])
  w = 1/(vals.in + 1)
  return(sum(w*(error^2)))
}

##### TO REPLICATE FITTING FROM INITIAL INPUTS
#pars = c(-12,-1.26,0.55,0.32,0.68,-0.15,-0.74,-0.96,-0.05,0.82,-2.94,-3.53)### for initializing ### 1e3 maxit on this set
#pars = c(-11.1264473,-1.1570564,0.4946110,0.3455874,0.6445062,-0.1626367,-0.7396088,-0.973156,-0.1169686,0.6921156,-2.9823241,-3.2603348) #results from first set
#pars = c(-2.95581247,-1.19720573,0.46297079,0.33415381,0.64976963,-0.11103997,-0.73384064,-0.95729825,-0.03456376,0.83948563,-2.93517025,-3.26592167) #second
pars = c(-2.961813861,-1.201787799,0.472963764,0.342473591,0.666359957,-0.134920864,-0.716061260,-0.905156718,0.003439328,0.930980029,-2.937777770,-3.266766357) #third and essentially last (converges, just need hessian)
tryopt = optim(pars,fn=sse.fn,contact=contact,init=init,vals.in=vals.1967,control=list(trace=2,maxit=5e2,abstol=1e-2),hessian=F)### 
exp(opt1967$par[12])/0.7
#opt1967 = tryopt
#save(opt1967,file='opt1967.Rdata')

#install.packages('numDeriv')
#library(numDeriv)
#hessTry = hessian(func=sse.fn,x=pars,init=init,vals.in=vals.1967,contact=contact,method='Richardson',method.args=list(zero.tol=1))

###########################################################
###########################################################
##### Extension to try with natural waning of protection ##
###########################################################
###########################################################

sirW.fn = function(t,y,parms,contact){
  
  alpha = 1/(365*c(1,4,5,5,5,5,5,10,25,15))
  
  b = parms[1:10]
  beta = parms[11]
  wan = parms[12]
  
  S = y[1:10]
  E = y[11:20]
  I = y[21:30]
  R = y[31:40]
  
  N = c()
  for (i in 1:10){
    N[i] = sum(y[seq(i,40,10)])
  }
  
  lambda = rep(0,10)
  for (i in 1:10){
    lambda[i] = beta*sum(contact[i,]*I/N)
  }
  
  mu = c(1/80,rep(0,9))/365
  #mu = c(1,rep(0,9))/365
  sigma = 1/17
  gamma = 1/5
  
  dS = mu*sum(N) - (lambda*b + alpha)*S + wan*R
  dS[2:10] = dS[2:10] + alpha[1:9]*S[1:9]
  
  dE = lambda*b*S - (alpha + sigma)*E
  dE[2:10] = dE[2:10] + alpha[1:9]*E[1:9]
  
  dI = sigma*E - (alpha + gamma)*I
  dI[2:10] = dI[2:10] + alpha[1:9]*I[1:9]
  
  dR = gamma*I - (alpha + wan)*R
  dR[2:10] = dR[2:10] + alpha[1:9]*R[1:9]
  
  return(list(c(dS,dE,dI,dR)))#,dY))) 
}

sseW.fn = function(pars,contact,init,vals.in){

  out = ode(y=init,times=(0:70)*365,func=sirW.fn,parms=exp(pars[c(1:11,13)]),contact=contact,method='ode45')[,2:41]
  ns = c()
  for (i in 1:10){
    ns[i] = sum(out[70,seq(i,40,10)])
  }
  
  prev = out[70,21:30]/ns
  inc = 1e5*prev/(5/365)
  
  error = vals.in - inc*exp(pars[12])
  w = 1/(vals.in + 1)
  return(sum(w*(error^2)))
}

#load('opt1967.Rdata')
#library(deSolve)
#parsW = c(opt1967$par,log(1/(20*365)))

#parsW = c(-2.4121,-0.28211,1.3906,1.1828,1.3562,0.4485,-0.1173,-0.67214,-0.2414,-0.15077,-3.85654,-3.28982,-11.51724)
parsW = c(-1.983,-0.217,1.448,1.2375,1.4179,0.4678,-0.166,-0.6248,-0.5586,-0.1243,-3.9299,-3.2800,-11.4670)
#tryoptW = optim(parsW,fn=sseW.fn,contact=contact,init=init,vals.in=vals.1967,control=list(trace=2,maxit=5e3,abstol=1e-2),hessian=F)### 
save(tryoptW,file='tryoptW.Rdata')
exp(-tryoptW$par[13])/365
cbind(tryopt$par,tryoptW$par)
pexp(15,1/266.396)
pexp(15,mean(lambda))

aicBase = 2*length(opt1967$par) + length(vals.in)*log(opt1967$value/length(vals.in)) ### -46.10088
aicWane = 2*length(tryoptW$par) + length(vals.in)*log(tryoptW$value/length(vals.in)) ### -43.38043



######################################################################################## deltaAIC = -2.647523; p = 0.03541323

mean(startR[3:4]/N[3:4])
mean(startR[5:6]/N[5:6])
mean(startR[7:8]/N[7:8])

###########################################################
###########################################################
mean(1-exp(est0)); quantile(1-exp(est0),c(0.025,0.975))

########### get nit posterior samples so you can exclude 1 each for 95%CI
nit = 500
#pars = mvrnorm(2e3,mu=opt1967$par,Sigma=matrix(0,12,12))
pars = opt1967$par
  
#################################################################
#### Translating this into starting conditions for model run ####
#################################################################

### for i in 1:nit, generate 
startStates = ode(y=init,times=(0:100)*365,func=sir.fn,parms=exp(pars[1:11]),contact=contact,method='ode45')[,2:41]
startS = startStates[100,1:10]
N = c()
for (i in 1:10){
  N[i] = sum(startStates[100,seq(i,40,10)])
}
startR = N - startStates[100,1:10]
startI = startStates[100,21:30]

b = exp(pars[1:10])
beta = exp(pars[11])

lambda = rep(0,10)
for (i in 1:10){
  lambda[i] = beta*sum(contact[i,]*startI/N)
}


lambda0 = lambda*b
beta0 = b*beta

####################################################################################
#### Getting other conditions (wr, protection) #### run meta-analysis code first ###
####################################################################################

x = seq(-1,5,0.01)
nit = 1e4
omega = c()
for (i in 1:nit){
  dat = (1-exp(est[i,]))/(1-exp(est0[i]))
  omega[i] = optim(0.5,fn=opt.fn,dat=dat,x=x)$par
  print(i)
}
q95fn(1/omega)
range(1-pv)
set.seed(1)
rr = exp(est0) ############### vaccine combined effect on infection and progression

pv = 2.7/runif(nit,7,10) ################### symptom probability if vaccinated
prot = 0.7*rr/pv ######## vaccine protection against infection (1-ve)
repCas = exp(opt1967$par[12])/0.7 ########## reporting of cases (as opposed to reporting of infections)
prot = 0.7*rr/pv
quantile(1 - pv/0.7,c(0.5,0.025,0.975))
pars.input = cbind(1,omega[1:nit],pv,prot[1:nit])
#pars.input = cbind(1,mean(omega),pv,mean(prot))

##################################################
#### Function to compute prevalence of states ####
#### (given FoI inputs) ##########################
##################################################


sr.fn = function(t,y,parms,v1,v2){
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
  Y = y[41:50]
  N = sum(y[1:40])
  
  mu = c(1/80,rep(0,9))/365
  
  dS = mu*N - (lambda + alpha)*S
  dS[2:10] = dS[2:10] + alpha[1:9]*S[1:9]*(1-v1[1:9])
  
  dR = lambda*(S + P) - alpha*R
  dR[2:10] = dR[2:10] + alpha[1:9]*R[1:9]
  
  dV = -(alpha + w)*V
  dV[2:10] = dV[2:10] + alpha[1:9]*(V[1:9] + nu*v1[1:9]*S[1:9] + rho*nu*v2[1:9]*P[1:9])
  
  dP = -(alpha + lambda)*P + w*V
  dP[2:10] = dP[2:10] + alpha[1:9]*(P[1:9]*(1 - rho*nu*v2[1:9]) + (1-nu)*v1[1:9]*S[1:9])
  
  dY = lambda*(p*S + pv*P)
  
  return(list(c(dS,dR,dV,dP,dY)))  
}

load('d1ts.Rdata'); d1ts = c(0,d1ts)
load('d2ts.Rdata'); d2ts = c(0,d2ts)

init = c(startS,startR,rep(0,30))
#setwd('~/Google drive/mumps/data')

v1 = rbind(d1ts,matrix(0,9,length(d1ts)))
v2 = rbind(rep(0,length(d2ts)),d2ts,matrix(0,8,length(d2ts)))

i = 1
out = ode(y=init,times=(0:1)*365,func=sr.fn,parms=c(pars.input[i,],lambda0),v1=v1[,1],v2=v2[,1])[2,2:51]
plot((out[41:50]/N)*1e5*repCas[1]) #### all good to go

init1 = out; init1[41:50] = 0

sqerr.fn = function(lambda.in,pars.input,init.in,v1.in,v2.in,vals.obs.in,N,rep){
  out = ode(y=init.in,times=(0:1)*365,func=sr.fn,parms=c(pars.input,exp(lambda.in)),v1=v1.in,v2=v2.in)[2,2:51]
  err = exp(vals.obs.in) - 1e5*(out[41:50]/N)*rep
  w = 1/(exp(vals.obs.in) + 1)
  return(sum(w*(err^2)))    
}

nit = 500
lambda = list()
init1 = list()

preds = S = R = V = P = array(NA,dim=c(nit,100,10))
S[,1,] = init[1:10]; R[,1,] = init[11:20]; V[,1,] = init[21:30]; P[,1,] = init[31:40]; preds[,1,] = init[41:50]/N
for (i in 1:nit){
  lambda[[i]] = list(); lambda[[i]][[1]] = opt1967; lambda[[i]][[1]]$par = log(lambda0);
  init1[[i]] = list(); init1[[i]][[1]] = init;
  for (t in 2:49){
    lambda[[i]][[t]] = optim(par=lambda[[i]][[t-1]]$par,fn=sqerr.fn,pars.input=pars.input[i,],
                             init.in=init1[[i]][[t-1]],v1.in=v1[,t],v2.in=v2[,t],
                             vals.obs.in=vals.obs[,7+t],N=N,rep=repCas,control=list(trace=0,maxit=1e5),hessian=F) 
    init1[[i]][[t]] = ode(y=init1[[i]][[t-1]],times=(0:1)*365,func=sr.fn,
                          parms=c(pars.input[i,],exp(lambda[[i]][[t]]$par)),v1=v1[,t],v2=v2[,t])[2,2:51]
    
    S[i,t,] = init1[[i]][[t]][1:10]
    R[i,t,] = init1[[i]][[t]][11:20]
    V[i,t,] = init1[[i]][[t]][21:30]
    P[i,t,] = init1[[i]][[t]][31:40]
    
    preds[i,t,] = init1[[i]][[t]][41:50]/N
    init1[[i]][[t]][41:50] = 0
    print(c(lambda[[i]][[t]]$value,i,t))
  }
  output = list(lambda,init1,S,R,V,P,preds)
  save(output,file='output.Rdata')
}

#####################################################################
#### Fit when attributing changes in part to decline in reporting ###
#### (fix at 1% decline per year) ###############################
#####################################################################
#set.seed(10101)
set.seed(20202)
rep1 = rep(repCas[1],49)
for (i in 2:49){
  rep1[i] = rep1[i-1]*0.98 ##99 for 1%
}

sqerr.fn = function(lambda.in,pars.input,init.in,v1.in,v2.in,vals.obs.in,N,rep){
  out = ode(y=init.in,times=(0:1)*365,func=sr.fn,parms=c(pars.input,exp(lambda.in)),v1=v1.in,v2=v2.in)[2,2:51]
  err = exp(vals.obs.in) - 1e5*(out[41:50]/N)*rep
  w = 1/(exp(vals.obs.in) + 1)
  return(sum(w*(err^2)))
}

lambda = list()
init1 = list()

preds = S = R = V = P = array(NA,dim=c(nit,100,10))
S[,1,] = init[1:10]; R[,1,] = init[11:20]; V[,1,] = init[21:30]; P[,1,] = init[31:40]; preds[,1,] = init[41:50]/N
for (i in 1:nit){
  lambda[[i]] = list(); lambda[[i]][[1]] = opt1967; lambda[[i]][[1]]$par = log(lambda0);
  init1[[i]] = list(); init1[[i]][[1]] = init;
  for (t in 2:49){
    lambda[[i]][[t]] = optim(par=lambda[[i]][[t-1]]$par,fn=sqerr.fn,pars.input=pars.input[i,],
                             init.in=init1[[i]][[t-1]],v1.in=v1[,t],v2.in=v2[,t],
                             vals.obs.in=vals.obs[,7+t],N=N,rep=rep1[t],control=list(trace=0,maxit=1e5),hessian=F) 
    init1[[i]][[t]] = ode(y=init1[[i]][[t-1]],times=(0:1)*365,func=sr.fn,
                          parms=c(pars.input[i,],exp(lambda[[i]][[t]]$par)),v1=v1[,t],v2=v2[,t])[2,2:51]
    
    S[i,t,] = init1[[i]][[t]][1:10]
    R[i,t,] = init1[[i]][[t]][11:20]
    V[i,t,] = init1[[i]][[t]][21:30]
    P[i,t,] = init1[[i]][[t]][31:40]
    
    preds[i,t,] = init1[[i]][[t]][41:50]*rep1[i]/N
    init1[[i]][[t]][41:50] = 0
    print(c(lambda[[i]][[t]]$value,i,t))
  }
  outputRep2 = list(lambda,init1,S,R,V,P,preds)
  save(outputRep2,file='outputRep2.Rdata')
}


load('outputRep.Rdata')
sum(is.na(outputRep[[3]][,2,2])==F)

##################################################################
##################################################################
#### Repeat with no waning of vaccine effectiveness ##############
##################################################################
##################################################################

sr.fn = function(t,y,parms,v1,v2){
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
  Y = y[41:50]
  N = sum(y[1:40])
  
  mu = c(1/80,rep(0,9))/365
  
  dS = mu*N - (lambda + alpha)*S
  dS[2:10] = dS[2:10] + alpha[1:9]*S[1:9]*(1-v1[1:9])
  
  dR = lambda*(S + P) - alpha*R
  dR[2:10] = dR[2:10] + alpha[1:9]*R[1:9]
  
  dV = -(alpha + w)*V
  dV[2:10] = dV[2:10] + alpha[1:9]*(V[1:9] + nu*v1[1:9]*S[1:9] + rho*nu*v2[1:9]*P[1:9])
  
  dP = -(alpha + lambda)*P + w*V
  dP[2:10] = dP[2:10] + alpha[1:9]*(P[1:9]*(1 - rho*nu*v2[1:9]) + (1-nu)*v1[1:9]*S[1:9])
  
  dY = lambda*(p*S + pv*P)
  
  return(list(c(dS,dR,dV,dP,dY)))  
}

load('d1ts.Rdata'); d1ts = c(0,d1ts)
load('d2ts.Rdata'); d2ts = c(0,d2ts)

init = c(startS,startR,rep(0,30))
#setwd('~/Google drive/mumps/data')

v1 = rbind(d1ts,matrix(0,9,length(d1ts)))
v2 = rbind(rep(0,length(d2ts)),d2ts,matrix(0,8,length(d2ts)))

i = 1
out = ode(y=init,times=(0:1)*365,func=sr.fn,parms=c(pars.input[i,],lambda0),v1=v1[,1],v2=v2[,1])[2,2:51]
plot((out[41:50]/N)*1e5*repCas[1]) #### all good to go

init1 = out; init1[41:50] = 0

sqerr.fn = function(lambda.in,pars.input,init.in,v1.in,v2.in,vals.obs.in,N,rep){
  out = ode(y=init.in,times=(0:1)*365,func=sr.fn,parms=c(pars.input,exp(lambda.in)),v1=v1.in,v2=v2.in)[2,2:51]
  err = exp(vals.obs.in) - 1e5*(out[41:50]/N)*rep
  w = 1/(exp(vals.obs.in) + 1)
  return(sum(w*(err^2)))    
}

nit = 500
lambda = list()
init1 = list()

preds = S = R = V = P = array(NA,dim=c(nit,100,10))
S[,1,] = init[1:10]; R[,1,] = init[11:20]; V[,1,] = init[21:30]; P[,1,] = init[31:40]; preds[,1,] = init[41:50]/N
for (i in 1:nit){
  lambda[[i]] = list(); lambda[[i]][[1]] = opt1967; lambda[[i]][[1]]$par = log(lambda0);
  init1[[i]] = list(); init1[[i]][[1]] = init;
  for (t in 2:49){
    lambda[[i]][[t]] = optim(par=lambda[[i]][[t-1]]$par,fn=sqerr.fn,pars.input=pars.input[i,],
                             init.in=init1[[i]][[t-1]],v1.in=v1[,t],v2.in=v2[,t],
                             vals.obs.in=vals.obs[,7+t],N=N,rep=repCas,control=list(trace=0,maxit=1e5),hessian=F) 
    init1[[i]][[t]] = ode(y=init1[[i]][[t-1]],times=(0:1)*365,func=sr.fn,
                          parms=c(pars.input[i,],exp(lambda[[i]][[t]]$par)),v1=v1[,t],v2=v2[,t])[2,2:51]
    
    S[i,t,] = init1[[i]][[t]][1:10]
    R[i,t,] = init1[[i]][[t]][11:20]
    V[i,t,] = init1[[i]][[t]][21:30]
    P[i,t,] = init1[[i]][[t]][31:40]
    
    preds[i,t,] = init1[[i]][[t]][41:50]/N
    init1[[i]][[t]][41:50] = 0
    print(c(lambda[[i]][[t]]$value,i,t))
  }
  outputNoWan = list(lambda,init1,S,R,V,P,preds)
  save(outputNoWan,file='outputNoWan.Rdata')
}
 

##################################################################
##################################################################


##################################################################
##################################################################
#### Calculate degree of vax mismatch for RE=1 at all times ######
##################################################################
##################################################################


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


init2 = c(startS,startR,rep(0,20))
curState = array(NA,dim=c(1000,49,40))
for (j in 1:1000){
  curState[j,1,] = init2
  for (t in 2:49){
    curState[j,t,] = ode(y=curState[j,t-1,],parms=c(1,omega[j],pv[j],rr[j],(1/repCas)*exp(vals.obs[,t+6])/(1e5*365)),times=(0:1)*365,
                         v1=v1[,t],v2=v2[,t],func=srVax.fn)[2,2:41]
  }
  print(j)
}

multNW = array(NA,dim=c(1000,49,10))
for (j in 1:500){
  for (t in 2:49){
    Lambda = (1/repCas)*exp(vals.obs[,t+6])/(1e5*365)
    
    lambda = Lambda*ns/(0.7*curState[j,t,1:10] + pv[j]*curState[j,t,31:40])
    prev = lambda*(curState[j,t,1:10] + curState[j,t,31:40])*5
    for (i in 1:10){
      multNW[j,t,i] = lambda[i]/sum(beta*exp(opt1967$par[i])*contact[i,]*prev/ns)
    }
  }
}

plot(multNW[j,,5],type='l',ylim=c(0,2))
for (j in 1:200){
  lines(multNW[j,,5]/multNW[j,2,5])
  lines(alphas[j,,5],col='red')
}




multNW[j,2,5]
#save(multNW,file='multNW.Rdata')
#contact[i,]*
#setwd('~/Google drive/mumps/data')
#load('multNW.Rdata')
# 
# betaMat = matrix(NA,10,10)
# for (i in 1:10){
#   betaMat[i,] = contact[i,]*beta*exp(opt1967$par[i])
# }
# r0base = max(Re(eigen(betaMat*5)$values))
# 
# mult = c(2.5,3,3.5,4,r0base,5)/r0base
# 
# escapeSolve = function(v,S,V,P,ns,betaMat){
#   suscepts = (S+P+exp(v)*V)/ns
#   inputMat = matrix(NA,10,10)
#   for (i in 1:10){
#     inputMat[i,] = betaMat[i,]*suscepts[i]
#   }
#   r0 = max(Re(eigen(inputMat)$values)*5)
#   return((r0-1)^2)
# }
# 
# r0input = v = array(NA,dim=c(6,200,49))
# for (j in 1:200) for (t in 2:49) for (a in 1:6){
#   S.in = curState[j,t,1:10]
#   P.in = curState[j,t,31:40]
#   V.in = curState[j,t,21:30]
#   
#   betaMat.in = betaMat*mult[a]
#   suscepts = S.in+P.in
#   inputMat = matrix(NA,10,10)
#   for (i in 1:10){
#     inputMat[i,] = betaMat.in[i,]*(suscepts[i]/ns[i])
#   }
#   
#   r0input[a,j,t] = max(Re(eigen(inputMat*5)$values))
#   if (r0input[a,j,t]<1){
#     v[a,j,t] = exp(optim(-0.65,fn=escapeSolve,S=S.in,V=V.in,P=P.in,ns=ns,betaMat=betaMat.in,control=list(maxit=1e5))$par)
#   } else{
#     v[a,j,t] = 0
#   }
#   print(c(a,j,t))
# }
# 
# a = 6
# plot(v[a,1,10:49],ylim=c(0,1),type='n')#; lines(v[2,1,10:49]); lines(v[3,1,10:49])
# for (j in 1:200){
#   lines(v[a,j,10:49])
# }
# 
# escape = array(NA,dim=c(6,3,49))
# for (a in 1:6) for (t in 2:49){
#   escape[a,,t] = quantile(v[a,,t],c(0.5,0.025,0.975))
# }


################################################
################################################
###### MODELING OUTBREAK IN VAX POP ############
################################################
################################################

sir.fn = function(t,y,parms,contact,v1,v2,mult,vaxrr){
  
  alpha = 1/(365*c(1,4,5,5,5,5,5,10,25,15))
  
  b = parms[1:10]
  beta = parms[11]
  rep = parms[12]
  rho = 1
  nu = 1 - parms[13]
  pv = parms[14]
  p = 0.7
  w = parms[15]/365
  
  S = y[1:10]
  E = y[11:20]
  I = y[21:30]
  R = y[31:40]
  V = y[41:50]
  P = y[51:60]
  Ep = y[61:70]
  Y = y[71:80]
  
  N = c()
  for (i in 1:10){
    N[i] = sum(y[seq(i,70,10)])
  }
  
  lambda = rep(0,10)
  for (i in 1:10){
    lambda[i] = mult[i]*b[i]*beta*sum(contact[i,]*I/N)
  }
  
  mu = c(1/80,rep(0,9))/365
  
  sigma = 1/17
  gamma = 1/5
  
  dS = mu*sum(N) - (lambda + alpha)*S
  dS[2:10] = dS[2:10] + alpha[1:9]*S[1:9]*(1-v1[1:9])
  
  dE = lambda*S - (alpha + sigma)*E
  dE[2:10] = dE[2:10] + alpha[1:9]*E[1:9]
  
  dI = sigma*(E + Ep) - (alpha + gamma)*I
  dI[2:10] = dI[2:10] + alpha[1:9]*I[1:9]
  
  dR = gamma*I - alpha*R
  dR[2:10] = dR[2:10] + alpha[1:9]*R[1:9]
  
  dV = -(alpha + w + vaxrr*lambda)*V
  dV[2:10] = dV[2:10] + alpha[1:9]*(V[1:9] + nu*v1[1:9]*S[1:9] + rho*nu*v2[1:9]*P[1:9])
  
  dP = -(alpha + lambda)*P + w*V
  dP[2:10] = dP[2:10] + alpha[1:9]*(P[1:9]*(1 - rho*nu*v2[1:9]) + (1-nu)*v1[1:9]*S[1:9])
  
  dEp = lambda*(P + vaxrr*V) - (alpha + sigma)*Ep
  dEp[2:10] = dEp[2:10] + alpha[1:9]*Ep[1:9]
  
  dY = sigma*(p*E + pv*Ep)
  
  return(list(c(dS,dE,dI,dR,dV,dP,dEp,dY)))
}

##### Reading in inputs

setwd('~/Google drive/mumps/data')
load('output.Rdata')
lambda = output[[1]]; init1 = output[[2]]; S = output[[3]]; R = output[[4]]; V = output[[5]]; P = output[[6]]; preds = output[[7]]

##### Getting alpha from Lambda and modeled prev

I = array(NA,dim=c(200,50,10))
for (j in 1:max(which(is.na(preds[,2,1])==F))) for (i in 1:10) for (t in 2:49){
  I[j,t,i] = exp(lambda[[j]][[t]]$par[i])*(S[j,t,i] + P[j,t,i])*5
}

alphas = array(NA,dim=c(200,50,10))
for (j in 1:max(which(is.na(preds[,2,1])==F))) for (i in 1:10) for (t in 2:49){
  alphas[j,t,i] = exp(lambda[[j]][[t]]$par[i])/(b[i]*beta*sum(contact[i,]*I[j,t,]/N))
}
#### alphas[j,t,] is the modifier on beta

v1 = rbind(d1ts,matrix(0,9,length(d1ts)))
v2 = rbind(rep(0,length(d2ts)),d2ts,matrix(0,8,length(d2ts)))
### v1[,t] and v2[,t] are the inputs

### t=39 for conditions at end of 2005
startS06 = S[,39,]
startR06 = R[,39,]
startV06 = V[,39,]#matrix(0,200,10)#V[,39,]
startP06 = P[,39,]

startS16 = S[,49,]
startR16 = R[,49,]
startV16 = V[,49,]#matrix(0,200,10)#V[,49,]
startP16 = P[,49,]

startS06[,5] = startS06[,5]*(1-1e-4)
startE06 = matrix(0,200,10)
startI06 = matrix(0,200,10); startI06[,5] = startS06[,5]*1e-4

startS16[,5] = startS16[,5]*(1-1e-4)
startE16 = matrix(0,200,10)
startI16 = matrix(0,200,10); startI16[,5] = startS16[,5]*1e-4

start06 = cbind(startS06,startE06,startI06,startR06,startV06,startP06,matrix(0,200,20))
alpha06 = alphas[,39,]
multNW06 = multNW[,39,]

start16 = cbind(startS16,startE16,startI16,startR16,startV16,startP16,matrix(0,200,20))
alpha16 = alphas[,49,]
multNW16 = multNW[,49,]

startS06nw = curState[,39,1:10]; startR06nw = curState[,39,11:20]; startV06nw = curState[,39,21:30]; startP06nw = curState[,39,31:40]
startS06nw[,5] = startS06nw[,5]*(1-1e-4)
startE06nw = matrix(0,200,10);
startI06nw = matrix(0,200,10); startI06nw[,5] = startS06nw[,5]*1e-4
start06nw = cbind(startS06nw,startE06nw,startI06nw,startR06nw,startV06nw,startP06nw,matrix(0,200,20))

startS16nw = curState[,49,1:10]; startR16nw = curState[,49,11:20]; startV16nw = curState[,49,21:30]; startP16nw = curState[,49,31:40]
startS16nw[,5] = startS16nw[,5]*(1-1e-4)
startE16nw = matrix(0,200,10);
startI16nw = matrix(0,200,10); startI16nw[,5] = startS16nw[,5]*1e-4
start16nw = cbind(startS16nw,startE16nw,startI16nw,startR16nw,startV16nw,startP16nw,matrix(0,200,20))
#### compare against predictions with partial vax protection against the introduced strain; with partial vax/nat protection against introduced strain
############# had immunity waned
############# had immunity not waned (use output from other analysis to initialize)
vaxrr = c(0.01,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
pred06rate = pred16rate = matrix(NA,200,10); pred06Mrate = pred16Mrate = array(NA,dim=c(length(vaxrr),200,10))

for (j in 1:200){
  out = ode(y=start06[j,],times=c(0:1)*365,func=sir.fn,parms=c(exp(opt1967$par),rr[j],pv[j],omega[j]),contact=contact,mult=alpha06[j,],v1=v1[,39],v2=v2[,39],vaxrr=0,method='ode45')
  pred06rate[j,] = out[2,71:80] - out[1,71:80]
  for (a in 1:length(vaxrr)){
    #### could also use multNW06 which reflects lambda from mod w/out waning, but is it better to show impact of just immune pressure rather than
    #### immune pressure + contact differences?
    out = ode(y=start06nw[j,],times=c(0:1)*365,func=sir.fn,parms=c(exp(opt1967$par),rr[j],pv[j],0),contact=contact,mult=multNW06[j,],v1=v1[,39],v2=v2[,39],vaxrr=vaxrr[a],method='ode45')
    pred06Mrate[a,j,] = out[2,71:80] - out[1,71:80]
  }
  
#  out = ode(y=start16[j,],times=c(0:1)*365,func=sir.fn,parms=c(exp(opt1967$par),rr[j],pv[j],omega[j]),contact=contact,mult=alpha16[j,],v1=v1[,39],v2=v2[,39],vaxrr=0,method='ode45')
#  pred16rate[j,] = out[2,71:80] - out[1,71:80]
#  for (a in 1:length(vaxrr)){
#    #### same as above goes for multNW16 
#    out = ode(y=start16nw[j,],times=c(0:1)*365,func=sir.fn,parms=c(exp(opt1967$par),rr[j],pv[j],0),contact=contact,mult=multNW16[j,],v1=v1[,39],v2=v2[,39],vaxrr=vaxrr[a],method='ode45')
#    pred16Mrate[a,j,] = out[2,71:80] - out[1,71:80]
#  }
  print(j)
}

#save(pred06rate,file='pred06rate.Rdata')
#save(pred06Mrate,file='pred06Mrate.Rdata')
#save(pred16rate,file='pred16rate.Rdata')
#save(pred16Mrate,file='pred16Mrate.Rdata')


#### stuff below is wrong and needs to be modified but gives gist of what to create
t = 48; lines(exp(vals.obs[,t])/max(exp(vals.obs[,t])),lwd=3,col='red'); round(1e5*pred06Mrate[a,10,]/ns,1)
#lines(apply(exp(vals.obs[,48]),1,mean)/max(apply(exp(vals.obs[,48]),1,mean)),col='red',lwd=3)
durs = c(1,4,rep(5,5),10,25,15)
agemeans = c(0.5,3,seq(7.5,27.5,5),35,52.5,72.5)
outAge = outInc = array(NA,dim=c(length(vaxrr),200))
for (a in 1:8) for (j in 1:200){
  outAge[a,j] = sum(pred06Mrate[a,j,]*durs*agemeans/ns)/sum(durs*pred06Mrate[a,j,]/ns)
  outInc[a,j] = sum(1e5*pred06Mrate[a,j,]*durs/(80*ns))
}
obsAges = totInc =
for (j in 1:200) for (t in 1:11){
  obsAges[j,t] = sum(preds[j,39+t-1,]*agemeans*durs)/sum(durs*preds[j,39,])
  totInc[j,t] = sum(1e5*preds[j,39+t-1,]*durs)/80
}
repAges = repInc = c()
for (t in 1:9){
  repAges[t] = sum(exp(vals.obs[,48+t-1])*durs*agemeans)/sum(durs*exp(vals.obs[,48+t-1]))
  repInc[t] = sum(exp(vals.obs[,48+t-1])*durs)/80
}



plot(1,xlim=c(10,30),ylim=c(-7,7))
for (j in 1:200) for (a in 1:8){
  points(x=outAge[a,j],y=log(repCas*outInc[a,j]),col=rgb(a/8,0,1-a/8,0.25),pch=16,cex=0.5)
}

for (j in 1:200) for (t in 1:11){
  points(x=obsAges[j,t],y=log(totInc[j,t]),col=rgb(154/255+runif(1,-0.1,0.1),205/255+runif(1,-0.1,0.1),50/255+runif(1,-0.1,0.1),0.25),pch=16,cex=0.5)
}

for (t in 1:9){
  points(x=repAges[t],y=log(repInc[t]))
}




obsAges = totInc = c()
for (t in 48:56){
  obsAges[t-47] = sum(exp(vals.obs[,t])*agemeans)/sum(exp(vals.obs[,t]))
  totInc[t-47] = sum(exp(vals.obs[,t])*c(1,4,rep(5,5),10,25,15))/80
}


for (i in 48:56){
  lines(exp(vals.obs[,i])/max(exp(vals.obs[,i])),col='red',lwd=3)
}




###### just recover true distributions of cases to compare by multiplying rates by population sizes (by age); can then plot densities side by side

pop06 = c((20071/5)*c(1,4),19606,21145,20730,20971,19561,20471,21052,23056,22123,19496,16490,12589,18463,12971,4860) ### from 2004
######## arranged as: 0-1, 1-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-74, 75-84, 85+

pop16 = c((22358/5)*c(1,4),21623,20984,20243,21810,22195,21858,20543,20250,20926,22376,21649,18761,15621,10987,7761,5600,6822)
######## as: 0-1,1-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70-74, 75-79, 80-84, 85+

set.seed(1)
ages06 = ages16 = obs06 = obs16 = c()
for (i in 1:200){
  ages06 = c(ages06,sample(1:length(pop06),5000,prob=pop06*pred06rate[i,c(1:7,8,8,9,9,9,9,9,10,10,10)]/ns[c(1:7,8,8,9,9,9,9,9,10,10,10)],replace=T))
  ages16 = c(ages16,sample(1:length(pop16),5000,prob=pop16*pred16rate[i,c(1:7,rep(8,2),rep(9,5),rep(10,5))]/ns[c(1:7,rep(8,2),rep(9,5),rep(10,5))],replace=T))
  obs06 = c(obs06,sample(1:length(pop16),5000,prob=pop16*(exp(vals.obs[,56]))[c(1:7,rep(8,2),rep(9,5),rep(10,5))],replace=T))
 # obs16 = c(obs16,sample(1:length(pop06),5000,prob=pop06*exp(vals.obs[,48])[c(1:7,8,8,9,9,9,9,9,10,10,10)],replace=T))
}

ages06M = ages16M = list()
for (a in 1:8){
  temp06 = c()
  temp16 = c()
  for (i in 1:200){
    temp06 = c(temp06,sample(1:length(pop06),5000,prob=pop06*pred06Mrate[a,i,c(1:7,8,8,9,9,9,9,9,10,10,10)]/ns[c(1:7,8,8,9,9,9,9,9,10,10,10)],replace=T))
    temp16 = c(temp16,sample(1:length(pop16),5000,prob=pop16*pred16Mrate[a,i,c(1:7,rep(8,2),rep(9,5),rep(10,5))]/ns[c(1:7,rep(8,2),rep(9,5),rep(10,5))],replace=T))
  }
  ages06M[[a]] = temp06
  ages16M[[a]] = temp16
  print(a)
}

plot(table(obs06),type='l',ylim=c(0,4e5)) #### YES THIS WAY
lines(table(ages16),col='grey',type='l')

lines(table(ages16M[[2]]),col='pink',type='l')
lines(table(ages16M[[3]]),col='red',type='l')
lines(table(ages16M[[4]]),col='dark red',type='l')
#lines(table(ages06M[5]]),col='purple',type='l')

1e5*pred06Mrate[4,1,]/ns

1e5*pred06Mrate[5,1,]/ns

relMax16
relMax16M

par(mfrow=c(2,1))

a = 5

for (i in 1:200){
  if (leaveout[i]==F){
    lines(log(1e5*pred06rate[i,]/ns),x=c(0.5,3,seq(7.5,27.5,5),35,52.5,72.5))
    lines(log(1e5*pred06Mrate[a,i,]/ns),x=c(0.5,3,seq(7.5,27.5,5),35,52.5,72.5),col='red') 
  }
}

plot(pred06rate[1,]/ns,type='l')
plot(exp(vals.obs[,48]),type='l',ylim=c(0,100))
lines(1e5*pred06rate[1,]/ns,col='red')
lines(1e5*pred06Mrate[2,1,]/ns,col='blue')
#########################################################################
#########################################################################
#########  Assume no waning, just mismatch ##############################
#########################################################################
#########################################################################

#save(curState,file='curState.Rdata')



#########################################################################
#########################################################################
#########  Simluating incidence assuming vax waning did occur: main difference is in the attack rate (very slight)
#########################################################################
#########################################################################



#save(curState,file='curState.Rdata')




##################################################################
##################################################################
#### Function to predict incidence if a new strain introduced ####
##################################################################
##################################################################

sirCrossNat.fn = function(t,y,parms,contact,v1,v2){
  
  alpha = 1/(365*c(1,4,5,5,5,5,5,10,25,15))
  
  b = parms[1:10]
  beta = parms[11]
  rep = parms[12]
  rho = 1
  nu = 1 - parms[13]
  pv = parms[14]
  p = 0.7
  w = parms[15]/365
  ant = parms[16]
  modvar = parms[17:26]
  
  S = y[1:10]
  E = y[11:20]
  I = y[21:30]
  R = y[31:40]
  V = y[41:50]
  P = y[51:60]
  
  Eb = y[61:70]
  Epb = y[71:80]
  Ib = y[81:90]
  Rb = y[91:100]
  
  Eab = y[101:110]
  Iab = y[111:120]
  Rab = y[121:130]
  
  Eba = y[131:140]
  Iba = y[141:150]
  
  Ep = y[151:160]
  Y = y[161:170]
  
  N = c()
  for (i in 1:10){
    N[i] = sum(y[seq(i,160,10)])
  }
  
  lambdaA = lambdaB = rep(0,10)
  for (i in 1:10){
    lambdaA[i] = b[i]*beta*sum(modvar[i]*contact[i,]*(I+Iab)/N)
    lambdaB[i] = b[i]*beta*sum(modvar[i]*contact[i,]*(Ib+Iba)/N)
  }
  
  mu = c(1/80,rep(0,9))/365
  
  sigma = 1/17
  gamma = 1/5
  
  dS = mu*sum(N) - (lambdaA + lambdaB + alpha)*S #### added in lambdaB
  dS[2:10] = dS[2:10] + alpha[1:9]*S[1:9]*(1-v1[1:9])
  
  dE = lambdaA*S - (alpha + sigma)*E
  dE[2:10] = dE[2:10] + alpha[1:9]*E[1:9]
  
  dI = sigma*(E + Ep) - (alpha + gamma)*I
  dI[2:10] = dI[2:10] + alpha[1:9]*I[1:9]
  
  dR = gamma*I - (alpha + ant*lambdaB)*R #### added in lambdaB*ant
  dR[2:10] = dR[2:10] + alpha[1:9]*R[1:9]
  
  dV = -(alpha + w + ant*lambdaB)*V #### added in lambdaB*ant
  dV[2:10] = dV[2:10] + alpha[1:9]*(V[1:9] + nu*v1[1:9]*S[1:9] + rho*nu*v2[1:9]*P[1:9])
  
  dP = -(alpha + lambdaA + lambdaB)*P + w*V #### added in lambdaB
  dP[2:10] = dP[2:10] + alpha[1:9]*(P[1:9]*(1 - rho*nu*v2[1:9]) + (1-nu)*v1[1:9]*S[1:9])
  
  dEp = lambdaA*P - (alpha + sigma)*Ep
  dEp[2:10] = dEp[2:10] + alpha[1:9]*Ep[1:9]
  
  ######
  
  dEb = lambdaB*S - (alpha + sigma)*Eb
  dEb[2:10] = dEb[2:10] + alpha[1:9]*Eb[1:9]
  
  dIb = sigma*Eb - (alpha + gamma)*Ib
  dIb[2:10] = dIb[2:10] + alpha[1:9]*Ib[1:9]
  
  dRb = gamma*Ib - (alpha + ant*lambdaA)*Rb  #### added in ant*lambdaA
  dRb[2:10] = dRb[2:10] + alpha[1:9]*Rb[1:9]
  
  #######
  
  dEab = ant*lambdaA*Rb - (alpha + sigma)*Eab
  dEab[2:10] = dEab[2:10] + alpha[1:9]*Eab[1:9]
  
  dIab = sigma*Eab - (alpha + gamma)*Iab
  dIab[2:10] = dIab[2:10] + alpha[1:9]*Iab[1:9]
  
  dRab = gamma*(Iab + Iba) - alpha*Rab
  dRab[2:10] = dRab[2:10] + alpha[1:9]*Rab[1:9]
  
  dEba = ant*lambdaB*R - (alpha + sigma)*Eba
  dEba[2:10] = dEba[2:10] + alpha[1:9]*Eba[1:9]
  
  dEpb = lambdaB*(P+ant*V) - (alpha + sigma)*Epb
  dEpb[2:10] = dEpb[2:10] + alpha[1:9]*Epb[1:9]
  
  dIba = sigma*(Eba + Epb) - (alpha + gamma)*Iba
  dIba[2:10] = dIba[2:10] + alpha[1:9]*Iba[1:9]
  
  ######
  
  dY = rep*sigma*(p*(E + Eb) + pv*(Eba + Eab + Ep + Epb))
  
  return(list(c(dS,dE,dI,dR,dV,dP,dEb,dEpb,dIb,dRb,dEab,dIab,dRab,dEba,dIba,dEp,dY)))
}

setwd('~/Google drive/mumps/data')
load('output.Rdata')
lambda = output[[1]]; init1 = output[[2]]; S = output[[3]]; R = output[[4]]; V = output[[5]]; P = output[[6]]; preds = output[[7]]

### get betas for inputs

I = array(NA,dim=c(200,50,10))
for (j in 1:max(which(is.na(preds[,2,1])==F))) for (i in 1:10) for (t in 2:49){
  I[j,t,i] = exp(lambda[[j]][[t]]$par[i])*(S[j,t,i] + P[j,t,i])*5
}

alphas = array(NA,dim=c(200,50,10))
for (j in 1:max(which(is.na(preds[,2,1])==F))) for (i in 1:10) for (t in 2:49){
  alphas[j,t,i] = exp(lambda[[j]][[t]]$par[i])/(b[i]*beta*sum(contact[i,]*I[j,t,]/N))
}


v1 = rbind(d1ts,matrix(0,9,length(d1ts)))
v2 = rbind(rep(0,length(d2ts)),d2ts,matrix(0,8,length(d2ts)))

init.last = c(startStates[100,],rep(0,130))
saved = matrix(NA,49,170)




j = 13
for (t in 2:49){
  out = ode(y=init.last,times=(1:2)*365,func=sir.fn,method='ode45',
            parms=c(exp(opt1967$par),
                    rr[j],pv[j],omega[j],0.05,alphas[j,t,]),
            v1=v1[,t],v2=v2[,t],contact=contact)[2,2:171]
  saved[t,] = out
  init.last = out
  init.last[161:170] = 0
  print(t)
} ### modify beta inputs

plot(apply(saved[,161:170],1,sum),type='l')


plot(apply(out[,1:160],1,sum))

for (t in 2:50){
  cas[t,] = out[t,161:170] - out[t-1,161:170]
}

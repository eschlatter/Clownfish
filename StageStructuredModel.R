setwd('C:/Users/eschlatter/Desktop/Code/Clownfish')

library(pse)
library(ppcor)

n=55  #should be >50; otherwise error from response variables
K=1000

F=1000
M1=0
M2=1000
JA=1000


# Main Function #################

ClownfishPopSim <- function(hf=0,hm=0,hj=0,bMF=.006854,bMM=.005781,bMJA=.005056,bMJB=.02,surv=.00021,
                            delay_f=6, m_growth=7,delay_j=5,mean_hatch=447.5,diff_hatch=205,
                            P2=1,resp='persist',plot=FALSE){
  
  #reproductive parameters
  # B=clutches_month*(clutch_mean+.5*clutch_diff)  #448
  # S=clutches_month*max(0,clutch_mean-.5*clutch_diff)  #328
  B=mean_hatch+.5*diff_hatch
  S=mean_hatch-.5*diff_hatch
  
  #mortality parameters adjusted for harvest
  MF=bMF+hf*(1-bMF)         
  MM=bMM+hm*(1-bMM)
  MJA=bMJA+hj*(1-bMJA)
  MJB=bMJB+hj*(1-bMJB)
  
  #initial conditions
  JB=JA*P2
  L=(F/K)*surv*(B*M1+S*M2)
  
  #population matrix
  pop=matrix(ncol=6,nrow=n)
  colnames(pop)=c('F','M1','M2','JA','JB','L')
  pop[1,]=c(F,M1,M2,JA,JB,L)
  
  for (i in 2:n){
    
    M=M1+M2
    
    #number of eggs hatched and survived to settlement (able to fill an open anemone spot)
    settlers=(F/K)*surv*(B*M1+S*M2)             
    
    #transition rates between classes
    TF=min((1/delay_f)*(K-F)/(M+.01),1)
    TM=min((1/delay_j)*(K-M)/(JA+.01),1)
    TJAB=min(P2*(K-JA)/(JB+.01),1)                                          #from JB to JA
    TJAL=min((1-P2)*(K-JA)/(settlers*(1-P2)+.01),1)                         #from L to JA -- number, not rate
    TJBL=min((K*P2-JB)/(settlers*P2+.01),1)                                 #from L to JB -- number, not rate
    
    #class populations
    F=max(F-MF*F+(1-MM)*TF*M1+(1-MM)*(1-(1/m_growth))*TF*M2,0)
    M1=max(M1-MM*M1+(1-MM)*(1/m_growth)*M2-(1-MM)*TF*M1,0)
    M2=max(M2-MM*M2-(1-MM)*(1/m_growth)*M2-(1-MM)*(1-(1/m_growth))*TF*M2+(1-MJA)*TM*JA,0)
    JA=max(JA-MJA*JA-(1-MJA)*TM*JA+(1-MJB)*TJAB*JB+TJAL*settlers*(1-P2),0)
    JB=max(JB-MJB*JB-(1-MJB)*TJAB*JB+TJBL*settlers*P2,0)
    L=max(settlers,0)        #number of eggs hatched and survived to settlement
    
    n_JA=(1-MJB)*TJAB*B+TJAL*settlers*(1-P2)
    n_M2=(1-MJA)*TM*JA
    
    new<-t(c(F,M1,M2,JA,JB,L))
    pop[i,]=new
    
  }
  
  pop=as.data.frame(pop)
  colnames(pop)=c('F','M1','M2','JA','JB','L')
  
  
  ##Optional dynamics plot
  
  if(plot==TRUE){
    gen<-c(1:n)
    totalF=pop[,1]
    totalM=pop[,2]+pop[,3]
    totalJA=pop[,4]
    totalJB=pop[,5]
    totalJ=totalJA+totalJB
    
    plot(gen,totalF, ylim=c(0,K), type="l", lwd=2, col=rgb(.5,0,.7),
         xlab='time (months)',ylab='class population')
    title(main='Population Dynamics')
    lines(gen,totalM,type="l",col=4, lwd=2)
    lines(gen,totalJA,type="l",col=2, lwd=2)
    lines(gen,totalJB,type="l",col=1, lwd=2)
    legend("topright",c('F','M','JA','JB'),col=c(rgb(.5,0,.7),4,2,1),lty=c(1,1,1,1),cex=1)
  }
  
  ##Response variables
  
  #change in population over last 50 timesteps
  delta=sum(pop[n,(1:4)])-sum(pop[(n-50),(1:4)])
  if(resp=='delta50'){
    return(delta)}
  
  #total population (minus larvae and J2) at last timestep
  finaltotal=sum(pop[n,(1:4)])+pop[n,5]/P2
  if(resp=='finaltotal'){
    return(finaltotal)}
  
  # returns binary value: does population persist (1) or go extinct (0)?
  # 1 if population is >10 and nondecreasing; 0 otherwise
  # i.e. 1=population persistence; 0=extinction or undetermined in the time of simulation
  if(resp=='persist'){
    if(delta>-.5) (if(finaltotal>10) y=1 else y=0) else y=0
    return(y)}
  
  # persistence measure incorporating undetermined cases
  if(resp=='persist_uncert'){
    if(pop[n,4]+pop[n,5]<1) y=-1
    else(if(delta>-.5) y=1 else y=0)
    return(y)
  }
  
  #returns binary value: is pop at equilibrium?
  if(resp=='equil'){
    if(delta< -.5) equil=0 else equil=1
    return(equil)}
  
  #returns full population dynamics matrix
  if(resp=='matrix'){
    return(pop)}
  
  #returns total of harvested individuals
  if(resp=='totalyield'){
    yield.F=hf*(1-bMF)*sum(pop[,1])
    yield.M=hm*(1-bMM)*(sum(pop[,2])+sum(pop[,3]))
    yield.J=hj*(1-bMJA)*sum(pop[,4])+hj*(1-bMJB)*sum(pop[,5])
    yield=c(yield.F,yield.M,yield.J)
    
    return(sum(yield))}
  
  #returns per-month yield over long term (corrected for %anemones w/2 juveniles)
  if(resp=='yield_monthly'){
    yield.F=hf*(1-bMF)*(sum(pop[(n-50):n,1]))/50
    yield.M=hm*(1-bMM)*(sum(pop[(n-50):n,2])+sum(pop[(n-50):n,3]))/50
    yield.J=hj*(1-bMJA)*sum(pop[(n-50):n,4])/50+hj*(1-bMJB)*sum(pop[(n-50):n,5])/(50*P2)
    yield=c(yield.F,yield.M,yield.J)
    return(sum(yield))}
  
  #returns 2-item vector with persistence and monthly yield
  if(resp=='p_and_y'){
    yield.F=hf*(1-bMF)*(sum(pop[(n-50):n,1]))/50
    yield.M=hm*(1-bMM)*(sum(pop[(n-50):n,2])+sum(pop[(n-50):n,3]))/50
    yield.J=hj*(1-bMJA)*sum(pop[(n-50):n,4])/50+hj*(1-bMJB)*sum(pop[(n-50):n,5])/(50*P2)
    yield=c(yield.F,yield.M,yield.J)
    y=sum(yield)
    
    if(delta>-.5) (if(finaltotal>10) p=1 else p=0) else p=0
    
    return(c(p,y))
  }
}

# Plateau plots ################
#investigate different ranges of harvest strategies (hf, hm, hj) while keeping bio params constant

Harvest <- function(hf=0,hm=0,hj=0){
  return(ClownfishPopSim(hf=hf, hm=hm, hj=hj))
}

HarvestRun<-function(params){
  response=mapply(Harvest,hf=params[,1],hm=params[,2],hj=params[,3])
  response=as.matrix(response)
  response=as.data.frame(cbind(params,response))
  colnames(response)=c('hf','hm','hj','yield')
  return(response)
}

HarvestThreshold <- function(harvest_params,thresh='p'){
  
  # # for yield curve turning point threshold
  if(thresh=='y'){
    yield=mapply(ClownfishPopSim,hf=harvest_params[1],hm=harvest_params[2],hj=harvest_params[3],bMF=.006854,
                 bMM=.005781,bMJA=.005056,bMJB=.02,surv=.00021,delay_f=5,m_growth=7,delay_j=6,mean_hatch=447.5,
                 diff_hatch=205,P2=1,resp='yield_monthly')
  }
  
  # for population persistence threshold
  if(thresh=='p'){
    yield=mapply(ClownfishPopSim,hf=harvest_params[,1],hm=harvest_params[,2],hj=harvest_params[,3],bMF=.006854,
                 bMM=.005781,bMJA=.005056,bMJB=.02,surv=.00021,delay_f=5,m_growth=7,delay_j=6,mean_hatch=447.5,
                 diff_hatch=205,P2=1,resp='persist')
  }
  
  # for population persistence range
  if(thresh=='p_range'){
    yield=mapply(ClownfishPopSim,hf=harvest_params[,1],hm=harvest_params[,2],hj=harvest_params[,3],bMF=.006854,
                 bMM=.005781,bMJA=.005056,bMJB=.02,surv=.00021,delay_f=5,m_growth=7,delay_j=6,mean_hatch=447.5,
                 diff_hatch=205,P2=1,resp='persist_uncert')
  }
  
  # store data with parameter values
  params=data.frame(hf=harvest_params[,1],hm=harvest_params[,2],hj=harvest_params[,3],bMF=.006854,bMM=.005781,
                    bMJA=.005056,bMJB=.02,surv=.00021,delay_f=5,m_growth=7,delay_j=6,mean_hatch=447.5,
                    diff_hatch=205,P2=1)
  results=cbind(params,yield)
  if(thresh=='test') return(results)
  
  #define threshold in terms of population persistence
  if(thresh=='p'){
    threshold=subset(results,yield<1)[1,]
  }
  
  #define range in terms of population persistence
  if(thresh=='p_range'){
    threshold_min=subset(results,yield<1)[1,]
    threshold_max=subset(results,yield<0)[1,]
    return(list(threshold_min,threshold_max))
  }
  
  #define threshold in terms of decreased yield w/increased harvest -- not harvest rate that results in extinction
  if(thresh=='y'){
    results=cbind(results,'decreasing'=NA)
    yield_range=max(results$yield)-min(results$yield)
    for(i in 2:nrow(params)){
      if((results$yield[i]-results$yield[(i-1)])< -(0.01*yield_range))results$decreasing[i]=1
      else results$decreasing[i]=0
    }
    threshold=subset(results,decreasing==1)[1,]
  }
  
  return(threshold[1,])
}

#Plots

par(mfrow=c(1,3))

###Change F harvest rate
fparams=t(rbind(seq(0,1,by=.001),rep(0,1001),rep(0,1001)))
fyield=HarvestRun(fparams)

#persistence threshold
f_threshold=HarvestThreshold(fparams,thresh='p_range')
f_thresh_min=f_threshold[[1]]$hf
f_thresh_max=f_threshold[[2]]$hf
MSY_f=max(subset(fyield,hf<f_thresh_min)$yield)

#plot
plot(fyield$hf,fyield$yield, xlim=c(0,0.15),type="l",xlab="",ylab="")
abline(v=f_thresh_min,col='red')
#abline(v=f_thresh_max,col='blue')
title(main="Female Harvest Yields")

###Change M harvest rate
mparams=t(rbind(0*c(0:1000),0.001*c(0:1000),0*c(0:1000)))
myield=HarvestRun(mparams)

#persistence threshold
m_threshold=HarvestThreshold(mparams,thresh='p_range')
m_thresh_min=m_threshold[[1]]$hm
m_thresh_max=m_threshold[[2]]$hm
MSY_m=max(subset(myield,hm<m_thresh_min)$yield)

plot(myield$hm,myield$yield, xlim=c(0,.15), type="l",xlab="",ylab="")
abline(v=m_thresh_min,col='red')
#abline(v=m_thresh_max,col='blue')
title(main="Male Harvest Yields")

###Change NB harvest rate
jparams=t(rbind(0*c(0:1000),0*c(0:1000),.001*c(0:1000)))
jyield=HarvestRun(jparams)

#persistence threshold
j_threshold=HarvestThreshold(jparams,thresh='p_range')
j_thresh_min=j_threshold[[1]]$hj
j_thresh_max=j_threshold[[2]]$hj
MSY_j=max(subset(jyield,hj<j_thresh_min)$yield)

plot(jyield$hj,jyield$yield, type="l",xlab="",ylab="")
abline(v=j_thresh_min,col='red')
#abline(v=j_thresh_max,col='blue')
title(main="Nonbreeder Harvest Yields")


# PRCC #####################

#define parameter distributions (from file Params.R)

findalpha <- function(mu,sigma2){
  return(mu^2/sigma2-mu^3/sigma2-mu)
}
findbeta <- function(mu,sigma2){
  return(mu^3/sigma2-2*mu^2/sigma2+(1+1/sigma2)*mu-1)
}

bMF=read.csv('bMF.csv')
bMF_a=findalpha(mean(bMF[,1]),var(bMF[,1]))
bMF_b=findbeta(mean(bMF[,1]),var(bMF[,1]))
par(mfrow=c(3,2))
plot(0.001*c(0:100),dbeta(0.001*c(0:100),bMF_a,bMF_b),type='l',ylab='p',xlab='bMF')

bMM=read.csv('bMM.csv')
bMM_a=findalpha(mean(bMM[,1]),var(bMM[,1]))
bMM_b=findbeta(mean(bMM[,1]),var(bMM[,1]))
plot(0.001*c(0:100),dbeta(0.001*c(0:100),bMM_a,bMM_b),type='l',ylab='p',xlab='bMM')

bMJA=read.csv('bMJA.csv')
bMJA_a=findalpha(mean(bMJA[,1]),var(bMJA[,1]))
bMJA_b=findbeta(mean(bMJA[,1]),var(bMJA[,1]))
plot(0.001*c(0:100),dbeta(0.001*c(0:100),bMJA_a,bMJA_b),type='l',ylab='p',xlab='bMJA')

bMJB=read.csv('bMJB.csv')
bMJB_a=findalpha(mean(bMJB[,1]),var(bMJB[,1]))
bMJB_b=findbeta(mean(bMJB[,1]),var(bMJB[,1]))
plot(0.001*c(0:100),dbeta(0.001*c(0:100),bMJB_a,bMJB_b),type='l',ylab='p',xlab='bMJB')

mean_hatch=read.csv('mean_hatch.csv')
mean_hatch_mu=mean(mean_hatch[,1])
mean_hatch_sigma=sd(mean_hatch[,1])
plot(c(0:1000),dnorm(c(0:1000),mean=mean_hatch_mu,sd=mean_hatch_sigma),type='l',ylab='p',xlab='mean_hatch')

#Run LHS and PRCC

#utility function for PRCC to test sensitivity of population size to bio parameters *no harvest values
BioTest=function(params){
  mapply(ClownfishPopSim,hf=0,hm=0,hj=0,bMF=params[,1],bMM=params[,2],bMJA=params[,3],bMJB=params[,4],
         surv=params[,5],delay_f=params[,6], m_growth=params[,7],delay_j=params[,8],mean_hatch=params[,9],
         diff_hatch=params[,10],P2=params[,11])
}

#for one bio param set, find a yield or persistence harvest parameter threshold value *range of harvest values
BioThreshold <- function(bio_params=c(.006854,.005781,.005056,.02,0.87*.00021,5,7,6,447.5,205,1),thresh='msy'){   
  #y:   returns hj value for which monthly yield is greatest
  #p:   returns lowest hj value that results in extinction (according to 'persist' criteria, 
  #        i.e. pop decreasing by more than 0.5 and at total <10)
  #msy: returns highest monthly yield w/p=1
  
  #run base function (get persist and yield output) for all values of hj
  yield=mapply(ClownfishPopSim,hf=0,hm=0,hj=seq(0,1,by=.01),bMF=bio_params[1],bMM=bio_params[2],
               bMJA=bio_params[3],bMJB=bio_params[4],surv=bio_params[5],delay_f=bio_params[6],
               m_growth=bio_params[7],delay_j=bio_params[8],mean_hatch=bio_params[9],diff_hatch=bio_params[10],
               P2=bio_params[11],resp='p_and_y')
  yield=t(yield)
  colnames(yield)=c('p','y')  #get a 2-column by 101-row matrix
  
  # store yield and persistence data with corresponding parameter values in results matrix
  params=data.frame(hf=0,hm=0,hj=seq(0,1,by=.01),bMF=bio_params[1],bMM=bio_params[2],
                    bMJA=bio_params[3],bMJB=bio_params[4],surv=bio_params[5],delay_f=bio_params[6],
                    m_growth=bio_params[7],delay_j=bio_params[8],mean_hatch=bio_params[9],
                    diff_hatch=bio_params[10],P2=bio_params[11])
  results=cbind(params,yield)
  
  #define threshold to return
  if(thresh=='p'){
    threshold=subset(results,p==0)[1,]  #first (lowest) value of hj that results in persistence=0
    return(threshold$hj)}
  
  if(thresh=='y'){
    threshold=results[which(results$y==max(results$y)),] #row of results matrix for which yield is maximized
    return(threshold$hj)}  #return the hj value for that row
  
  if(thresh=='msy'){
    threshold=subset(results,p==1)
    if(nrow(threshold)==0) return(0)  #in case population always goes extinct
    else return(max(threshold$y))}    #return max y for which p=1
}

#utility function for PRCC to test sensitivity of threshold *range of harvest values
ThresholdTest<-function(params){  #takes a matrix (like generated by LHS) with each row a set of bio parms
  params=as.list(as.data.frame(t(params)))
  return(mapply(BioThreshold,params))
}

factors<-c("bMF","bMM","bMJA","bMJB","surv","delay_f","m_growth","delay_j","mean_hatch","diff_hatch","P2")
q<-c("qbeta","qbeta","qbeta","qbeta","qunif","qunif","qunif","qunif","qnorm","qunif","qunif")
q.arg<-list(list(shape1=bMF_a,shape2=bMF_b),list(shape1=bMM_a,shape2=bMM_b),list(shape1=bMJA_a,shape2=bMJA_b),
            list(shape1=1.7,shape2=120), list(min=0,max=0.00042),list(min=1,max=11),list(min=1,max=13),
            list(min=1,max=9),list(mean=mean_hatch_mu,sd=mean_hatch_sigma),list(min=0,max=410),
            list(min=0,max=1))
factor_info=rbind(q,q.arg)
colnames(factor_info)=factors

myLHS<-LHS(BioTest,factors,500,q,q.arg, nboot=50)

#PRCC significance

tlow=qt(.025,488)
thigh=qt(.975,488)
N=500  #number of PRCC runs
p=10   #df (number of params minus one)
alpha1=tlow/sqrt(N-2-p+tlow^2)
alpha2=thigh/sqrt(N-2-p+thigh^2)

#Plots
par(mfrow=c(1,1))

plotecdf(myLHS)
plotscatter(myLHS)

#PRCC plot w/significance bars (make sure alpha1 and alpha2 are correct for situation: i.e. N=500, p=10)
plotprcc(myLHS,ylab='PRCC')
lines(x=c(0,12),y=rep(alpha2,2))
lines(x=c(0,12),y=rep(alpha1,2))

LHSparams=get.data(myLHS)  #show parameter values generated by LHS
LHSresults=get.results(myLHS)
PRCCoutput=cbind(LHSparams,LHSresults)
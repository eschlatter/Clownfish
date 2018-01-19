library(Matrix)

bmF=.006854
bmM=.005781
bmNBA=.005056
bmNBB=.02
mean_hatch=447.5
diff_hatch=205
surv=.00021
delay_f=6
delay_j=5
m_growth=7
hf=0
hm=0
hnb=0

#To do: create a vector lowest to store the lowest rank on each anemone

n=55
ranks=5
K=1000
pop_init=data.frame('F'=1000,'M'=1000,'NBA'=1000,'NBB'=1000,'L'=0)

mF=hf+bmF
mM=hm+bmM
mNBA=hnb+bmNBA
mNBB=hnb+bmNBB

#reproductive parameters
B=mean_hatch+.5*diff_hatch
S=max(0,mean_hatch-.5*diff_hatch)

#initialize matrices to store data
dynamics=array(data=0,dim=c(n,ranks,K),dimnames=list(NULL,c('F','M','NBA','NBB','L'),NULL))
harvest=matrix(nrow=n,ncol=ranks-1)
colnames(harvest)=c('F','M','NBA','NBB')

#allocate initial population to anemones -- fix to require saturation from above?
new_f=sample(c(1:K),pop_init$F,replace=FALSE)
dynamics[1,1,new_f]=1

new_m=sample(c(1:K),pop_init$M,replace=FALSE)
dynamics[1,2,new_m]=1

new_nba=sample(c(1:K),pop_init$NBA,replace=FALSE)
dynamics[1,3,new_nba]=1

new_nbb=sample(c(1:K),pop_init$NBB,replace=FALSE)
dynamics[1,4,new_nbb]=1

#run simulation over time
for(i in 2:n){
  
  fmort=sample(c(1:K),size=mF*K)
  mmort=sample(c(1:K),size=mM*K)
  nbamort=sample(c(1:K),size=mNBA*K)
  nbbmort=sample(c(1:K),size=mNBB*K)
  
  for(anem in 1:K){
    
    #mortality (incl. harvest) +growth
    # if(rbinom(1,1,prob=mF)==1) dynamics[i,'F',anem]=0
    # else if (dynamics[(i-1),'F',anem]!=0) dynamics[i,'F',anem]=dynamics[(i-1),'F',anem]+1
    # 
    # if(rbinom(1,1,prob=mM)==1) dynamics[i,'M',anem]=0
    # else if (dynamics[(i-1),'M',anem]!=0) dynamics[i,'M',anem]=dynamics[(i-1),'M',anem]+1
    # 
    # if(rbinom(1,1,prob=mNBA)==1) dynamics[i,'NBA',anem]=0
    # else if (dynamics[(i-1),'NBA',anem]!=0) dynamics[i,'NBA',anem]=dynamics[(i-1),'NBA',anem]+1
    # 
    # if(rbinom(1,1,prob=mNBB)==1) dynamics[i,'NBB',anem]=0
    # else if (dynamics[(i-1),'NBB',anem]!=0) dynamics[i,'NBB',anem]=dynamics[(i-1),'NBB',anem]+1
    
    if(anem %in% fmort) dynamics[i,'F',anem]=0
    else if (dynamics[(i-1),'F',anem]!=0) dynamics[i,'F',anem]=dynamics[(i-1),'F',anem]+1

    if(anem %in% mmort) dynamics[i,'M',anem]=0
    else if (dynamics[(i-1),'M',anem]!=0) dynamics[i,'M',anem]=dynamics[(i-1),'M',anem]+1

    if(anem %in% nbamort) dynamics[i,'NBA',anem]=0
    else if (dynamics[(i-1),'NBA',anem]!=0) dynamics[i,'NBA',anem]=dynamics[(i-1),'NBA',anem]+1

    if(anem %in% nbbmort) dynamics[i,'NBB',anem]=0
    else if (dynamics[(i-1),'NBB',anem]!=0) dynamics[i,'NBB',anem]=dynamics[(i-1),'NBB',anem]+1
    
    #transition (change the order to prevent cascades in one timestep)
    
    #M to F
    if(dynamics[i,'F',anem]==0){
      if(dynamics[i,'M',anem]>delay_f-1){
        dynamics[i,'F',anem]=1
        dynamics[i,'M',anem]=0
      }
    }
    
    #NBA to M
    if(dynamics[i,'M',anem]==0){
      if(dynamics[i,'NBA',anem]>delay_j-1){
        dynamics[i,'M',anem]=1
        dynamics[i,'NBA',anem]=0
      }
    }
    
    #NBB to NBA
    if(dynamics[i,'NBA',anem]==0){
      if(dynamics[i,'NBB',anem]>0){
        dynamics[i,'NBA',anem]=1
        dynamics[i,'NBB',anem]=0
      }
    }
    
    #reproduction
    if(dynamics[i,1,anem]>0){       #if there's a F in the anemone
      if(dynamics[i,2,anem]>0){     #if there's a M in the anemone
        if(dynamics[i,2,anem]>m_growth-1){
          dynamics[i,'L',anem]=B
        }
        else dynamics[i,'L',anem]=S
      }
    }
  }
  
  #larval survival and settlement
  larvae=surv*sum(dynamics[i,'L',])
  
  #Way 1: each larva only has one chance to settle: if it picks an occupied anem, it dies
  # if(larvae>=1){
  #   for(larva in 1:larvae){
  #     home=sample(c(1:K),1)
  #     if(dynamics[i,'NBB',home]==0) dynamics[i,'NBB',home]=1
  #   }
  # }
  
  #Way 2: fill as many empty spaces as there are larvae (i.e. a larva wanders til it finds an empty space)
  empty=which(dynamics[i,'NBB',]==0)
  homes=sample(empty,min(larvae,length(empty)))
  dynamics[i,'NBB',homes]=1
}

#create matrix pop to collate rank totals at each timestep (and match stage-structured model output)
pop=data.frame('F'=rep(NA,n),'M'=rep(NA,n),'NBA'=rep(NA,n),'NBB'=rep(NA,n))
for(i in 1:n){
  pop[i,'F']=nnzero(dynamics[i,'F',1:K])
  pop[i,'M']=nnzero(dynamics[i,'M',1:K])
  pop[i,'NBA']=nnzero(dynamics[i,'NBA',1:K])
  pop[i,'NBB']=nnzero(dynamics[i,'NBB',1:K])
}


gen<-c(1:n)
totalF=pop[,1]
totalM=pop[,2]
totalJA=pop[,3]
totalJB=pop[,4]
totalJ=totalJA+totalJB

plot(gen,totalF, ylim=c(0,K), type="l", lwd=2, col=rgb(.5,0,.7),
     xlab='time (months)',ylab='class population')
title(main='Population Dynamics: no harvest')
lines(gen,totalM,type="l",col=4, lwd=2)
lines(gen,totalJA,type="l",col=2, lwd=2)
lines(gen,totalJB,type="l",col=1, lwd=2)
legend("topright",c('F','M','JA','JB'),col=c(rgb(.5,0,.7),4,2,1),lty=c(1,1,1,1),cex=1)
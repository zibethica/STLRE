#################################################################
####### Stochastic LTRE - SLTRE - LEAF HARVEST AND FIRE ####### 
###############################################################
#Georgia Hart-Fredeluces
#University of Hawaii/Idaho State University
#hartgeo2@isu.edu
#August 2020
##########################

##PARTS####
#This code has two similar and sequential parts. Part 1 is code to run STLRE that was used for the leaf harvest simulations. Part 2 is code to run STLRE that was used to compare fire severity and fire frequency.

#Attributions#####
## Original code name: Lomatium1.m
## Original MATLAB code by Hal Caswell, Woods Hole Oceanographic Institution (hcaswell@whoi.edu), dated November 2009, licensed under the Creative Commons Attribution-Noncommercial-Share Alike 3.0 United States License. 

## MATLAB code (mostly) translated into R by
## Orou Gaoue, formerly University of Hawaii at Manoa
## Februrary 14, 2015
## Current affiliation: University of Tennessee Knoxville, ogaoue@utk.edu

## R code translation continued and code adjusted for analysis for beargrass by G. Hart-Fredeluces on October 31st 2018, support from T. Ticktin and L. Bialic-Murphy fall 2018

#IPM and bigmat functions modified from code provide by Lisa Mandle, Natural Capital Project, Stanford University

#set working directory
setwd("")

#R version used in this script: 3.5.2 (2018-12-20) -- "Eggshell Igloo"

#################################################################
####### PART I - LEAF HARVEST ####### 
###############################################################

####clear all####
rm(list=ls(all=TRUE))  # Clear all

#####packages needed for this script###########
#install.packages("truncnorm")
#install.packages("msm")
library(truncnorm)
library(msm)


#### Objects needed #########
p.vec<-readRDS(file="pvec_vitalrates.rds") #loading a matrix of regression coefficients, rates and size distributions for the vital rate regression models


##### The IPM FUNCTIONS: s(x), g(y,x), f(y,x), c(y,x), p(y,x), K(y,x)#####
##### A. SURVIVAL function s(x) SEEDLINGS + NONSEEDLINGS####
#survival function for ramets
sx.ramet<-function(x,pvec,canopy,harvest,burn,year) {
  #below functon has no harvest
  xbeta<-pvec[1,2,burn,year]+(pvec[2,2,burn,year])*x+pvec[11,2,burn,year]*canopy+pvec[16,2,burn,year]*canopy*x+pvec[14,2,burn,year]*harvest*0+pvec[19,2,burn,year]*harvest*x*0;#everyone no harvest
  #below function is with harvest
  xbeta_h<-pvec[1,2,burn,year]+(pvec[2,2,burn,year])*x+pvec[11,2,burn,year]*canopy+pvec[16,2,burn,year]*canopy*x+pvec[14,2,burn,year]*harvest+pvec[19,2,burn,year]*harvest*x;#everyone with harvest
  s_nh<-(exp(xbeta))/(1+exp(xbeta)) #no flowering, no harvest
  s_h<-(exp(xbeta_h))/(1+exp(xbeta_h)) #conversion nonseedling, no flowering, with harvest 
  xbetafl<-pvec[1,3,burn,year]+pvec[2,3,burn,year]*x+pvec[11,3,burn,year]*canopy; #prob of flowering equation
  p.cap<-(exp(xbetafl)/(1+exp(xbetafl)))#converting prob of flowering
  sp_f<-s_nh*(1-(p.cap)) #multiplying survival prob by the prob of not flowering bc flowering is another way to die; flowering possible here, but not harvest
  sp_h_f<-s_h*(1-(p.cap)) #for plants above flowering and harvest size
  sp_all<-c(s_nh[x<=log(7)],sp_f[x>log(7)&x<log(8)],sp_h_f[x>=log(8)])
  #changed 17 to 8
  return(sp_all)
}

#survival function for seedlings
sx.single<-function(x,burn,year) {
  sdlgs_By1=0.05+0*x#seedling-size plants (no flowering or harvest of course)
  sdlgs_SCy1=0.1+0*x#seedling-size plants (no flowering or harvest of course)
  sdlgs_UBy1=0.56+0*x#seedling-size plants (no flowering or harvest of course)#0.56 is the mean overall observed survival rate
  sdlgs_By23=0.56+0*x#seedling year 2 or 3
  sdlgs_SCy23=0.56+0*x#seedling year 2 or 3
  sdlgs_UBy23=0.56+0*x#seedling year 2 or 3
  sp_all=c(sdlgs_By1[year==1&burn==1],sdlgs_SCy1[year==1&burn==2],sdlgs_UBy1[year==1&burn==3],#year 1 seedlings
           sdlgs_By23[year==2&burn==1],sdlgs_SCy23[year==2&burn==2],sdlgs_UBy23[year==2&burn==3],#year 2 seedlings
           sdlgs_By23[year==3&burn==1],sdlgs_SCy23[year==3&burn==2],sdlgs_UBy23[year==3&burn==3])#year 3 seedlings 
  return(sp_all)
}

#test functions:
#sx.ramet(c(log(0.6),log(3.5),log(6.9),log(7.1),log(8),log(20),log(100)),p.vec,55,0,2,2)
#sx.single(c(log(0.6),log(3.5),log(6.9),log(7.1),log(8),log(20),log(100)),2,1)

##### B. GROWTH FUNCTIONS g(y,x) #####
covarident_2017<-c(1.000, 1.2322, 1.1311) #calculated through model selection process in R script "XETE vital rates [date]"

gyx.ramet<-function(y, x, pvec,soilmois_e,burn,year) {
  if(year==3){
    mux<-((pvec[1,1,burn,year])+(pvec[2,1,burn,year])*x+(pvec[12,1,burn,year])*soilmois_e+pvec[17,1,burn,year]*soilmois_e*x);
    sigmax2<-pvec[4,1,burn,year]*exp(2*-0.3488*mux)*covarident_2017[burn]#the growth variance (pvec[4,1,burn,year], calculated previously) multiplied by the variance covariate which is the exponential covariance. -0.3488 is the value of the variance covariate, shown by running this: varExp2_17[2], then multiplied by the identity covariance variate obtained from from varExp2_17[1] that varies with plot type (burn category)
    sigmax<-sqrt(sigmax2)
    g<-dtnorm(y, mux, sigmax, lower=min(y), upper=max(y))
    return(g)}
  
  if(year==2){mux<-((pvec[1,1,burn,year])+(pvec[2,1,burn,year])*x);#
  muxp<-mux#ifelse(mux<=0,0.01,mux)
  sigmax2<-pvec[4,1,burn,year]*exp(2*pvec[5,1,burn,year]*(muxp)) #the variance, pvec[4,1,burn,year], calculated previously is multiplied by the exponential variance covariate with teh value -0.383 calculcated from varExp2_16, and stored as pvec[5,1,burn,year]
  sigmax<-sqrt(sigmax2)
  g<-dtnorm(y, muxp, sigmax, lower=min(y), upper=max(y))
  return(g)}
  if(year==1){mux<-((pvec[1,1,burn,year])+(pvec[2,1,burn,year])*x);
  muxp<-mux
  sigmax2<-pvec[4,1,burn,year]*exp(2*pvec[5,1,burn,year]*(muxp)) #the variance, pvec[4,1,burn,year], calculated previously is multiplied by the exponential variance covariate with teh value -0.383 calculcated from varExp2_16, and stored as pvec[5,1,burn,year]
  sigmax<-sqrt(sigmax2)
  g<-dtnorm(y, muxp, sigmax, lower=min(y), upper=max(y))
  return(g)}
}


gyx.single<-function(y,x) {
  mux.single=0.26487+0.7234*x #linear model estimated from seedling growth 
  sigmax2<-0.05791761 #variance of the above linear model wtih respect to our data
  sigma<-sqrt(sigmax2)
  g=dtnorm(y, mux.single, sigma,lower=min(y),upper=max(y))
  return(g)
}

##### C.The SURVIVAL-GROWTH function P(y, x) #####
pyx.ramet<-function(y,x,pvec,canopy,soilmois_e,harvest,burn,year) { 
  p<-sx.ramet(x,pvec,canopy,harvest,burn,year)*gyx.ramet(y,x,pvec,soilmois_e,burn,year)
  return(p) 
}

pyx.single<-function(y,x,burn,year) { 
  p<-sx.single(x,burn,year)*gyx.single(y,x)
  return(p) 
}

##### D. FERTILITY (capsules) function f(x,y) #####
fyx<-function(y, x, pvec,canopy,burn,year) {
  xbeta<-pvec[1,3,burn,year]+pvec[2,3,burn,year]*x+pvec[11,3,burn,year]*canopy;
  p.cap<-ifelse(x<log(7),0,(exp(xbeta)/(1+exp(xbeta)))) #conversion function for binomial
  no.cap<-pvec[1,4,burn,year]+pvec[2,4,burn,year]*x; ## normal distribution
  no.sdls.cap<-pvec[10,4,burn,year] #calculated separately for each burn treatment
  #dtnorm; lower=0.45, upper=max(y))
  #scd<-dtruncnorm(y,log(0.4),max(y),pvec[8,1,burn,year], pvec[9,1,burn,year]) #an alternative using a truncated normal distribution
  scd<-dnorm(y,pvec[8,1,burn,year], pvec[9,1,burn,year])
  f<-p.cap*no.cap*no.sdls.cap*scd
  #fn<-ifelse(f>0,f,0)
  return(f)
}

##### CLONAL FERTILITY#####
cyx<-function(y, x, pvec,canopy,harvest,burn,year) { 
  #the probability of flowering function sets the value of the flowering parameter
  xbeta_f<-pvec[1,3,burn,year]+pvec[2,3,burn,year]*x+pvec[11,3,burn,year]*canopy;
  p.cap<-ifelse(x<log(7),0,(exp(xbeta_f)/(1+exp(xbeta_f)))) 
  
  xbeta<-(pvec[1,5,burn,year]+pvec[2,5,burn,year]*x+pvec[13,5,burn,year]*p.cap+pvec[14,5,burn,year]*harvest*0+pvec[18,5,burn,year]*p.cap*x)#no harvest effect
  
  xbeta_h<-(pvec[1,5,burn,year]+pvec[2,5,burn,year]*x+pvec[13,5,burn,year]*p.cap+pvec[14,5,burn,year]*harvest+pvec[18,5,burn,year]*p.cap*x)#with harvest 

  p.veg<-(exp(xbeta))/(1+exp(xbeta))
  p.veg_h<-(exp(xbeta_h))/(1+exp(xbeta_h))
  no.veg<-(pvec[1,6,burn,year]+pvec[2,6,burn,year]*x+pvec[18,6,burn,year]*p.cap) 
  scd.v<-dtruncnorm(y,min(y),max(y),pvec[6,6,burn,year],pvec[7,6,burn,year])
  v<-p.veg*no.veg*scd.v
  v_h<-p.veg_h*no.veg*scd.v
  v_o<-p.veg*no.veg*scd.v*0
  v_sp<-c(v_o[x<log(8)],v[x>=log(8)&x<log(17)],v_h[x>=log(17)])
  return(v_sp)
}


#####Survival-Growth Kernals#####

## Ramet kernal: K(y,x)= p(y,x)+cyx(y,x) (survial, growth and clonal reproduction)
ramet.Kyx<-function(y, x, pvec,canopy,soilmois_e,harvest,burn,year) {#Flowering
  ramet.k<-pyx.ramet(y, x, pvec,canopy,soilmois_e,harvest,burn,year)+cyx(y, x, pvec,canopy,harvest,burn,year) 
  return(ramet.k) 
}#fertility kept separately


##SEEDLINGS.yx (growth plus survival)
single.Kyx<-function (y,x, burn,year) {
  single.k<-pyx.single(y,x, burn,year)
  return(single.k) 
}






############## Stochastic Life Table Response Experiment: SLTRE ##########

#modifed from the following:
## Caswell (2010) Life table response experiment analysis of a stochastic growth rate. Journal of Ecology 98(2): 324-333 

## Program for stochastic LTRE analysis for Lomatium
## example 1: two sites with identical environmental dynamics but different vital rate responses. All functions called by this program are included within this file

## Original code name: Lomatium1.m
## Original MATLAB code by Hal Caswell, Woods Hole Oceanographic Institution (hcaswell@whoi.edu), dated November 2009, licensed under the Creative Commons Attribution-Noncommercial-Share Alike 3.0 United States License. 

## MATLAB code (mostly) translated into R by
## Orou Gaoue, formerly University of Hawaii at Manoa
## Februrary 14, 2015
## Current affiliation: University of Tennessee Knoxville, ogaoue@utk.edu

## R code translation continued and code adjusted for analysis for beargrass by G. Hart-Fredeluces on October 31st 2018, support from T. Ticktin and L. Bialic-Murphy fall 2018

#IPM and bigmat functions modified from code provide by Lisa Mandle, Natural Capital Project, Stanford University

#the following code simulates populations with and without harvest, comparing stochastic growth rate and decomposing the causes of any differences in stochastic population growth rate. In the harvest simulation, harvest has a 33.33% chance of occurring each year when there is not fire. Further, in the harvest simualation, when a fire occurs, leaf harvest will automatically not occur one or two years post-fire (because our experiment did not include harvest in these years post-fire), but will occur in the third year after fire. This simulation is performed across our fire scenarios of interest (BAU, PRCUF(or INDFG) and NF). It compares the environment with and without harvest for each scenario.

###############################################################################
## PART I: CREATING DATA: Markov chain, array of population matrices
###############################################################################

#BOLTED IPM matrix####
bigmat<-function(bigM, pvec,Canopy,Soilmois_e,Harvest,burn,year){#Flowering,Soilmois_l
  ## Set matrix size and convergence tolerance 
  min.size<-log(0.9*0.5)                       
  max.size<-log(1.1*120)
  
  # Compute meshpoints iteration matrix KD 
  h=(max.size-min.size)/(bigM+1);#size/spacing of bins
  y.ramets = (seq(min.size, max.size, length=bigM) + seq(min.size+h, max.size+ h, length=bigM))/2; #midpoints 
  y.single<-y.ramets[y.ramets<log(3)] #seedlings/singletons defined to be less than 3
  
  #covariates
  #ramets only since seedling kernels do not include these
  canopy.ramet<-rnorm(length(y.ramets),mean=Canopy,sd=0)
  soilmois_e.ramet<-rnorm(n=length(y.ramets),mean=Soilmois_e, sd=0)
  harvest.ramet<-rnorm(n=length(y.ramets),mean=Harvest,sd=0)
  
  ## Apply Kyx funct for y and y for seedlings and ramets separately
  K.single=outer(y.ramets,y.single, single.Kyx, burn, year);
  KD.single=h*K.single; #matrix for singleton transitions
  
  K.ramet=outer(y.ramets,y.ramets, ramet.Kyx, pvec,canopy.ramet,soilmois_e.ramet,harvest.ramet,burn,year);
  KD.ramet=h*K.ramet
  
  F.sdlg=outer(y.single,y.ramets, fyx, pvec,canopy.ramet, burn, year);
  FD.sdlg=h*F.sdlg
  FD.sdlg[,1:sum(y.ramets<log(7))]<-0  #no seedlings made by ramets <7
  
  #Bolting -- put together seedlings, ramets and fertility 
  mat.dim<-length(y.ramets)+length(y.single)
  full.matrix<-matrix(nrow=mat.dim, ncol=mat.dim) #pre-make matrix, should be faster than using rbind/cbind
  full.matrix[1:length(y.single),1:length(y.single)]<-KD.single[1:length(y.single),] #seedlings staying seedlings from KD.single
  full.matrix[(length(y.single)+1):(2*length(y.single)),1:length(y.single)] <- matrix(0,length(y.single),length(y.single)) #seedlings can't become ramets<3 (by definition), so assign 0
  full.matrix[(2*length(y.single)+1):mat.dim,1:length(y.single)]<-KD.single[(length(y.single)+1):length(y.ramets),] #seedlings becoming ramets from KD.single
  full.matrix[1:length(y.single),(length(y.single)+1):mat.dim]<-FD.sdlg
  full.matrix[(length(y.single)+1):mat.dim,(length(y.single)+1):mat.dim]<-KD.ramet
  return(full.matrix);  ## This is the full, bolted-together matrix
}


#cell size####
cs<-6

bcs<-cs*(4/3) #adjustment of cell size for bolted matrix

#generating matrices across years since fire and fire severity and with and without harvest in the third year at mean covariate values by year and fire severity####
As<-array(0, c(bcs,bcs,3,3,2)) #change to cs for unbolted option

harvestmat<-array(c(0,0,0,0,0,0,0,0,0,0,0,0.3786,0,0,0.3152,0,0,0.2731),dim=list(3,3,2))

canmat<-c(47.65659,40.15312,21.4834)

soilmat_e<-matrix(data=c(14.04873,14.04873,13.22821,19.78298,19.78298,20.21367,16.92545,16.92545,13.05035),nrow=3,ncol=3,byrow=FALSE,dimnames = list(c("2015","2016","2017"),c("burned","scorched","unburned")) )

for (j in 1:3){#high sev, low sev, unburned
  for (k in 1:3){#2015, 2016, and 2017
    for (h in 1:2){#2 means has been harvested
      As[,,j,k,h]<-bigmat(cs, pvec=p.vec,Canopy=canmat[j],Soilmois_e=soilmat_e[k,j],Harvest=harvestmat[k,j,h],burn=j,year=k) #using means by burn class and year for covariates
    }}}

B_hs1<-As[,,1,1,1]#setting j (year) to 1 and k (burn class) to 1 (high severity)
B_hs2<-As[,,1,2,1]#setting j (year) to 2 and k (burn class) to 1 (high severity)
B_hs3<-As[,,1,3,1]#setting j (year) to 3 and k (burn class) to 1 (high severity)
B_hs3H<-As[,,1,3,2]#setting j (year) to 3 and k (burn class) 
B_ls1<-As[,,2,1,1]#setting j (year) to 1 and k (burn class) to 2 (low severity)
B_ls2<-As[,,2,2,1]#setting j (year) to 2 and k (burn class) 
B_ls3<-As[,,2,3,1]#setting j (year) to 3 and k (burn class) 
B_ls3H<-As[,,2,3,2]#setting j (year) to 3 and k (burn class) 
UB_1<-As[,,3,1,1]#setting j (year) to 1 and k (burn class) to 3 (no fire)
UB_2<-As[,,3,2,1]#setting j (year) to 2 and k (burn class) to 3 (no fire)
UB_3<-As[,,3,3,1]#setting j (year) to 3 and k (burn class) to 3 (no fire)
UB_3H<-As[,,3,3,2]#setting j (year) to 3 and k (burn class) to 3 (no fire)

#for harvest or not comparison
B_wharv<-array(c(B_hs1,B_hs2,B_hs3H,B_ls1,B_ls2,B_ls3H,UB_1,UB_2,UB_3H), c(bcs,bcs,9)) #change to cs, cs for unbolted
B_noharv<-array(c(B_hs1,B_hs2,B_hs3,B_ls1,B_ls2,B_ls3,UB_1,UB_2,UB_3), c(bcs,bcs,9))#change to cs, cs for unbolted
B_all<-array(c(B_wharv,B_noharv),c(bcs,bcs,9,2))#change to cs, cs for unbolted

#time limit for iterations####
tl<-100

###############################################################################
## PART II: Vec operator or c()
###############################################################################
vec<-function(x){
  y=c(x)
  return(y) #added this return part because the function wasn't returning y
}

###############################################################################
### PART III: Creating a function to generate envir sequence from a markov chain implemented with a times-since-fire component set by the weibull hazard function
###############################################################################
#modifying for the case of weibull (probabilities of fire change with time since fire)
#Generating a sequence of fire states####
#modifying the sequences function (created by O. Gaoue) to draw from a different probability transition matrix for each year since fire
#env=sequences_weibull(tl, P[,,,1],fscen) this is the correct input. A set of matrices for a given fire scenario, across all years.For P[,,,1]- dim are 12,12,100.
#fscen is the fire scenario, defined as 1=BAU, 2=PRCUF, or 3=No Fire

#this has been corrected for this harvest/no harvest comparison
#we have 9 possible environmental states
sequences_weibull<-function(tl,A,fscen){
  B<-array(0,c(9,9,tl))
  for(i in 1:tl){ #number of years
    B[,,i]<-apply(A[,,i],2,cumsum) # cumulative sums of each column calculated separately for each time step dimension of the array and for each fire scenario.
  }
  #manually calculating culmulative sums for columns 9-12
  B[2:3,7:9,]<-0 #removing cumsums for impossible transitions
  B[5:6,7:9,]<-0 #removing cumsums for impossible transitions
  
    #for(j in 1:tl){#for each year
     # if(fscen!=3){#if this is NOT the unburned scenario
     # B[7,7:9,j]<-((1-B[4,7,j])/3)+B[4,7,j] #recalculating cumsums for not burning 
     # B[8,7:9,j]<-2*((1-B[4,7,j])/3)+B[4,7,j]
     # B[9,7:9,j]<-3*((1-B[4,7,j])/3)+B[4,7,j]
   # }}
 
  A; B;
  states=numeric(tl+1); # A vector of "nt" 0s to receive the sequence.
  rd=numeric(tl) # nt=tl
  rd=runif(tl);  	 # Random uniform distr of norm "nt"
  states[1] = 7; 		 # Start in unburned state
  for(i in 1:tl) {
    burn_years<-ifelse(1%in%states|4%in%states,which(states==1|states==4),1)#if a burn has occurred there is a 1 or a 4 in the states, this will return the positions of that state,if there is no burn yet, it will return the number 1 
    tsf<-min(i-burn_years)+1#calculating minimum time since the last fire in the year i used to predict fire in year i+1
    f=B[,states[i],tsf]; #cumulative probabilities for current state, using the sub-array for the years since fire (tsf) calculated above, and using the state in the current year (i) to project the next year
    states[i+1]=min(grep(TRUE,rd[i]<f))#taking the smallest row value or row vector index which is larger than the random draw
  }
  states2=states[-1]#removing the state 4 at the start of the string
  
  return("env"= states2) ## The function will yield a seq of states
}

###############################################################################
## PART IV: Creating function "stochgrow": To calculate the stochastic growth rate frow a markov chain environment
###############################################################################

stochgrow<-function(env,b){
  ## routine to calculate the stochastic growth rate for a markov chain environment
  
  ## input
  #   env=sequence of environmental states (for me this will be years), and is tl 
  #   b=array of environment-specific projection matrices (s x s x number of  environments)
  
  ## output
  #   loglams=stochastic growth rate
  
  ## number of stages
  s=dim(b[,,1])[1]
  
  #alternatively, starting with stable stage distribution for UB2
  ## setting initial stage distribution based on stable stage distribution for unburned plots in year 2
#  egA1<-eigen(bigmat(cs, pvec=p.vec,Canopy=canmat[3],Soilmois_e=soilmat_e[2,3],Soilmois_l=soilmat_l[2,3],Flowering=flwrmat[2,3],Harvest=harvestmat[2,3],burn=3,year=2)) #using stable stage distribution for the unburned plots in year 2
#  meig<-which(Re(egA1$values)==max(Re(egA1$values)))#position of dominant eigenvalue (though I think should always be 1)
#  W<-egA1$vectors
#  w<-abs(Re(W[,meig])) ## SSD: will be used as initial vectors n
#  n<-w /sum(w) #fir
  
  #starting with simple stage distribution
  #nzero=matrix(1,s,1)/s  
  #n=nzero   ## a s x 1 vector 
  #this section was changed from Orou's version to match Lisa's code
  nzero=matrix(1,s,1)/s  
  n0=array(0,dim=c(length(nzero),tl+1))
  n0[,1]<-nzero  
  N=array(0,c(tl+1))
  N[1]<-sum(n0[,1])
  r=array(0,c(tl+1))
  r[1]<-1
  
  for (t in 1:tl){
    n0[,t+1]=b[,,env[t]]%*%n0[,t]    
    N[t+1]=sum(n0[,t+1])
    r[t+1]=log(N[t+1]/N[t]) 
  } 
  
  loglams=exp(mean(r[2:tl+1]));

  #adding confidence intervals
  loglams_SE=sqrt(var(r[2:(tl+1)]/tl)) #standard error of r

  return(list(loglams,loglams_SE))
}
  

###############################################################################
## PART V: Creating function "stocsens_e": To calculates environment-specific sensitivity of stochastic growth rate
###############################################################################

#function to calculate sensitivity of a given environment in a given time step
#frame 8
sens_t<-function(t, wout, vout, dvecBdt, env, r){
  sens=array(0, dim=c(1, bcs^2)) #change to cs for unbolted
  sens=kronecker(t(wout[,t]),t(vout[,t+1]))%*%dvecBdt[,,env[t]]/(as.numeric(r[t]*t(vout[,t+1])%*%wout[,t+1]))
  return(sens)
}

#function to determine the mean sensitivity of an environment across the environmental sequence
#frame 5
stoch_sens<-function(ienv, env, wout, vout, dvecBdt, r){
  #determine years that match the given environment
  indices<-which(env==ienv)#ifelse(ienv%in%env,which(env==ienv),NA)
  #create an array to store the series of sensitivities across indices
  sens_save=array(0, dim=c(1, bcs^2, length(indices))) #change to cs for ubolted
  #save a simple array to create the list to hold lapply output
  sens_simple=array(0, dim=c(1, bcs^2)) #change to cs for ubolted
  #use this simple array to create the empty list to hold output
  sens_save_list<-replicate(length(indices), list(sens_simple))
  #apply the sens_subfunc to each index if indices exist
  sens_save_list[length(indices)>0]=lapply(indices,FUN=function(t){sens_t(t, wout, vout, dvecBdt,env,r)})
  #convert the outcome list back to an array
  sens_save[length(indices)>0]=array(as.numeric(unlist(sens_save_list)), dim=c(1, bcs^2, length(indices)))
  #create an array to store the sum of sensitivities across indices
  sens_save_ov<-array(0, dim=c(1, bcs^2))
  #take the sum of sensitivities across indices
  sens_save_ov[length(indices)>0]=apply(sens_save,1:2,sum)#does this work?
  #set up an array to hold final output for this environment
  dlogldt.p<-array(0,dim=c(1,bcs^2))
  #populate the array with the mean sensitivity of each element for the environment divided by the number of time steps in the simulation
  dlogldt.p=sens_save_ov/tl 
  #return this value from the function
  return(dlogldt.p)
}

#function to calculate the sensitivity of all environments across the environmental sequence
#frame 2
stocsens_e<-function(P,B,dvecBdt,tl,fscen){
  ## calculates environment-specific sensitivity of stochastic growth rate
  
  ## input
  #   P=markov transition matrices (ne x ne x tl)
  #   B=array of projection matrices (s x s x ne)
  #   dvecBdt=array of derivatives of the B_i to parameters (s^2 x p x ne)
  #   tl=time limit for iterations
  
  ## output
  #   dlogldt=array of derivatives of log lambda_s to parameters (1 x p x ne)
  
  ## number of stages, number of environments
  s=dim(B_all)[1] #number of stages is life stages
  ne=dim(B_all)[3] #number of environments (is the 12 states we can chose from)B1,B2,etc.)
  
  ## squeeze P if necessary
  #P=drop(P)  ## changed "squeeze" (matlab) to "drop" (R)
  
  ## generate environmental sequence from P
  env=sequences_weibull(tl,P,fscen) #changed P to P[,,,fscen]
  
  ## calculate stochastic growth rate
  cat('loglams')
  
  #this section was changed to match Lisa's code
  nzero=matrix(1,s,1)/s  
  n0=array(0,dim=c(length(nzero),tl+1))
  n0[,1]<-nzero  
  N=array(0,c(tl+1))
  N[1]<-sum(n0[,1])
  r=array(0,c(tl+1))
  r[1]<-1
  
  for (t in 1:tl){
    n0[,t+1]=B[,,env[t]]%*%n0[,t]    #B<-B_mean for error checking
    N[t+1]=sum(n0[,t+1])
    r[t+1]=log(N[t+1]/N[t]) 
  } 
  
  loglams=mean(r[2:tl+1]);
  
  ## calculate sensitivity, using forward and backward iteration
  ## see MPM, Section 14.4.1.1
  cat(' forward iteration') 
  
  ## forward iteration for stage distribution
  wout=array(0,c(s,tl+1))
  w=array(1,c(s,1))/s #create an array of vectors of equal frequencies
  wout[,1]=w #fill w into wout in the first column position
  
  for (it in 1:tl){ #iterate forward
    w=B[,,env[it]]%*%w #right side matrix multiplication of w and the appropriate matrix from the set of mean matrices (B)
    r[it]=sum(w) #save the population size to r
    w=w/sum(w) #make w proportional
    wout[,it+1]=w #store w in the wout object
  } 
  
  cat(' backward iteration') #### 
  
  ## backward iteration for repro value
  vout=array(0,c(s,tl+1))  
  v=array(1,c(1,s))/s #create a matrix of equal initial frequencies
  vout[,tl+1]=t(v) #transpose the vector and place it as the final column in vout
  for (it in tl:1){ #for the backwards iteration
    v=v%*%B[,,env[it]] #multiple v by the appropriate matrix in the sequence (env), this is left multiplication
    v=v/sum(v); #redefine v from above as a proportion
    vout[,it]=t(v) #fill the transpose of this vector into vout
  } 
  
  ## now calculate the environment-specific sensitivities; eqn 16
  cat('env-specific sensitivities') 
  
  #dlogldt.p<-array(0,dim=c(1,cs^2))
  #dlogldt<-array(0,dim=c(1,cs^2,ne))
  
  #returning to lappy for HPC use
  dlogdt_result<-lapply(1:ne,FUN=function(ienv){stoch_sens(ienv,env,wout, vout,dvecBdt,r)})
  
  #option for parallelizing, Macintosh only
  #dlogdt_result<-mclapply(1:ne,FUN=function(ienv){stoch_sens(ienv,env)},mc.cores=4,mc.preschedule = TRUE)
  
  #converting the list into an array (could have just used sapply)
  dlogldt<-array(as.numeric(unlist(dlogdt_result)), dim=c(1, bcs^2, ne))
  
  return(dlogldt)
} ## for ienv


################################################################################
## PART VI: Creating function "beargrass_1" to do SLTRE
#################################################################################

beargrass_1<-function(B_wharv, B_noharv){
  ## environment-specific matrices, averaged between treatments
  B_mean=(B_wharv + B_noharv)/2  ## this is an n by n by l array that creates average matrices between B_fb and B_rp
  
  ## Time limit for iterations
  tl=100
  
  ## number of stages (s), environments (ne), and treatments (nt)
  ## 100 stages, 3 envir states (years) and 3 treatments (HS and LS)
  s=dim(B_mean)[1] 
  ne=dim(B_mean)[3]
  #nt=dim(B_all)[4]
  
  #arrays to hold outputs
  nf=3
  loglams_wharv<-array(0,c(1,nf))
  loglams_noharv<-array(0,c(1,nf))
  loglams_SE_wharv<-array(0,c(1,nf))
  loglams_SE_noharv<-array(0,c(1,nf))
  contmats<-array(0,c(bcs^2,ne,nf)) #change to cs for unbolted
  contmats_env<-array(0,c(ne,nf))
  contmats_aij<-array(0,c(bcs,bcs,nf)) #change to cs for unbolted
  #dlogldt_save<-array(0,dim=c(1,cs^2,ne,nf))
  
  ## in all calculations
  ## envt 1 = with harvest (harvest occurs the third year following fire, and every third year for the unburend scenario)
  ## envt 2 = without harvest (unburned third year kernal is always without harvest)
  
  ## derivatives of the matrices B wrt the parameters in the paper, the parameters are the elements of B
  dvecBdt<-array(0, c(s^2, s^2, ne))
  for(i in 1:ne){
    dvecBdt[,,i]<-diag(s^2);
  } # for i
  
  #P here is created for the specific harvest/no harvest comparisons
  #create P####
  #setting of probability of fire matrices
  #hazard function
  FRP <- function(tl,wb,wc){ 
    fp<-((wc*tl)^(wc-1))/(wb^wc)
    return(fp)
  }
  wb<-180#180 #fire recurrence interval that will be exceeded 36.8% of the time for weibull function
  wc<-1.5 #dependence of fire prob on fuel accumulation, c=1 means that we are actually using a negative exponential distribution, not a weibull and that there is no fuel-dependency 
  
  #creating an array of zeros of appropriate dimensions
  P<-array(0,c(9,9,tl,nf),dimnames=list(c("B1","B2","B3","L1","L2","L3","UB1","UB2","UB3"),c("B1","B2","B3","L1","L2","L3","UB1","UB2","UB3"),1:tl,c("BAU","PRCUF","NF"))) 
  #setting unburned matrix probabilities
  ub1_p<-0.17
  ub2_p<-0.415#1-(0.10+0.333)
  ub3_p<-0.415#0.333
    
  
  #filling in fire probabilities across the array from weibull equation (BAU scenario)
  for(i in 1:tl){
    burnprob<-FRP(i,wb,wc)
    P[1,7:9,i,1]<-burnprob*0.58#58% chance high severity whenit burns
    P[4,7:9,i,1]<-burnprob*0.42#42% chance low severity whenit burns
  }
  #this will calculate the cumulative probability of fire over time, with a given fire-return interval(wb) and fuel dependency (wc),and fill in the cumulative fire probability across each sub array by time step. It will do so for the first fire scenario, "BAU."
  
  #fill in the corresponding probability that it will NOT burn
  for(i in 1:tl){ #for each time step
    P[7,7:9,i,1]<-(1-(P[1,7,i,1]+P[4,7,i,1]))*ub1_p
    P[8,7:9,i,1]<-(1-(P[1,7,i,1]+P[4,7,i,1]))*ub2_p
    P[9,7:9,i,1]<-(1-(P[1,7,i,1]+P[4,7,i,1]))*ub3_p
  }#the probability of not burning in each year is given by: 1 - (the prob of burning). I multiple by 0.10 for UB1 and 0.45 for UB2 and UB3 unburned states having equal probability
  
  #filling in fire probabilities across the array for probability of burning hs (0.1) and low sev (0.9) (PRCUF scenario) every 10 years
  #HS
  P[1,7,,2]<-rep_len(c(0,0,0,0,0,0,0,0,0,0.1),length.out=tl)
  P[1,8,,2]<-rep_len(c(0,0,0,0,0,0,0,0,0,0.1),length.out=tl)
  P[1,9,,2]<-rep_len(c(0,0,0,0,0,0,0,0,0,0.1),length.out=tl)
  #LS
  P[4,7,,2]<-rep_len(c(0,0,0,0,0,0,0,0,0,0.9),length.out=tl)
  P[4,8,,2]<-rep_len(c(0,0,0,0,0,0,0,0,0,0.9),length.out=tl)
  P[4,9,,2]<-rep_len(c(0,0,0,0,0,0,0,0,0,0.9),length.out=tl)
  
  #fill in the corresponding probability that it will NOT burn (PRCUF scenario)
  for(i in 1:tl){#for each time step
    P[7,7:9,i,2]<-(1-(P[4,7,i,2]+P[1,7,i,2]))*ub1_p
    P[8,7:9,i,2]<-(1-(P[4,7,i,2]+P[1,7,i,2]))*ub2_p
    P[9,7:9,i,2]<-(1-(P[4,7,i,2]+P[1,7,i,2]))*ub3_p
  }#the probability of not burning in each year is given by: 1 - (the prob of burning). I multiple by 0.10 for UB1 and 0.45 for UB2 and UB3 unburned states having equal probability
  #No fire (NF scenario) - ~45% chance of harvest in any year
  for(i in 1:tl){#for each time step
    P[7,7:9,i,3]<-ub1_p
    P[8,7:9,i,3]<-ub2_p
    P[9,7:9,i,3]<-ub3_p 
  }
  
  #adding guaranteed transitions from B1 to B2 to B3 to UB1,2 or 3
  #high sev
  for(i in 1:tl){#for each time step
    P[2,1,i,1:2]<-1 #the prob of B1 to B2
  }
  for(i in 1:tl){ #for each time step
    P[3,2,i,1:2]<-1 #the prob of B2 to B3, or B3H for harvested scenario
  }
  for(i in 1:tl){ #for each time step
    P[7,3,i,1:2]<-ub1_p
    P[8,3,i,1:2]<-ub2_p #automatic return to unburned state, set probs
    P[9,3,i,1:2]<-ub3_p
  }
  #low sev
  for(i in 1:tl){#for each time step
    P[5,4,i,1:2]<-1 #the prob of L1 to L2
  }
  for(i in 1:tl){ #for each time step
    P[6,5,i,1:2]<-1 #the prob of L2 to L3 or L3H
  }
    for(i in 1:tl){ #for each time step
    P[7,6,i,1:2]<-ub1_p
    P[8,6,i,1:2]<-ub2_p#automatic return to unburned state
    P[9,6,i,1:2]<-ub3_p
  }
  
  #end P####
  
  # evaluate stochastic growth rates
  env_save<-array(0,c(tl,nf))
  dlogldt<-array(0,dim=c(1,bcs^2,ne)) #change to cs for unbolted
  
  for(i in 1:nf){ 
    print(c("In progress:","firescen:",i, "time:", date()), quote=FALSE)
    # generate environmental state sequences
    env<-sequences_weibull(tl, P[,,,i],i) ## new function
    
    env_save[,i]<-env
    # growth rate with harvest for each fire scenario
    loglams_wharv[i]<-stochgrow(env,B_wharv)[[1]]
    
    # growth rate without harvest for each fire scenario
    loglams_noharv[i]<-stochgrow(env,B_noharv)[[1]]
  
    # environment-specific sensitivities to the vital rates
    dlogldt<-stocsens_e(P[,,,i],B_mean,dvecBdt,tl,i)
    # dlogldt dimensions 1 x number of parametrs x number of envts
    
    # save the sensitivities for all fire freqs
    #dlogldt_save[,,,i]<-dlogldt 
    
    # contributions from vital rate differences
    for (ienv in 1:ne){ # contmats dim 1 x no parameters x no envts x no. of freq
    contmats[,ienv,i]<-t(vec(B_wharv[,,ienv]- B_noharv[,,ienv]))*dlogldt[,,ienv]              #the vec function vectorizes the matrix, column by column. After being transposed, its dimentions are 1 X cs^2           
   } # for ienv
    
    # sum contributions over matrix elements - for each i, which is fire frequency
    #contmats_env=array(0,c(ne,length(freqarray)))
  contmats_env[,i]<-apply(contmats[,,i],2,sum) 
    
    #sum contributions over envts
  contmats_aij[,,i]<-matrix(apply(contmats[,,i],1,sum),bcs,bcs) #change to cs for unbolted
  }
  return(list(loglams_wharv,loglams_noharv,contmats, contmats_env, contmats_aij,env_save))
}    

## This function returns:
## loglams_rp (growth rate for Hs at each fire freq)
## loglams_fb (growth rate for Ls at each fire freq)
## contmats (contributions from vital rate differences)
## contmats_env (sum contributions over matrix elements) 
## contmats_aij (sum contributions over envts)
#} #for i this completes the loop starting at P

## with and without harvest
SLTRE_list_harvest<-beargrass_1(B_wharv, B_noharv) 

names(SLTRE_list_harvest)<-c("loglams_wharv",
                     "loglams_noharv",
                     "contmats",
                     "contmats_env",
                     "contmats_aij",
                     "env_save"
                     )

#save results 
saveRDS(object=SLTRE_list_harvest, file="SLTRE_harvest_output.rds")

#################################################################
####### PART II - FIRE ####### 
###############################################################

####clear all####
rm(list=ls(all=TRUE))  # Clear all

#####IPM functions needed for this script###########
#install.packages("truncnorm")
#install.packages("msm")
library(truncnorm)
library(msm)

#### Objects needed #########
p.vec<-readRDS(file="pvec_vitalrates.rds")#martrix of regression coefficients, rates and size distributions for vital rate regression functions

##### The IPM FUNCTIONS: s(x), g(y,x), f(y,x), c(y,x), p(y,x), K(y,x)#####
##### A. SURVIVAL function s(x) SEEDLINGS + NONSEEDLINGS####
#survival function for ramets
sx.ramet<-function(x,pvec,canopy,harvest,burn,year) {
  #below functon has no harvest
  xbeta<-pvec[1,2,burn,year]+(pvec[2,2,burn,year])*x+pvec[11,2,burn,year]*canopy+pvec[16,2,burn,year]*canopy*x+pvec[14,2,burn,year]*harvest*0+pvec[19,2,burn,year]*harvest*x*0;#everyone no harvest
  #below function is with harvest
  xbeta_h<-pvec[1,2,burn,year]+(pvec[2,2,burn,year])*x+pvec[11,2,burn,year]*canopy+pvec[16,2,burn,year]*canopy*x+pvec[14,2,burn,year]*harvest+pvec[19,2,burn,year]*harvest*x;#everyone with harvest
  s_nh<-(exp(xbeta))/(1+exp(xbeta)) #no flowering, no harvest
  s_h<-(exp(xbeta_h))/(1+exp(xbeta_h)) #conversion nonseedling, no flowering, with harvest 
  xbetafl<-pvec[1,3,burn,year]+pvec[2,3,burn,year]*x+pvec[11,3,burn,year]*canopy; #prob of flowering equation
  p.cap<-(exp(xbetafl)/(1+exp(xbetafl)))#converting prob of flowering
  sp_f<-s_nh*(1-(p.cap)) #multiplying survival prob by the prob of not flowering bc flowering is another way to die; flowering possible here, but not harvest
  sp_h_f<-s_h*(1-(p.cap)) #for plants above flowering and harvest size
  sp_all<-c(s_nh[x<=log(7)],sp_f[x>log(7)&x<log(8)],sp_h_f[x>=log(8)])
  #changed 17 to 8
  return(sp_all)
}

#survival function for seedlings
sx.single<-function(x,burn,year) {
  sdlgs_By1=0.05+0*x#seedling-size plants (no flowering or harvest of course)
  sdlgs_SCy1=0.1+0*x#seedling-size plants (no flowering or harvest of course)
  sdlgs_UBy1=0.56+0*x#seedling-size plants (no flowering or harvest of course)#0.56 is the mean overall observed survival rate
  sdlgs_By23=0.56+0*x#seedling year 2 or 3
  sdlgs_SCy23=0.56+0*x#seedling year 2 or 3
  sdlgs_UBy23=0.56+0*x#seedling year 2 or 3
  sp_all=c(sdlgs_By1[year==1&burn==1],sdlgs_SCy1[year==1&burn==2],sdlgs_UBy1[year==1&burn==3],#year 1 seedlings
           sdlgs_By23[year==2&burn==1],sdlgs_SCy23[year==2&burn==2],sdlgs_UBy23[year==2&burn==3],#year 2 seedlings
           sdlgs_By23[year==3&burn==1],sdlgs_SCy23[year==3&burn==2],sdlgs_UBy23[year==3&burn==3])#year 3 seedlings 
  return(sp_all)
}

##### B. GROWTH FUNCTIONS g(y,x) #####
covarident_2017<-c(1.000, 1.2322, 1.1311) #calculated through model selection process in "XETE vital rates [date]" file

gyx.ramet<-function(y, x, pvec,soilmois_e,burn,year) {
  if(year==3){
    mux<-((pvec[1,1,burn,year])+(pvec[2,1,burn,year])*x+(pvec[12,1,burn,year])*soilmois_e+pvec[17,1,burn,year]*soilmois_e*x);
    sigmax2<-pvec[4,1,burn,year]*exp(2*-0.3488*mux)*covarident_2017[burn]#the growth variance (pvec[4,1,burn,year], calculated previously) multiplied by the variance covariate which is the exponential covariance. -0.3488 is the value of the variance covariate, shown by running this: varExp2_17[2], then multiplied by the identity covariance variate obtained from from varExp2_17[1] that varies with plot type (burn category)
    sigmax<-sqrt(sigmax2)
    g<-dtnorm(y, mux, sigmax, lower=min(y), upper=max(y))
    return(g)}
  
  if(year==2){mux<-((pvec[1,1,burn,year])+(pvec[2,1,burn,year])*x);#
  muxp<-mux#ifelse(mux<=0,0.01,mux)
  sigmax2<-pvec[4,1,burn,year]*exp(2*pvec[5,1,burn,year]*(muxp)) #the variance, pvec[4,1,burn,year], calculated previously is multiplied by the exponential variance covariate with teh value -0.383 calculcated from varExp2_16, and stored as pvec[5,1,burn,year]
  sigmax<-sqrt(sigmax2)
  g<-dtnorm(y, muxp, sigmax, lower=min(y), upper=max(y))
  return(g)}
  if(year==1){mux<-((pvec[1,1,burn,year])+(pvec[2,1,burn,year])*x);
  muxp<-mux#ifelse(mux<=0,0.01,mux)
  sigmax2<-pvec[4,1,burn,year]*exp(2*pvec[5,1,burn,year]*(muxp)) #the variance, pvec[4,1,burn,year], calculated previously is multiplied by the exponential variance covariate with teh value -0.383 calculcated from varExp2_16, and stored as pvec[5,1,burn,year]
  sigmax<-sqrt(sigmax2)
  g<-dtnorm(y, muxp, sigmax, lower=min(y), upper=max(y))
  return(g)}
}


gyx.single<-function(y,x) {
  mux.single=0.26487+0.7234*x#linear model estimated from seedling growth 
  sigmax2<-0.05791761 #variance of the above linear model wtih respect to our data
  sigma<-sqrt(sigmax2)
  g=dtnorm(y, mux.single, sigma,lower=min(y),upper=max(y))
  return(g)
}

##### C.The SURVIVAL-GROWTH function P(y, x) #####
pyx.ramet<-function(y,x,pvec,canopy,soilmois_e,harvest,burn,year) { 
  p<-sx.ramet(x,pvec,canopy,harvest,burn,year)*gyx.ramet(y,x,pvec,soilmois_e,burn,year)
  return(p) 
}

pyx.single<-function(y,x,burn,year) { 
  p<-sx.single(x,burn,year)*gyx.single(y,x)
  return(p) 
}

##### D. FERTILITY (capsules) function f(x,y) #####
fyx<-function(y, x, pvec,canopy,burn,year) {
  xbeta<-pvec[1,3,burn,year]+pvec[2,3,burn,year]*x+pvec[11,3,burn,year]*canopy;
  p.cap<-ifelse(x<log(7),0,(exp(xbeta)/(1+exp(xbeta)))) #conversion function for binomial
  no.cap<-pvec[1,4,burn,year]+pvec[2,4,burn,year]*x; ## normal distribution
  #sigmax2<-??
  no.sdls.cap<-pvec[10,4,burn,year] #calculated separately for each burn treatment
  scd<-dnorm(y,pvec[8,1,burn,year], pvec[9,1,burn,year])
  f<-p.cap*no.cap*no.sdls.cap*scd
  #fn<-ifelse(f>0,f,0)
  return(f)
}

##### CLONAL FERTILITY#####
cyx<-function(y, x, pvec,canopy,harvest,burn,year) { 
  #the probability of flowering function sets the value of the flowering parameter
  xbeta_f<-pvec[1,3,burn,year]+pvec[2,3,burn,year]*x+pvec[11,3,burn,year]*canopy;
  p.cap<-ifelse(x<log(7),0,(exp(xbeta_f)/(1+exp(xbeta_f)))) 
  
  xbeta<-(pvec[1,5,burn,year]+pvec[2,5,burn,year]*x+pvec[13,5,burn,year]*p.cap+pvec[14,5,burn,year]*harvest*0+pvec[18,5,burn,year]*p.cap*x)#no harvest effect
  
  xbeta_h<-(pvec[1,5,burn,year]+pvec[2,5,burn,year]*x+pvec[13,5,burn,year]*p.cap+pvec[14,5,burn,year]*harvest+pvec[18,5,burn,year]*p.cap*x)#with harvest 
  
  p.veg<-(exp(xbeta))/(1+exp(xbeta))
  p.veg_h<-(exp(xbeta_h))/(1+exp(xbeta_h))
  no.veg<-(pvec[1,6,burn,year]+pvec[2,6,burn,year]*x+pvec[18,6,burn,year]*p.cap) 
  scd.v<-dtruncnorm(y,min(y),max(y),pvec[6,6,burn,year],pvec[7,6,burn,year])
  v<-p.veg*no.veg*scd.v
  v_h<-p.veg_h*no.veg*scd.v
  v_o<-p.veg*no.veg*scd.v*0
  v_sp<-c(v_o[x<log(8)],v[x>=log(8)&x<log(17)],v_h[x>=log(17)])
  return(v_sp)
}


#####Survival-Growth Kernals#####

## Ramet kernal: K(y,x)= p(y,x)+cyx(y,x) (survial, growth and clonal reproduction)
ramet.Kyx<-function(y, x, pvec,canopy,soilmois_e,harvest,burn,year) {#Flowering
  ramet.k<-pyx.ramet(y, x, pvec,canopy,soilmois_e,harvest,burn,year)+cyx(y, x, pvec,canopy,harvest,burn,year) 
  return(ramet.k) 
}#fertility kept separately


##SEEDLINGS.yx (growth plus survival)
single.Kyx<-function (y,x, burn,year) {
  single.k<-pyx.single(y,x, burn,year)
  return(single.k) 
}





#################################################################
####### Stochastic LTRE - SLTRE ####### 
###############################################################
############## Stochastic Life Table Response Experiment: SLTRE ##########

#modifed from the following:
## Caswell (2010) Life table response experiment analysis of a stochastic growth rate. Journal of Ecology 98(2): 324-333 

## Program for stochastic LTRE analysis for Lomatium
## example 1: two sites with identical environmental dynamics but different vital rate responses. All functions called by this program are included within this file

## Original code name: Lomatium1.m
## Original MATLAB code by Hal Caswell, Woods Hole Oceanographic Institution (hcaswell@whoi.edu), dated November 2009, licensed under the Creative Commons Attribution-Noncommercial-Share Alike 3.0 United States License. 

## MATLAB code translated into R by
## Orou Gaoue, University of Hawaii at Manoa 
##currrent affiliation, University of Tennessee Knoxville, ogaoue@utk.edu
## Februrary 14, 2015

## R code adjusted for analysis for beargrass by G. Fredeluces Hart on October 31st 2018, support from T. Ticktin and L. Bialik-Murphy fall 2018

##The following code is modified to explore general responses of beargrass to combinations of fire frequency and severity. Fire severities are considered as the separate environments (these were sites in Caswell's example). Matrices for each fire severity are derived from plant measurements of separate populations in low and high severity fire areas (see Hart and Ticktin 2019 for details). Each matrix was constructed from an IPM kernal that was built for a given year for three populations across different sites. Effects of harvest are not included in this simulation.

###############################################################################
## PART I: CREATING DATA: Markov chain, array of population matrices
###############################################################################

#BOLTED IPM matrix####
bigmat<-function(bigM, pvec,Canopy,Soilmois_e,Harvest,burn,year){#Flowering,Soilmois_l
  ## Set matrix size and convergence tolerance 
  min.size<-log(0.9*0.5)                       
  max.size<-log(1.1*120)
  
  # Compute meshpoints iteration matrix KD 
  h=(max.size-min.size)/(bigM+1);#size/spacing of bins
  y.ramets = (seq(min.size, max.size, length=bigM) + seq(min.size+h, max.size+ h, length=bigM))/2; #midpoints 
  y.single<-y.ramets[y.ramets<log(3)] #seedlings/singletons defined to be less than 3
  
  #covariates
  #ramets only since seedling kernels do not include these
  canopy.ramet<-rnorm(length(y.ramets),mean=Canopy,sd=0)
  soilmois_e.ramet<-rnorm(n=length(y.ramets),mean=Soilmois_e, sd=0)
  harvest.ramet<-rnorm(n=length(y.ramets),mean=Harvest,sd=0)
  
  ## Apply Kyx funct for y and y for seedlings and ramets separately
  K.single=outer(y.ramets,y.single, single.Kyx, burn, year);
  KD.single=h*K.single; #matrix for singleton transitions
  
  K.ramet=outer(y.ramets,y.ramets, ramet.Kyx, pvec,canopy.ramet,soilmois_e.ramet,harvest.ramet,burn,year);
  KD.ramet=h*K.ramet
  
  F.sdlg=outer(y.single,y.ramets, fyx, pvec,canopy.ramet, burn, year);
  FD.sdlg=h*F.sdlg
  FD.sdlg[,1:sum(y.ramets<log(7))]<-0  #no seedlings made by ramets <7
  
  #Bolting -- put together seedlings, ramets and fertility 
  mat.dim<-length(y.ramets)+length(y.single)
  full.matrix<-matrix(nrow=mat.dim, ncol=mat.dim) #pre-make matrix, should be faster than using rbind/cbind
  full.matrix[1:length(y.single),1:length(y.single)]<-KD.single[1:length(y.single),] #seedlings staying seedlings from KD.single
  full.matrix[(length(y.single)+1):(2*length(y.single)),1:length(y.single)] <- matrix(0,length(y.single),length(y.single)) #seedlings can't become ramets<3 (by definition), so assign 0
  full.matrix[(2*length(y.single)+1):mat.dim,1:length(y.single)]<-KD.single[(length(y.single)+1):length(y.ramets),] #seedlings becoming ramets from KD.single
  full.matrix[1:length(y.single),(length(y.single)+1):mat.dim]<-FD.sdlg
  full.matrix[(length(y.single)+1):mat.dim,(length(y.single)+1):mat.dim]<-KD.ramet
  return(full.matrix);  ## This is the full, bolted-together matrix
}



#setting up cell size####
cs<-6 #cell size
bms<-cs+(cs/3) #bolted matrix size

#generating matrices across years since fire and fire severity. These are created using mean covariate values by year and fire severity####
#this part of the code depends on running the regressions and IPM in the file "XETE IPM [date].R", which 
As<-array(0, c(bms,bms,3,3))

harvestmat<-array(c(0,0,0,0,0,0,0,0,0,0,0,0.3786,0,0,0.3152,0,0,0.2731),dim=list(3,3,2))

canmat<-c(47.65659,40.15312,21.4834)

soilmat_e<-matrix(data=c(14.04873,14.04873,13.22821,19.78298,19.78298,20.21367,16.92545,16.92545,13.05035),nrow=3,ncol=3,byrow=FALSE,dimnames = list(c("2015","2016","2017"),c("burned","scorched","unburned")) )


for (j in 1:3){ #burn severity
  for (k in 1:3){ #year
    for (h in 1:2){#2 means has been harvested
      As[,,j,k]<-bigmat(cs, pvec=p.vec,Canopy=canmat[j],Soilmois_e=soilmat_e[k,j],Harvest=harvestmat[k,j,h],burn=j,year=k) #using means by burn class and year for covariates
    }}}

B_hs1<-As[,,1,1]#setting j (burn) to 1 (high severity) and k (year) to 1 
B_hs2<-As[,,1,2]#setting j (burn) to 1 (high severity) and k (year) to 2 
B_hs3<-As[,,1,3]#setting j (burn) to 1  (high severity)and k (year) to 3 
UB_1<-As[,,3,1]#setting j (burn) to 3 (unburned) and k (year) to 1
UB_2<-As[,,3,2]#setting j (burn) to 3 and k (year) to 2
UB_3<-As[,,3,3]#setting j (burn) to 3 and k (year) to 3
B_hs=array(c(B_hs1,B_hs2,B_hs3,UB_1,UB_2,UB_3), c(bms,bms,6))

B_ls1<-As[,,2,1]#setting j (burn) to 2 (low severity) and k  1 (year)
B_ls2<-As[,,2,2]#
B_ls3<-As[,,2,3]#
UB_1<-As[,,3,1]#setting j (burn) to 3 (unburned) and k (year) to 1 
UB_2<-As[,,3,2]#
UB_3<-As[,,3,3]#
B_ls=array(c(B_ls1,B_ls2,B_ls3,UB_1,UB_2,UB_3), c(bms,bms,6))

## create an array with the matrices for high severity, low severity, and all environmental states (years)
B<-array(c(B_hs, B_ls), c(bms,bms,6,2))

#time limit for iterations####
tl<-100

###############################################################################
## PART II: Vec operator or c()
###############################################################################
vec<-function(x){
  y=c(x)
  return(y) #added this return part because the function wasn't returning y
}

###############################################################################
### PART III: Creating a function to generate envir sequence from a markov chain
###############################################################################

## Function "marksim" from Caswell code is missing. O. Gaoue developed this new function "sequences" to replace it.

sequences<-function(nt,A){			   
  B = apply(A,2,cumsum); # cumulative sums of each column
  A; B;
  
  states=numeric(nt+1); # A vector of "nt" 0 to receive the sequence
  rd=runif(nt);		 # Random uniform distr of norm "nt"
  states[1] = 4; 		 # Start in open state
  
  for(i in 1:nt) {
    b=B[,states[i]]; #cumulative probabilities for current state
    states[i+1]=sum(rd[i]>b)+1 # based on current state
  }
  return("env"=states) ## The function will yield a seq of states
}


###############################################################################
## PART IV: Creating function "stochgrow": To calculate the stochastic growth rate frow a markov chain environment
###############################################################################

stochgrow<-function(env,b){
  ## routine to calculate the stochastic growth rate for a markov chain environment
  
  ## input
  #   env=sequence of environmental states 
  #   b=array of environment-specific projection matrices (s x s x number of  environments)
  
  ## output
  #   loglams=stochastic growth rate
  
  ## number of stages
  s=dim(b[,,1])[1]
  
  #this section was changed to match Lisa Mandle's code
  nzero=matrix(1,s,1)/s  
  n0=array(0,dim=c(length(nzero),tl+1))
  n0[,1]<-nzero  
  N=array(0,c(tl+1))
  N[1]<-sum(n0[,1])
  r=array(0,c(tl+1))
  r[1]<-1
  
  for (t in 1:tl){
    n0[,t+1]=b[,,env[t]]%*%n0[,t]    
    N[t+1]=sum(n0[,t+1])
    r[t+1]=log(N[t+1]/N[t]) 
  } 
  
  loglams=exp(mean(r[2:tl+1]));
  return(list(loglams,r))
}


###############################################################################
## PART V: Creating function "stocsens_e": To calculates environment-specific sensitivity of stochastic growth rate
###############################################################################

stocsens_e<-function(P,B,dvecBdt,tl){
  ## calculates environment-specific sensitivity of stochastic growth rate
  
  ## input
  #   P=markov transition matrices (ne x ne x tl )
  #   B=array of projection matrices (s x s x ne)
  #   dvecBdt=array of derivatives of the B_i to parameters (s^2 x p x ne)
  #   tl=time limit for iterations
  
  ## output
  #   dlogldt=array of derivatives of log lambda_s to parameters (1 x p x ne)
  
  ## number of stages, number of environments
  s=dim(B)[1] #number of stages is life stages
  ne=dim(B)[3] #number of environments (is the 6 years we can chose from)B1,B2,etc.)
  
  ## squeeze P if necessary
  #P=drop(P)  ## changed "squeeze" (matlab) to "drop" (R)
  
  
  dlogldt<-array(0,dim=c(1,bms^2,ne))
  #dlogldt_save<-array(0,dim=c(length(freqarray),bms^2,ne))
  
  ## generate environmental sequence from P
  #env=sequences(tl, P) #this is for time invariant probabilities
  env=sequences(tl,P)
  
  ## calculate stochastic growth rate
  cat('loglams') #### ***********
  
  n=matrix(1,s,1)/s #create a one-column matrix with equal frequencies for each life stage
  r<-array(0,c(tl))
  
  for (it in 1:tl){
    n=B[,,env[it]]%*%n #B is actually a set of mean matrices
    r[it]=sum(n) #store the population size
    n=n/sum(n) #create a proportion
  }
  loglams=mean(log(r)) #save the stochastic growth rate for each fire scenario
  
  ## calculate sensitivity, using forward and backward iteration
  ## see MPM, Sectioni 14.4.1.1
  cat('forward iteration') #### 
  
  ## forward iteration for stage distribution
  wout<-array(0,c(s,tl+1))
  w=array(1,c(s,1))/s #create an array of vectors of equal frequencies
  wout[,1]=w #fill w into wout in the first column position
  
  for (it in 1:tl){ #iterate forward
    w=B[,,env[it]]%*%w #right side matrix multiplication of w and the appropriate matrix from the set of mean matrices (B)
    r[it]=sum(w) #save the population size to r
    w=w/sum(w) #make w proportional
    wout[,it+1]=w #stor w in the wout object
  } 
  
  cat('backward iteration') #### 
  
  ## backward iteration for repro value
  vout<-array(0,c(s,tl+1))  
  v=array(1,c(1,s))/s #create a matrix of equal initial frequencies
  vout[,tl+1]=t(v) #transpose the vector and place it as the final column in vout
  for (it in tl:1){ #for the backwards iteration
    v=v%*%B[,,env[it]] #multiple v by the appropriate matrix in the sequence (env), this is left multiplication
    v=v/sum(v); #redefine v from above as a proportion
    vout[,it]=t(v) #fill the transpose of this vector into vout
  } 
  
  ## now calculate the environment-specific sensitivities; eqn 16
  cat('env-specific sensitivities') #### 
  
  for (ienv in 1:ne){ #for each environmental state (B1, B2, etc)
    J=as.numeric(env%in%ienv) ## identifies places in sequence with given env value
    sens=J[1]*kronecker(t(wout[,1]),t(vout[,1+1]))%*%dvecBdt[,,env[1]]/(as.numeric(r[1]*t(vout[,1+1])%*%wout[,1+1]))
    
    for (t in 2:tl){
      sens= sens + J[t]*kronecker(t(wout[,t]),t(vout[,t+1]))%*%dvecBdt[,,env[t]]/(as.numeric(r[t]*t(vout[,t+1])%*%wout[,t+1])) 
    } ## for t
    
    dlogldt[,,ienv]=sens/tl #here I think sens is defined as the sum of all sensitivities for each year across the given environment (ne)
  }
  return(dlogldt)
} ## for ienv


#################################################################################
## PART VI: Creating function "beargrass_1" to do SLTRE
#################################################################################

beargrass_1<-function(B_ls, B_hs){
  ## environment-specific matrices, averaged between treatments
  B_mean=(B_ls + B_hs)/2  ## this is an n by n by l array that create average matrices between B_fb and B_rp
  
  ## Time limit for iterations
  tl=100#10000
  
  ## number of stages (s), environments (ne), and treatments (nt)
  s=dim(B)[1] #side length of the input  matrix
  ne=dim(B)[3] #number of environmental stages (B1, UB3, etc.)
  nt=dim(B)[4]#number of treatments, here high and low severity fire 
  
  ## in all calculations
  ## envt 1 = high severity fire (G)
  ## envt 2 = low severity fire (H)
  ## (alphabetical order; it's easier to remember)
  
  ## derivatives of the matrices B wrt the parameters in the paper, the parameters are the elements of B
  dvecBdt=array(0, c(s^2, s^2, ne))
  for(i in 1:ne){
    dvecBdt[,,i]=diag(s^2);
  } # for i
  
  ## specify frequencies to examine, and autocorrelations
  freqarray=seq(0.01, 0.99, length.out=12)
  auto=matrix(0,1,length(freqarray))
  
  #arrays to hold outputs
  loglams_ls<-array(0,c(1,length(freqarray))) 
  loglams_hs<-array(0,c(1,length(freqarray)))
  growth_hs<-array(0,c(length(freqarray),tl))
  growth_ls<-array(0,c(length(freqarray),tl))
  contmats<-array(0,c(bms^2,ne,length(freqarray)))
  contmats_env<-array(0,c(ne,length(freqarray)))
  contmats_aij<-array(0,c(bms,bms,length(freqarray)))
  dlogldt<-array(0,dim=c(1,bms^2,ne))
  #dlogldt_save<-array(0,dim=c(length(freqarray),bms^2,ne))
  
  ## decomposition for each frequency
  P=array(0, c(6,6,length(freqarray)))
  
  for(i in 1:length(freqarray)){
    
    # create environmental transition matrices
    f=freqarray[i]
    rho=auto[i]
    q=f*(1-rho) #q and p are the probability of burning with no autocorrelation, 1-p is the probability of not burning
    p=q+rho 
    
    P[,,i]=matrix(c(p, q, q, q, q, q,
                    1-p, 0, 0, 0, 0, 0,
                    0, 1-q, 0, 0, 0, 0,
                    0, 0, (1-q)/3, (1-q)/3, (1-q)/3, (1-q)/3,
                    0, 0, (1-q)/3, (1-q)/3, (1-q)/3, (1-q)/3,
                    0, 0, (1-q)/3, (1-q)/3, (1-q)/3, (1-q)/3
    ),6,6, byrow=T) #when a fire occurs, the sequence cycles through 1, 2 and 3 years post-fire matrices (unless another fire occurs) and then to the unburned state. All unburned matrices have equal probability of occuring. q is the probability of fire if it did not burn in the previous year, 1-q is the probability of not burning. 
    # evaluate stochastic growth rates
    
    # generate environmental state sequences
    env=sequences(tl, P[,,i]) ## new function
    
    # growth rate for HS at each fire freq
    loglams_hs[i]<-stochgrow(env,B_hs)[1]
    
    # growth rate for LS at each fire freq
    loglams_ls[i]<-stochgrow(env,B_ls)[1]
    
    growth_hs[i]<-stochgrow(env,B_hs)[2]
    growth_ls[i]<-stochgrow(env,B_ls)[2]
    
    # environment-specific sensitivities to the vital rates
    dlogldt=stocsens_e(P[,,i],B_mean,dvecBdt,tl)
    # dlogldt dimensions 1 x number of parametrs x number of envts
    
    # save the sensitivities for all fire freqs
    #dlogldt_save[,,i]<-drop(dlogldt) #this doesn't work/not sure I need it
    
    # contributions from vital rate differences
    for (ienv in 1:ne){ # contmats dim 1 x no parameters x no envts x no. of freq
      contmats[,ienv,i]<-t(vec(B_ls[,,ienv]-B_hs[,,ienv]))*dlogldt[,,ienv]              #the vec function vectorizes the matrix, column by column. After being transposed, its dimentions are 1 X bms^2. I subtracted B_hs from B_ls because B_ls is larger.           
    } # for ienv
    
    # sum contributions over matrix elements - for each i, which is fire frequency
    contmats_env[,i]<-apply(contmats[,,i],2,sum) 
    
    #sum contributions over envts
    contmats_aij[,,i]<-matrix(apply(contmats[,,i],1,sum),bms,bms) ## fix reshape ***
  }
  return(list(loglams_hs, loglams_ls,  contmats, contmats_env, contmats_aij))
}    
## This function returns:
## loglams_hs (growth rate for Hs at each fire freq)
## loglams_ls (growth rate for Ls at each fire freq)
## contmats (contributions from vital rate differences)
## contmats_env (sum contributions over matrix elements) 
## contmats_aij (sum contributions over envts)

#comparing high and low severity fire
## with and without harvest
SLTRE_list_fire<-beargrass_1(B_ls, B_hs) 

names(SLTRE_list_fire)<-c("loglams_hs",
                          "loglams_ls",
                          "contmats",
                          "contmats_env",
                          "contmats_aij")

#save results
saveRDS(object=SLTRE_list_fire, file="SLTRE_fire_output.rds")



##Author: Holly K. Kindsvater, 2016

    
#Model parameters:
Npop=100  #initial pop size (not important)
Tmax=500 #max years
Tfishing=200 #start fishing
Trecovery=400 #stop fishing 
  
    #Life history parameters      	
   recruitsize=30
    Amax=18 #max age
 	amat=5 #age at maturation - this varies depending on the matration function chosen
	Linf=299.5 #max size
	 fishingsize = .75*Linf #size-selective fishery threshold size 
        
       k=.12   # Growth function param	 	 
	   v = 50 #mass-at-size coef
	    wexp = 3  #mass-at-size exp
	    h=1 #sperm fertilization effectiveness
 	    N0 = 20000       
  	    alpha = .9  	
  	    q=.8 #maturation ogive steepness (1 = very steep)
  
   mu_f=0.109 #fishing mortality
     natmort=0.087
 
      
    
    #set up vectors and matrices
pmat=rep(0, Amax) 
select=rep(0, Amax)
pups=rep(0, Amax)
mu=rep(0, Amax)
N=matrix(nrow=Amax, ncol=Tmax, data=0)
 
W=matrix(nrow=Amax, ncol=1, data=1)
L=matrix(nrow=Amax, ncol=1, data=1)
mu=matrix(nrow=Amax, ncol=1, data=0)
E=matrix(nrow=Amax, ncol=Tmax, data=0)
S=matrix(nrow=1, ncol=Tmax, data=0)
P=matrix(nrow=1, ncol=Tmax)
B=matrix(nrow=1, ncol=Tmax)
 
 
  
 N[, 1]=Npop #initial population size
L[1]=recruitsize #larval size at recruitment (birth)
   
#First: Define the age-specific growth, mortality, and maturation function
  
  	for (a in 1:(Amax-1)) {
  	
  	
	#GROWTH 		 
	L[a+1]= Linf*(1-exp(-k)) + (L[a])*exp(-k) #length at age
     W[a]= v*L[a]^wexp #weight at age
     
  
	# AGE-DEP MATURATION:
       pmat[a]= 1/(1+exp(-q*(a-amat))) 
           
     
 	#FISHING SELECTIVITY # knife edge after age 2
      select[a] = ifelse(a < 2, 0, 1)     
           
 	#FECUNDITY
       pups[a] = 0.5  #mobulid fecundity is on average 1 female per year
       
    #MORTALITY
       mu[a] = natmort #could vary with age but assumed constant here
 
   } #end 1st age loop
 
 #With these age specific functions, simulate pop dynamics through time       
 for(t in 1:(Tmax-1)) {
	 
   
	E[,t]=N[, t]*pmat*pups  #assuming reproduction occurs between 1 t and the next
 
	P[t]= sum(E[,t])  #assuming fertilization is 100%, this is the year class of progeny produced in year t
  
  	beta= (alpha*P[t]-1)/(N0*P[t])  #derive the recruitment function from the slope near the origin
    
 	N[1,t+1]= alpha*P[t]/(1+beta*P[t]) #this is the recruit class that enters in the next time step.. 
 
     for (age in 1:(Amax-1)) { #now calculate age-specific probability of fishing mortality
      
   Fishing= if (Tfishing < t & Trecovery > t) {
   				Fishing = select[age]*mu_f } else {
   				Fishing = 0
                     }
    
 		N[age+1,t+1] = N[age, t]*exp(-mu[age]-Fishing) #numbers in each age class in the next time step 
 		#assuming mortality happens after spawning

 		B[t]=sum(N[,t]*W) #biomass
  			
   		 } #end second age loop
    
   
       } #end t loop
 
 
 
 steadystate=N[,Tmax]
 fishedsteady=N[ ,Trecovery -2]
   
quartz()
par(mfrow=c(2, 2),  las=1)
 
barplot(steadystate, names.arg=c(1:Amax ), main="Unfished Age Structure", ylab="Abundance" )
barplot(fishedsteady, names.arg=c(1:Amax), col=1, main="Constant selectivity for A > 2") 

barplot(log(steadystate), names.arg=c(1:Amax), ylab="Ln(Abundance)", xlab="Age")
barplot(log(fishedsteady), names.arg=c(1:Amax), col=1, xlab="Age") 

######Plot life history functions and population dynamics 
 quartz()
 par(mfrow=c(2,2))
   plot(1:Amax, L, type="l", main="von Bert growth", ylab="Disc Width (cm)", xlab="Age", ylim=c(0, 300)) #growth curve
    legend("bottomright", legend=c( paste("Linf = ", Linf), paste("von Bert k = ", k)))

  
plot(1:(Amax-1), pmat[-Amax], type="l", main="Maturation ogive", ylab = "Probability Mature", xlab="Age")
 legend("bottomright", legend=c(paste("50% mature at age", amat)))
 plot(1:(Amax-1), pups[-Amax]*pmat[-Amax], type="l", main="Age-specific Fecundity", ylab="Number of female progeny", xlab= "Age")
    legend("bottomright", legend=c(paste("50% mature at age", amat)))

   
  P_= 1:N0
  R=alpha*P_/(1+beta*P_)
   plot(P_, R, type="l", main="Beverton Holt Recruitment", xlab="Female progeny[t]", ylab="Recruits[t+1]")
    legend("bottomright", legend=c(paste("alpha = ", round(alpha,6)), paste("beta = ", round(beta, 6))))


     # plot(1:Tmax, colSums(N), type="l", lwd=1.5, main="Population Dynamics", ylab="Population Size", xlab="Time")
   # # plot(sums[1:(Tmax-1)], N[amat,2:Tmax], main="Stock Recruitment Curve", xlab="Spawners", ylab="1st time spawners", type="l", lwd=1.5)
   
   
   
  
       
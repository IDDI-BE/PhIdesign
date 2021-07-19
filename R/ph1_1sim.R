
#' @title Simulation of one Phase I study
#' @description Simulation of one Phase I study with classical 3+3 design, 
#'     BOIN design or Keyboard design
#' @param phi target toxicity rate (needed for all designs)
#' @param phi1 for BOIN design: the highest DLT rate that is deemed subtherapeutic (i.e., underdosing), 
#'     such that dose escalation should be made, see \code{\link{get_BOIN_rules}}
#' @param phi2 for BOIN design: the lowest DLT rate that is deemed overly toxic (i.e., overdosing), 
#'     such that dose de-escalation is required, see \code{\link{get_BOIN_rules}}
#' @param halfkey for Keyboard design: half length of key in keyboard design, see \code{\link{get_Keyboard_rules}}
#' @param maxtox P(DLT rate>phi|c,n)<=maxtox as value for maxprob in \code{\link{get_Elim_rules}}
#' @param N number of patients (not needed for design="3+3")
#' @param truerate scenario of true DLT rates for each dose level. This also defines the number of dose levels
#' @param cohortsize cohort size, default is 3
#' @param maxNretain If N(patients) at next dose level >=maxNretain, then stop algorithm, default is 9
#' @param acc_tit logical indicator if algorithm should start with accelerated titration=algorithm starts at dose level with first DLT
#' @param design "BOIN" or "Keyboard"
#' @param MTD_safer logical indicator for ...
#' @param seed define seed
#' @param sim if "NO", calculate decision rules, define "YES" if done within simulation program for multiple studies
#' @param env parent environment to pass objects from \code{\link{ph1_sim_OC}}
#' @return a data.frame with elements
#' \itemize{
#' \item dose_level: dose level
#' \item dose: actual dose to be used
#' \item incr_ratio: The ratio of successive doses
#' }
#' @references Liu S, Yuan Y. Bayesian Optimal Interval Designs for Phase I Clinical Trials. 
#'     J R Stat Soc Ser C Appl Stat. 2015;64:507–23
#'     Yan F, Mandrekar SJ, Yuan Y. Keyboard: A Novel Bayesian Toxicity Probability Interval Design for Phase I Clinical Trials. Clin Cancer Res. 2017; 23: 3994–4003
#' @export
#' @examples
#' ph1_1sim(phi=0.3, phi1=0.6*0.3, phi2=1.4*0.3, maxtox=0.95, N=15, 
#' truerate=c(0.20 ,0.25 ,0.30 ,0.40 ,0.60 ,0.75), cohortsize=3, maxNretain=9, acc_tit=0, 
#' design="BOIN", MTD_safer=TRUE)
#' ph1_1sim(phi=0.3, phi1=0.6*0.3, phi2=1.4*0.3, maxtox=0.95, N=30, 
#' truerate=c(0.20 ,0.25 ,0.30 ,0.40 ,0.60 ,0.75), cohortsize=3, maxNretain=9, acc_tit=1, 
#' design="Keyboard", MTD_safer=TRUE, halfkey=0.05)
#' ph1_1sim(truerate=c(0.20 ,0.25 ,0.30 ,0.40 ,0.60 ,0.75),acc_tit=1,design="3+3")

# phi=args$phi;phi1=args$phi_1;phi2=args$phi2;maxtox=args$maxtox; N=args$N;truerate=args$truerate;cohortsize=args$cohortsize;maxNretain=args$maxNretain;acc_tit=args$acc_tit; design=args$design; MTD_safer=args$MTD_safer; halfkey=args$halfkey; sim=args$sim

ph1_1sim <- function (phi, phi1=NULL, phi2=NULL, maxtox=NULL, N=NULL,truerate=NULL,cohortsize=NULL,maxNretain=NULL, acc_tit, 
                      design, MTD_safer=NULL, seed = NULL, halfkey=NULL, sim="NO", env=parent.frame()){ # MTD_safer: MTD should be <lambda_d
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  #-----------------------------------------------------------#
  # First dose: Create empty vectors and initialize variables #
  #-----------------------------------------------------------#

  ndose<-length(truerate)
  
  dose_in_the_running <- dose <- 1:ndose
  npt  <- ndlt <- MTD  <- rep(0,ndose)
  p_iso<-rep(NA,ndose)
  
  dose_inv <- 1 # initialize "dose under investigation"
  n        <- 0 # number of patients in study, for first dose
  mtd      <- NA # initializing scalar for dose level equal to MTD
  
  #--------------------------------------------------------------------------------------------------#
  # if sim=="NO" (standard option, if code not run within simulation function), determine boundaries #
  #--------------------------------------------------------------------------------------------------"
  
  if (sim=="NO"){
    
    if (! design=="3+3"){
      DLT_STOP<-get_Elim_rules(Nmax=N,phi=phi,maxprob=maxtox)$c_STOP
    }
    
    if (design=="Keyboard"){
      key_dec <-get_Keyboard_rules(N=N,phi=phi,halfkey=halfkey)
    }
    
    if (design=="BOIN"){
      
      BOIN_thres <- get_BOIN_rules(phi=phi,phi1=phi1,phi2=phi2)
      lambda_e <- BOIN_thres$lambda_e
      lambda_d <- BOIN_thres$lambda_d
    }
  }
  
  if (sim=="YES"){
    
    if (! design=="3+3"){
      DLT_STOP<-env$DLT_STOP
    }
    
    if (design=="Keyboard"){
      key_dec <-env$key_dec
    }
    
    if (design=="BOIN"){
      lambda_e<-env$lambda_e
      lambda_d<-env$lambda_d
    }

  }
  
  #-----------------------------------------#
  # First dose: Accelerated Titration Phase #
  #-----------------------------------------#
  
  if (acc_tit==1){
    repeat{
      
      n              <- n+1
      dlt            <- rbinom(n=1,size=1,truerate[dose_inv])
      npt[dose_inv]  <- 1
      ndlt[dose_inv] <- dlt
      
      if (dlt==0){
        if (dose_inv!= ndose){dose_inv <- dose_inv + 1} else # Stay in accelerated titration phase
          if (dose_inv== ndose){break}                         # Go to cohort-wise dose-escalation algorithm
      }
      if (dlt==1){break}  # Go to cohort-wise dose-escalation algorithm
    }
  }
  
  #print(cbind(dose,npt,ndlt,MTD)) # for debugging
  
  #------------#
  # 3+3 design #
  #------------#
  
  if (design=="3+3"){
    
    repeat{
      
      #--------------------#
      # Get data by cohort #
      #--------------------#
      
      if (npt[dose_inv]==1){  # scenario of accelerated titration with DLT or no DLT at last dose
        npt[dose_inv]  <- npt[dose_inv]  + 2
        ndlt[dose_inv] <- ndlt[dose_inv] + rbinom(n=1,size=2,truerate[dose_inv])
      } else {
        npt[dose_inv]  <- npt[dose_inv]  + 3
        ndlt[dose_inv] <- ndlt[dose_inv] + rbinom(n=1,size=3,truerate[dose_inv])
      }
      
      p_hat<- ndlt/npt
      
      #---------------------------------------------------#
      # Get correct decision for dose under investigation #
      #---------------------------------------------------#
      
      phat<-p_hat[dose_inv]
      
      # phat<1/3
      #---------
      
      if (phat< 1/3 &   (dose_inv+1) %in% dose_in_the_running) {dose_inv <- dose_inv+1} else
        if (phat< 1/3 & ! (dose_inv+1) %in% dose_in_the_running) {
          
          if (npt[dose_inv]<=3) {dose_inv <- dose_inv} else
            if (npt[dose_inv]==6) {mtd<- dose_inv
            MTD[dose_inv]<-1
            break}
        } else
          
          # phat>1/3
          #---------
      
      if (phat> 1/3 & (dose_inv-1) %in% dose_in_the_running){
        
        if (npt[dose_inv-1]<=3) {dose_in_the_running <- dose_in_the_running[dose_in_the_running < dose_inv]
        dose_inv <- dose_inv-1} else
          if (npt[dose_inv-1]==6) {mtd<- dose_inv-1
          MTD[dose_inv-1]<-1
          break}
        
      } else
        
        if (phat> 1/3 & ! (dose_inv-1) %in% dose_in_the_running) {break} else
          
          # phat=1/3
          #---------
      
      if (phat==1/3 & npt[dose_inv]==3)  {dose_inv <- dose_inv} else
        if (phat==1/3 & npt[dose_inv]==6 & (dose_inv-1) %in% dose_in_the_running){
          
          if (npt[dose_inv-1]<=3){dose_in_the_running <- dose_in_the_running[dose_in_the_running < dose_inv]
          dose_inv <- dose_inv-1} else
            
            if (npt[dose_inv-1]==6)  {mtd<- dose_inv-1
            MTD[dose_inv-1]<-1
            break}
          
        } else
          
          if (phat==1/3 & npt[dose_inv]==6 & ! (dose_inv-1) %in% dose_in_the_running) {break}
      
      #print(cbind(dose,npt,ndlt,MTD,p_hat)) # for debugging
    }
    
    #cbind(dose,npt,ndlt,MTD,p_hat)
    
  }
  
  #--------------#
  # BOIN design  #
  #--------------#
  
  if (design=="BOIN"){
    
    while ((N-n)>= cohortsize) {
      
      #--------------------#
      # Get data by cohort #
      #--------------------#
      
      if (npt[dose_inv]==1){  # scenario of accelerated titration with DLT or no DLT at last dose
        npt[dose_inv]  <- npt[dose_inv]  + cohortsize-1
        ndlt[dose_inv] <- ndlt[dose_inv] + rbinom(n=1,size=(cohortsize-1),truerate[dose_inv])
        n              <- n+(cohortsize-1)
      } else {
        npt[dose_inv]  <- npt[dose_inv]  + cohortsize
        ndlt[dose_inv] <- ndlt[dose_inv] + rbinom(n=1,size=(cohortsize  ),truerate[dose_inv])
        n              <- n+(cohortsize)
      }
      
      p_hat<- ndlt/npt
      
      #----------------------------------------------------------------------------#
      # Check if current dose and doses above should be eliminated due to toxicity #
      #----------------------------------------------------------------------------#
      
      if(ndlt[dose_inv]>=DLT_STOP[npt[dose_inv]]){dose_in_the_running <- dose_in_the_running[dose_in_the_running < dose_inv] }
      
      #---------------------------------------------------#
      # Get correct decision for dose under investigation #
      #---------------------------------------------------#
      
      phat <-p_hat[dose_inv]
      
      if (phat<=lambda_e)               {decision<-"E"}  # Escalate
      if (phat>=lambda_d)               {decision<-"D"}  # De-escalate
      if (lambda_e<phat & phat<lambda_d){decision<-"R"}  # Retain
      
      
      if (decision=="E") {
        if ( (dose_inv+1) %in% dose_in_the_running){dose_inv <- dose_inv+1} else
          if (!(dose_inv+1) %in% dose_in_the_running){decision<-"R" } # If no higher dose available: RETAIN DOSE
      }
      
      if (decision=="D") {
        if ( (dose_inv-1) %in% dose_in_the_running){dose_inv <- dose_inv-1} else
          if (!(dose_inv-1) %in% dose_in_the_running){break  } # If no lower dose available: STOP
      }
      
      if (decision=="R") {
        if ( dose_inv %in% dose_in_the_running){dose_inv <- dose_inv} else {break}
      }
      
      if (npt[dose_inv]>=maxNretain){break} # if number of patients on "dose under investigation">= max specified (regardless decision E/D/R): STOP
      
    }
    
    #----------------------------------------------------------------#
    # Calculate isotonic estimate for p(DLT) (function BOIN package) #
    #----------------------------------------------------------------#
    
    if (MTD_safer==T){
      
      p_iso[npt>0]  <- as.numeric(as.character(BOIN::select.mtd(target=phi, npts=npt[npt>0], ntox=ndlt[npt>0])$p_est$phat)) # Get isotonic estimates
      # which lower than lambda_d
      safe_choice <- which(p_iso<lambda_d)
      # of those, which closest to phi
      if (length(safe_choice)!=0){
        distance<- abs(p_iso[safe_choice]-phi)
        mtd<-max(which(distance==min(distance))) # If more than two (under and above MTD), take highest.
      }
      
      if (length(mtd)!=0 & ! is.na(mtd)){
        MTD[mtd]<-1 # MTD: isotonic estimate < lambda_e
      }
      
    } else {
      mtd<-BOIN::select.mtd(target=phi, npts=npt[npt>0], ntox=ndlt[npt>0], cutoff.eli=maxtox)$MTD  # Note if no MTD, function puts it at '99'
      if (!(mtd==99)){MTD[mtd]<-1} else
        if ( (mtd==99)){mtd=NA}
    }
    
    #cbind(dose,npt,ndlt,p_hat,p_iso,MTD)
  }
  
  #------------------#
  # Keyboard design  #
  #------------------#
  
  if (design=="Keyboard"){
    
    while ((N-n)>= cohortsize) {
      
      #--------------------#
      # Get data by cohort #
      #--------------------#
      
      if (npt[dose_inv]==1){  # scenario of accelerated titration with DLT or no DLT at last dose
        npt[dose_inv]  <- npt[dose_inv]  + cohortsize-1
        ndlt[dose_inv] <- ndlt[dose_inv] + rbinom(n=1,size=(cohortsize-1),truerate[dose_inv])
        n              <- n+(cohortsize-1)
      } else {
        npt[dose_inv]  <- npt[dose_inv]  + cohortsize
        ndlt[dose_inv] <- ndlt[dose_inv] + rbinom(n=1,size=(cohortsize  ),truerate[dose_inv])
        n              <- n+(cohortsize)
      }
      
      p_hat<- ndlt/npt # not needed for algorithm
      
      #----------------------------------------------------------------------------#
      # Check if current dose and doses above should be eliminated due to toxicity #
      #----------------------------------------------------------------------------#
      
      if(ndlt[dose_inv]>=DLT_STOP[dose_inv]){dose_in_the_running <- dose_in_the_running[dose_in_the_running < dose_inv] }
      
      #---------------------------------------------------#
      # Get correct decision for dose under investigation #
      #---------------------------------------------------#
      
      decision<-key_dec[key_dec$n==npt[dose_inv] & key_dec$x==ndlt[dose_inv],"decision"]
      
      if (decision=="E") {
        if ( (dose_inv+1) %in% dose_in_the_running){dose_inv <- dose_inv+1} else
          if (!(dose_inv+1) %in% dose_in_the_running){decision<-"R" } # If no higher dose available: RETAIN DOSE
      }
      
      if (decision=="D") {
        if ( (dose_inv-1) %in% dose_in_the_running){dose_inv <- dose_inv-1} else
          if (!(dose_inv-1) %in% dose_in_the_running){break  } # If no lower dose available: STOP
      }
      
      if (decision=="R") {
        if ( dose_inv %in% dose_in_the_running){dose_inv <- dose_inv} else {break}
      }
      
      if (npt[dose_inv]>=maxNretain){break} # if number of patients on "dose under investigation">= max specified (regardless decision E/D/R): STOP
      
    }
    
    #----------------------------------------------------------------#
    # Calculate isotonic estimate for p(DLT) (function BOIN package) #
    #----------------------------------------------------------------#
    
    if (MTD_safer==T){
      
      p_iso[npt>0]  <- as.numeric(as.character(BOIN::select.mtd(target=phi, npts=npt[npt>0], ntox=ndlt[npt>0])$p_est$phat)) # Get isotonic estimates
      # which lower than lambda_d
      safe_choice <- which(p_iso<phi+halfkey)
      # of those, which closest to phi
      if (length(safe_choice)!=0){
        distance<- abs(p_iso[safe_choice]-phi)
        mtd<-max(which(distance==min(distance))) # If more than two (under and above MTD), take highest.
      }
      
      if (length(mtd)!=0 & ! is.na(mtd)){
        MTD[mtd]<-1 # MTD: isotonic estimate < lambda_e
      }
      
    } else {
      mtd<-BOIN::select.mtd(target=phi, npts=npt[npt>0], ntox=ndlt[npt>0], cutoff.eli=maxtox)$MTD  # Note if no MTD, function puts it at '99'
      if (!(mtd==99)){MTD[mtd]<-1} else
        if ( (mtd==99)){mtd=NA}
    }
    
    #cbind(dose,npt,ndlt,p_hat,p_iso,MTD)
  }
  
  #----------------#
  # Return results #
  #----------------#
  
  return(as.data.frame(cbind(dose,npt,ndlt,p_hat,p_iso,MTD)))
  
}



#' @title Simulation of Phase I study operating characteristics
#' @description Simulation of Phase I study operating characteristics 
#' @param nsim number of simulations
#' @param scenarios list of scenario's (each vector: set of true DLT rate for each dose level). 
#' Examples of scenarios can be found in a vignette accompanying this package. 
#' For more details see the help vignette: \code{vignette("scenarios", package = "PhIdesign")}
#' @param env parent environment (global) to pass scenario counter to global environment to print progress
#' @param OD_cut % of patients > OD_cut with overdose for overdose_1 and overdose_2 outcomes; should be 
#' vector of 2. Default value is c(0.6,0.8)
#' @param OD_delta dose level with toxicity >phi+OD_delta is considered as overdose. Default value is zero. 
#' @param ...  see \code{\link{ph1_1sim}}
#
#' @return a list containing of a matrix (results by dose level) and a vector (overall results)
#' \itemize{
#' \item [[1]]$dose_level: dose level
#' \item [[1]]$truerate: true DLT rate
#' \item [[1]]$sel_pct_dose: selection proportion at each dose level
#' \item [[1]]$avg_SS_dose: average number of patients treated by dose level
#' \item [[1]]$nDLT_dose: average number of DLT's observed at dose level
#' \item [[2]]$nsim: number of simulations
#' \item [[2]]$design: design used
#' \item [[2]]$sel_at_phi_pct: percent simulated trials in which the correct dose is selected, i.e., the dose with DLT rate equal to phi
#' \item [[2]]$sel_at_phi_rng_pct: percent simulated trials in which the selected dose has a DLT rate that lies in the interval [phi-0.05, phi+0.05]
#' \item [[2]]$sel_below_phi_pct: percent simulated trials in which the selected dose has a DLT rate < phi-0.05
#' \item [[2]]$sel_above_phi_pct: percent simulated trials in which the selected dose has a DLT rate > phi+0.05
#' \item [[2]]$avg_npt_MTD_sel: average number of patients treated at the actual selected MTD for each simulated study (studies without MTD selection not taken into account)
#' \item [[2]]$avg_pct_pts_at_phi: average percent patients in each simulated trial treated with dose with DLT rate =phi
#' \item [[2]]$avg_pct_pts_at_phi_rng: average percent patients in each simulated trial treated with dose with DLT rate that lies in the interval [phi-0.05, phi+0.05]
#' \item [[2]]$avg_pct_pts_below_phi: average percent patients in each simulated trial treated with dose with DLT rate < phi-0.05
#' \item [[2]]$avg_pct_pts_above_phi: average percent patients in the simulated trials treated with dose with DLT rate > phi+0.05
#' \item [[2]]$noMTD: percent simulated trials without MTD selection
#' \item [[2]]$avg_nDLT: average number of observed DLTs
#' \item [[2]]$avg_SS: average sample size of all simulated trials
#' \item [[2]]$max_SS: maximum sample size of all simulated trials
#' \item [[2]]$overdose_1: percent simulated trials in which >=OD_cut[1] of patients is assigned to a dose with DLT rate >phi
#' \item [[2]]$overdose_2: percent simulated trials in which >=OD_cut[2] of patients is assigned to a dose with DLT rate >phi
#' \item [[3]] labels for item[[2]]
#' }
#' @references Liu S, Yuan Y. Bayesian Optimal Interval Designs for Phase I Clinical Trials. 
#'     J R Stat Soc Ser C Appl Stat. 2015;64:507–23
#'     Yan F, Mandrekar SJ, Yuan Y. Keyboard: A Novel Bayesian Toxicity Probability Interval Design for Phase I Clinical Trials. Clin Cancer Res. 2017; 23: 3994–4003
#' @export
#' @examples
#' result<-ph1_sim_OC(nsim=5,  scenarios=list(c(0.40 ,0.50 ,0.60 ,0.70 ,0.80 ,0.90),
#'                                            c(0.30 ,0.40 ,0.50 ,0.60 ,0.70 ,0.80),
#'                                            c(0.30 ,0.35 ,0.40 ,0.60 ,0.75 ,0.90),
#'                                            c(0.20 ,0.30 ,0.40 ,0.50 ,0.60 ,0.70),
#'                                            c(0.20 ,0.30 ,0.50 ,0.65 ,0.80 ,0.90),
#'                                            c(0.20 ,0.25 ,0.30 ,0.40 ,0.60 ,0.75),
#'                                            c(0.10 ,0.20 ,0.30 ,0.40 ,0.50 ,0.60),
#'                                            c(0.10 ,0.15 ,0.25 ,0.30 ,0.50 ,0.75),
#'                                            c(0.05 ,0.10 ,0.15 ,0.20 ,0.25 ,0.30),
#'                                            c(0.05 ,0.10 ,0.20 ,0.30 ,0.50 ,0.65),
#'                                            c(0.025,0.05 ,0.10 ,0.20 ,0.30 ,0.40)), 
#'    phi=0.3, phi1=0.6*0.3, phi2=1.4*0.3, maxtox=0.95, N=18, 
#'    cohortsize=3, maxN=9, maxNdec="STAY",acc_tit=0, design="BOIN", MTD_safer=TRUE)

# nsim=100;scenarios=list(c(0, 0.05 ,0.15));env=parent.frame() # for debug
ph1_sim_OC<-function(nsim,scenarios,env=parent.frame(),OD_cut=c(0.6,0.8),OD_delta=0,...){ 
  
  res<-list()
  
  for (s in 1:length(scenarios)){

    DLT_STOP<-key_dec<-lambda_e<-lambda_d<-NULL
    
    # Create scenario counter
    #------------------------
    
    if (exists("scenario_sim")==F){
      env$scenario_sim <-0
    } # for printing progress in case multiple scenarios are tested (needs to be in global environment)
  
    env$scenario_sim <- env$scenario_sim+1
    
    # Get ... arguments to be used by 'ph1_1sim' funtion
    #-------------------------------------------------
    # args<-list(phi=0.1,phi1=0.2*0.1,phi2=2.2*0.1, maxtox=0.95, N=15, cohortsize=3, maxN=9, acc_tit=0, design="3+3", MTD_safer=TRUE);s=1 # for debug
    args<-list(...) # these are the arguments for the 'ph1_1sim' function (see function arguments above, and use with "args" here below)
    args$sim<-"YES" # Don't calculate thresholds in ph1_1sim function, rather only once in sim function
    
    # Get one set of true DLT rates
    #------------------------------
    args$truerate<-scenarios[[s]] 
    
    # Set empty matrix and vector to fill with simulation results
    #------------------------------------------------------------
    
    npt<-ndlt<-MTD<-matrix(nrow=nsim,ncol=length(args$truerate))
    npt_MTDselect<-overdose_pct<-n<-rep(NA,nsim)
    
    # Get decision thresholds
    #------------------------
    
    if (! args$design=="3+3"){
      DLT_STOP <-get_Elim_rules(Nmax=args$N,phi=args$phi,maxprob=args$maxtox)$c_STOP
    }
    
    if (args$design=="Keyboard"){
      key_dec <-get_Keyboard_rules(N=args$N,phi=args$phi,halfkey=args$halfkey)
    }
    
    if (args$design=="i3+3"){
      i33_dec <-get_i33_rules(N=args$N,phi=args$phi,phi1=args$phi1,phi2=args$phi2)
    }
    
    if (args$design=="BOIN"){
      BOIN_thres <- get_BOIN_rules(phi=args$phi,phi1=args$phi1,phi2=args$phi2)
      lambda_e <- BOIN_thres$lambda_e
      lambda_d <- BOIN_thres$lambda_d
    }
    
    # Run nsim simulations
    #---------------------
    
    for (i in 1:nsim){
      print(paste0("scenario ",env$scenario_sim," sim:",i,"/",nsim))
      
      # run one simulation

      result <- ph1_1sim(phi=args$phi,phi1=args$phi1,phi2=args$phi2,maxtox=args$maxtox, N=args$N,truerate=args$truerate,cohortsize=args$cohortsize,maxN=args$maxN,maxNdec=args$maxNdec,
                         acc_tit=args$acc_tit, dose_no_titr=args$dose_no_titr, BOIN_add33_rule=args$BOIN_add33_rule,design=args$design, MTD_safer=args$MTD_safer, halfkey=args$halfkey, sim=args$sim)
      result
      npt [i,]         <- result[,"npt"]
      ndlt[i,]         <- result[,"ndlt"]
      MTD [i,]         <- result[,"MTD"]
      overdose_pct [i] <- sum(result[args$truerate>args$phi+OD_delta,"npt"])/sum(result[,"npt"]) # %(patients) in one trial which is overdosed (i.e. treated with dose level with toxicity rate above target)
      
      if ((sum(result[,"MTD"])!=0)) {npt_MTDselect[i] <- result[result[,"MTD"]==1,"npt"]} else {npt_MTDselect[i]<-0} # Number of patients treated at selected MTD (not target rate)
    }
    
    # Summarize simulation results
    #-----------------------------
    
    # Summary by dose level
    sel_pct_dose <- apply(MTD ,2,mean)
    avg_SS_dose  <- apply(npt ,2,mean)
    nDLT_dose    <- apply(ndlt,2,mean)
    
    # Summary overall
    if (args$phi %in% args$truerate){
      sel_at_phi_pct     <- sel_pct_dose[args$truerate==args$phi]
      avg_pct_pts_at_phi <- mean(npt[,   args$truerate==args$phi]/apply(npt,1,sum))} else {
      sel_at_phi_pct     <- NA  # % Selection at target
      avg_pct_pts_at_phi <- NA} # Average %(patients) at target
    
    if (any((args$phi-0.05) <= args$truerate & (args$phi+0.05)>=args$truerate)){
      sel_at_phi_rng_pct     <- sum(sel_pct_dose[args$truerate>=(args$phi-0.05) & args$truerate<=(args$phi+0.05)])
      avg_pct_pts_at_phi_rng <- mean(apply(npt[, args$truerate>=(args$phi-0.05) & args$truerate<=(args$phi+0.05),drop=FALSE],1,sum)/apply(npt,1,sum))} else {
      sel_at_phi_rng_pct     <- NA  # % Selection at target +/-0.05
      avg_pct_pts_at_phi_rng <- NA} # Average %(patients) at target +/-0.05
    
    if (any(args$truerate<args$phi-0.05)){
      sel_below_phi_pct     <- sum(sel_pct_dose[args$truerate<(args$phi-0.05)])
      avg_pct_pts_below_phi <- mean(apply(npt[, args$truerate<(args$phi-0.05),drop=FALSE],1,sum)/apply(npt,1,sum))} else {
      sel_below_phi_pct     <- NA  # % Selection at < target -0.05
      avg_pct_pts_below_phi <- NA} # Average %(patients) at < target -0.05
    
    if (any(args$truerate>args$phi+0.05)){
      sel_above_phi_pct     <- sum(sel_pct_dose[args$truerate>(args$phi+0.05)])
      avg_pct_pts_above_phi <- mean(apply(npt[, args$truerate>(args$phi+0.05),drop=FALSE],1,sum)/apply(npt,1,sum))} else {
      sel_above_phi_pct     <- NA  # % Selection at > target +0.05
      avg_pct_pts_above_phi <- NA} # Average %(patients) at > target +0.05   
    
    noMTD            <- sum(apply(MTD,1,sum)==0)/nsim # % No MTD selection
    avg_SS           <- mean(apply(npt,1,sum))        # Average sample size
    max_SS           <- max(apply(npt,1,sum))         # Max sample size
    avg_npt_MTD_sel  <- sum(npt_MTDselect[npt_MTDselect>0])/length(npt_MTDselect[npt_MTDselect>0]) # Average N(patients at selected MTD)
    avg_nDLT         <- mean(apply(ndlt,1,sum))       # Average N(patients with DLT)
    overdose_1       <- sum(overdose_pct>=OD_cut[1])/nsim
    overdose_2       <- sum(overdose_pct>=OD_cut[2])/nsim
    
    # Collect simulation results in one list
    #---------------------------------------
    
    res[[s]]<-list(data.frame(dose_level=1:length(args$truerate),truerate=args$truerate,sel_pct_dose=sel_pct_dose,avg_SS_dose=avg_SS_dose,nDLT_dose=nDLT_dose),
                c(design=args$design,nsim=nsim,
                  avg_SS=avg_SS,max_SS=max_SS,
                  avg_npt_MTD_sel=avg_npt_MTD_sel,avg_nDLT=avg_nDLT,
                  sel_at_phi_pct=sel_at_phi_pct,sel_at_phi_rng_pct=sel_at_phi_rng_pct,sel_below_phi_pct=sel_below_phi_pct,sel_above_phi_pct=sel_above_phi_pct,noMTD=noMTD,
                  avg_pct_pts_at_phi=avg_pct_pts_at_phi,avg_pct_pts_at_phi_rng=avg_pct_pts_at_phi_rng,avg_pct_pts_below_phi=avg_pct_pts_below_phi,avg_pct_pts_above_phi=avg_pct_pts_above_phi,
                  overdose_1=overdose_1,overdose_2=overdose_2),
                 c("design","nsim",
                   "Average sample size","Max sample size",
                   "Average N(patients) at selected MTD","Average N(patients) with DLT",
                   "% Selected MTD=at target","% Selection at target +/-0.05","% Selection at < target -0.05","% Selection at > target +0.05","% No MTD selection",
                   "Average %(patients) at target","Average %(patients) at target +/-0.05","Average %(patients) at <target -0.05","Average %(patients) at >target +0.05",
                   paste0("Risk of overdosing >=",OD_cut[1]*100,"% of patients"),paste0("Risk of overdosing >=",OD_cut[2]*100,"% of patients")))
    
  }
  
  # Return results
  #---------------
  
  return(res)
  
}
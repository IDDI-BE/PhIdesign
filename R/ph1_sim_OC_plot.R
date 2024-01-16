
#' @title Plot operating characteristics of phase 1 design
#' @description Plot operating characteristics of phase 1 design. First design is "blue", second one"red",
#'     and third one "darkgreen"
#' @param ph1_sim_OC_output list of output objects of \code{\link{ph1_sim_OC}}
#' @param ylim_N max N for plots with number of patients (typically max sample size of 3+3 design)
#' @param save_as where to save output
#' @param outcomes select outcomes from ph1_sim_OC_output (item 2)
#' @param layout layout for panel of plots, default is c(5,3) for 15 outcomes (default for 'outcomes' parameter)
#' @param linecol vector of line colors for plot. Needs to have same length as ph1_sim_OC_output
#' @param title add title to plot
#' @return plot
#' @export
#' @examples
#' \dontrun{
#' BOIN<-ph1_sim_OC(nsim=10, scenarios=list(c(0.40 ,0.50 ,0.60 ,0.70 ,0.80 ,0.90),
#'                                          c(0.30 ,0.40 ,0.50 ,0.60 ,0.70 ,0.80),
#'                                          c(0.30 ,0.35 ,0.40 ,0.60 ,0.75 ,0.90),
#'                                          c(0.20 ,0.30 ,0.40 ,0.50 ,0.60 ,0.70),
#'                                          c(0.20 ,0.30 ,0.50 ,0.65 ,0.80 ,0.90),
#'                                          c(0.20 ,0.25 ,0.30 ,0.40 ,0.60 ,0.75),
#'                                          c(0.10 ,0.20 ,0.30 ,0.40 ,0.50 ,0.60),
#'                                          c(0.10 ,0.15 ,0.25 ,0.30 ,0.50 ,0.75),
#'                                          c(0.05 ,0.10 ,0.15 ,0.20 ,0.25 ,0.30),
#'                                          c(0.05 ,0.10 ,0.20 ,0.30 ,0.50 ,0.65),
#'                                          c(0.025,0.05 ,0.10 ,0.20 ,0.30 ,0.40)), 
#'                  phi=0.3, phi1=0.6*0.3, phi2=1.4*0.3, maxtox=0.95, N=18, 
#'                  cohortsize=3, maxNretain=9, acc_tit=0, design="BOIN", MTD_safer=TRUE)
#' 
#' three<-ph1_sim_OC(nsim=10, scenarios=list(c(0.40 ,0.50 ,0.60 ,0.70 ,0.80 ,0.90),
#'                                           c(0.30 ,0.40 ,0.50 ,0.60 ,0.70 ,0.80),
#'                                           c(0.30 ,0.35 ,0.40 ,0.60 ,0.75 ,0.90),
#'                                           c(0.20 ,0.30 ,0.40 ,0.50 ,0.60 ,0.70),
#'                                           c(0.20 ,0.30 ,0.50 ,0.65 ,0.80 ,0.90),
#'                                           c(0.20 ,0.25 ,0.30 ,0.40 ,0.60 ,0.75),
#'                                           c(0.10 ,0.20 ,0.30 ,0.40 ,0.50 ,0.60),
#'                                           c(0.10 ,0.15 ,0.25 ,0.30 ,0.50 ,0.75),
#'                                           c(0.05 ,0.10 ,0.15 ,0.20 ,0.25 ,0.30),
#'                                           c(0.05 ,0.10 ,0.20 ,0.30 ,0.50 ,0.65),
#'                                           c(0.025,0.05 ,0.10 ,0.20 ,0.30 ,0.40)), 
#'                phi=0.3, design="3+3",acc_tit=1)
#' ph1_sim_OC_plot(ph1_sim_OC_output=list(BOIN,three), 
#' ylim_N=36,save_as="C:/Users/kdhollander/Desktop/test", outcomes=c(3:17), 
#' linecol=c("blue","red"),title="Test")
#' }

ph1_sim_OC_plot<- function(ph1_sim_OC_output, ylim_N, save_as, outcomes=c(3:17), layout=c(5,3),linecol, title){
  
  png(filename=paste0(save_as,"_",Sys.Date(),".png"),width=1000,height=1000)
  
    par(oma = c(1,1,1,1))
    par(mfrow=layout, mar=c(1.5,1.5,4,1),lwd=2, cex.axis=1.5, cex.main=2)
    mtext(title, outer=TRUE, cex = 1.5)
    
    ds         <- ph1_sim_OC_output
    N_design   <- length(ds)
    N_scen     <- length(ds[[1]])
    N_outcomes <- length(outcomes)
    
    for (i in 1:N_outcomes){
      var<-outcomes[i]
      
      for (j in 1:N_design){
        
        y_results <- sapply(ds[[j]],"[[",2)[var,]
        y_labs    <- sapply(ds[[j]],"[[",3)[var,1]
        
        varname   <- names(sapply(ds[[j]],"[[",2)[var,1])
        
        if (j==1){
          
          if (varname  %in% c("avg_SS","max_SS","avg_npt_MTD_sel","avg_nDLT")){
            plot (1:N_scen,y_results,type="b",col=linecol[j], xaxt="n",xlab="",ylab="",main=paste0(i,")",y_labs),ylim=c(0,ylim_N))
          } else {
            plot (1:N_scen,y_results,type="b",col=linecol[j], xaxt="n",xlab="",ylab="",main=paste0(i,")",y_labs),ylim=c(0,1))
          }
          
          axis(1,at = seq(1,N_scen, by = 1), las=1)
        }
        
        if (j>1){
          lines(1:N_scen,y_results,type="b",col=linecol[j])
        }
        
        grid()
        
      }
    }
  dev.off()
}

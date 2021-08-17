

#' @title Plot scenarios of true DLT rates
#' @description Plot scenarios of true DLT rates
#' @param scenarios scenarios of true DLT rates for each dose level. This also defines the number of dose levels
#' @param save_as where to save output
#' @return plot (red: first design; blue: second design; green: third design)
#' @export
#' @examples
#' plot_scenarios(scenarios=list(c(0.40 ,0.50 ,0.60 ,0.70 ,0.80 ,0.90),
#'                               c(0.30 ,0.40 ,0.50 ,0.60 ,0.70 ,0.80)),
#'                               save_as="C:/Users/kdhollander/Desktop/test")

plot_scenarios<-function(scenarios,save_as){
  png(filename=paste0(save_as,"_",Sys.Date(),".png"),width=600,height=400)
  par(mar=c(4,4,0,2),lwd=1.2, cex.axis=1, cex.main=1)
  plot(1:length(scenarios[[1]]),scenarios[[1]],type="b",ylim=c(0,1),col=1,ylab="True DLT rate",xlab="Dose level")
  for (i in 2:length(scenarios)){
    lines(1:length(scenarios[[1]]),scenarios[[i]],type="b",col=i)
  }
  abline(h=0.3,lty=2,col="red",lwd=2)
  dev.off()
}
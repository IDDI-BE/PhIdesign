% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ph1_sim_OC_plot.R
\name{ph1_sim_OC_plot}
\alias{ph1_sim_OC_plot}
\title{Plot operating characteristics of phase 1 design}
\usage{
ph1_sim_OC_plot(
  ph1_sim_OC_output,
  ylim_N,
  save_as,
  outcomes = c(3:17),
  layout = c(5, 3),
  linecol,
  title
)
}
\arguments{
\item{ph1_sim_OC_output}{list of output objects of \code{\link{ph1_sim_OC}}}

\item{ylim_N}{max N for plots with number of patients (typically max sample size of 3+3 design)}

\item{save_as}{where to save output}

\item{outcomes}{select outcomes from ph1_sim_OC_output (item 2)}

\item{layout}{layout for panel of plots, default is c(5,3) for 15 outcomes (default for 'outcomes' parameter)}

\item{linecol}{vector of line colors for plot. Needs to have same length as ph1_sim_OC_output}

\item{title}{add title to plot}
}
\value{
plot
}
\description{
Plot operating characteristics of phase 1 design. First design is "blue", second one"red",
    and third one "darkgreen"
}
\examples{
\dontrun{
BOIN<-ph1_sim_OC(nsim=10, scenarios=list(c(0.40 ,0.50 ,0.60 ,0.70 ,0.80 ,0.90),
                                         c(0.30 ,0.40 ,0.50 ,0.60 ,0.70 ,0.80),
                                         c(0.30 ,0.35 ,0.40 ,0.60 ,0.75 ,0.90),
                                         c(0.20 ,0.30 ,0.40 ,0.50 ,0.60 ,0.70),
                                         c(0.20 ,0.30 ,0.50 ,0.65 ,0.80 ,0.90),
                                         c(0.20 ,0.25 ,0.30 ,0.40 ,0.60 ,0.75),
                                         c(0.10 ,0.20 ,0.30 ,0.40 ,0.50 ,0.60),
                                         c(0.10 ,0.15 ,0.25 ,0.30 ,0.50 ,0.75),
                                         c(0.05 ,0.10 ,0.15 ,0.20 ,0.25 ,0.30),
                                         c(0.05 ,0.10 ,0.20 ,0.30 ,0.50 ,0.65),
                                         c(0.025,0.05 ,0.10 ,0.20 ,0.30 ,0.40)), 
                 phi=0.3, phi1=0.6*0.3, phi2=1.4*0.3, maxtox=0.95, N=18, 
                 cohortsize=3, maxNretain=9, acc_tit=0, design="BOIN", MTD_safer=TRUE)

three<-ph1_sim_OC(nsim=10, scenarios=list(c(0.40 ,0.50 ,0.60 ,0.70 ,0.80 ,0.90),
                                          c(0.30 ,0.40 ,0.50 ,0.60 ,0.70 ,0.80),
                                          c(0.30 ,0.35 ,0.40 ,0.60 ,0.75 ,0.90),
                                          c(0.20 ,0.30 ,0.40 ,0.50 ,0.60 ,0.70),
                                          c(0.20 ,0.30 ,0.50 ,0.65 ,0.80 ,0.90),
                                          c(0.20 ,0.25 ,0.30 ,0.40 ,0.60 ,0.75),
                                          c(0.10 ,0.20 ,0.30 ,0.40 ,0.50 ,0.60),
                                          c(0.10 ,0.15 ,0.25 ,0.30 ,0.50 ,0.75),
                                          c(0.05 ,0.10 ,0.15 ,0.20 ,0.25 ,0.30),
                                          c(0.05 ,0.10 ,0.20 ,0.30 ,0.50 ,0.65),
                                          c(0.025,0.05 ,0.10 ,0.20 ,0.30 ,0.40)), 
               phi=0.3, design="3+3",acc_tit=1)
ph1_sim_OC_plot(ph1_sim_OC_output=list(BOIN,three), 
ylim_N=36,save_as="C:/Users/kdhollander/Desktop/test", outcomes=c(3:17), 
linecol=c("blue","red"),title="Test")
}
}

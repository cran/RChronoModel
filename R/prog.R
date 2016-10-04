
#library(KernSmooth)
#library(hdrcde)

#####################################################
#     Estimation of Credible interval        #
#####################################################

#' Importing a CSV file 
#'
#' Importing a CSV file containing the output of the MCMC algorithm from any software
#'
#' @details 
#' @param file the name of the CSV file containing the output of the MCMC algorithm 
#' @param dec the character used in the file for decimal points for the use of read.csv()
#' @param sep the field separator character for the use of read.csv()
#' @param comment.char a character vector of length one containing a single character or an empty string for the use of read.csv()
#' @param header a character vector of length one containing a single character or an empty string for the use of read.csv()
#' @return A data frame (data.frame) containing a representation of the data in the file.
#' @export
#'
ImportCSV <- function(file, dec='.', sep=',', comment.char = '#', header = TRUE){
  
  # importing the CSV file
  data = read.csv(file, dec = dec, sep=sep, comment.char = comment.char, header=header)
  
}

#####################################################
#     Anteriority / posteriority Probability        #
#####################################################

#' Bayesian test for anteriority / posteriority
#'
#' Tests that a_chain < b_chain
#'
#' @param a_chain : numeric vector containing the output of the MCMC algorithm for the parameter a
#' @param b_chain : numeric vector containing the output of the MCMC algorithm for the parameter b
#' @return The baysesian probability that a < b knowing the data
#' @export
MarginalProba <- function(a_chain,b_chain){

  mean(ifelse(a_chain < b_chain, 1, 0))   ## bayesion test : a < b

}


#####################################################
#     Estimation of Credible interval        #
#####################################################

#' Bayesian credible interval
#'
#' Estimation of the shorest credible interval of the output of the MCMC algorithm for the parameter a
#'
#' @details A 100*level % credible interval is an interval that keeps N*(1-level) elements of the sample outside the interval
#' The 100*level % credible interval is the shortest of all those intervals.
#' @param a_chain numeric vector containing the output of the MCMC algorithm for the parameter a
#' @param level probability corresponding to the level of confidence used for the credible interval and the highest density region
#' @return The endpoints of the shortest credible interval
#' @export
#'
CredibleInterval <- function(a_chain, level=0.95){

  sorted_sample <- sort(a_chain)     # ordering the sample
  N = length(a_chain)                # calculation of the sample size
  OutSample = N * (1-level)          # calculation of the number of data to be outside the interval

  I =  cbind(sorted_sample[1:(OutSample+1)] , sorted_sample[(N-OutSample):N])    #   combinasion of all credible intervals

  l = I[,2]-I[,1]   # length of intervals
  i <- which.min(l) # look for the shortest interval

  c(level = level, CredibleIntervalInf=I[i,1],CredibleIntervalSup=I[i,2])   # returns the level and the endpoints

}




#####################################################
#            Marginal  Statistics                   #
#####################################################
#' Summary statistics
#'
#' Estimation of all usual statistics
#'
#' @param a_chain numeric vector containing the output of the MCMC algorithm for the parameter a
#' @param level probability corresponding to the level of confidence used for the credible interval and the highest density region
#' @param max_decimal maximum number of decimal
#' @param title title of the summary statistics
#' @return A list of values corresponding to all the following statistics
#' @export
MarginalStatistics <- function(a_chain, level=0.95, max_decimal=0){

  # Position
  mean = round(mean(a_chain), max_decimal)
  hdr = hdr(a_chain, prob = c(level * 100))
  map = round(hdr$mode, max_decimal)
  quantiles = round(quantile(a_chain, c(0.25,0.5,0.75)), max_decimal)

  # Dispersion
  sd = round(sd(a_chain), max_decimal)            # standard deviation using the 'sd' function
  CI = c(round(CredibleInterval(a_chain, level)[2], max_decimal), round(CredibleInterval(a_chain, level)[3], max_decimal))           # Credible Interval using the function 'CredibleInterval' from the package 'Rchronomodel'
  HPDR = round(hdr$hdr, max_decimal)              # Highest posterior density function region using the function 'hdr' from the package 'hdrcde'

  # Resulted
  res = c(mean, map, sd, quantiles[1], quantiles[2], quantiles[3], level, CI[1], CI[2], HPDR)
  Mat = matrix(nrow=length(res), ncol=1)
  Mat[,1] = res
  
  nom=c()
  for( k in (1: (length(HPDR)/2) ) ) {
    nom=c(nom,paste("HPDRInf",k))
    nom=c(nom,paste("HPRDSup",k))
  }
  names1 = c("mean", "MAP", "sd", "Q1", "median", "Q2", "level", "CredibleInterval Inf", "CredibleInterval Sup")
  rownames(Mat) = c(names1, nom)
  return(Mat)
}




#####################################################
#          Marginal posterior Density               #
#####################################################
#' Marginal posterior Density
#'
#' Plots the density of a_chain + statistics (mean, CI, HPDR)
#'
#' @param a_chain numeric vector containing the output of the MCMC algorithm for the parameter a
#' @param level probability corresponding to the level of confidence used for the credible interval and the highest density region
#' @param title label of the title
#' @param colors if TRUE  -> use of colors in the graph
#' @param GridLength length of the grid used to estimate the density
#' @return a plot with the density of a_chain + +CI + mean + HDR
#' @export
MarginalPlot <- function(a_chain, level=0.95, title="Marginal posterior density", colors=T, GridLength=1024){

  step <- max(density(a_chain, n=GridLength)$y) /50   # used to draw CI and mean above the curve
  maxValuex <- max(density(a_chain, n=GridLength)$x)
  minValuex <- min(density(a_chain, n=GridLength)$x)
  middleValuex <- minValuex + ( maxValuex - minValuex ) / 2
  P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
  P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4
  maxValuey <- max(density(a_chain, n=GridLength)$y)
  middleValuey <- maxValuey /2

  if (colors==T){
    par(mfrow=c(1,1))
    plot(density(a_chain, n=GridLength), main = title, xlab = "Date", axes = F, ylim=c(0,max(density(a_chain, n=GridLength)$y) + step))
    # abscissa axis
    axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor(minValuex), floor(P1Valuex), floor(middleValuex), floor(P3Valuex), floor(maxValuex )))
    # ordinate axis
    axis(2, at=c(0, middleValuey , maxValuey),labels =c(0, round(middleValuey, 3), round(maxValuey, 3)) )

    segments(CredibleInterval(a_chain, level)[2], 0, CredibleInterval(a_chain, level)[3], 0, lwd=6, col = 4)
    points(mean(a_chain), 0 , lwd=6, col = 2)

    # legend
    legend(P3Valuex, maxValuey, c("Density", "Credible Interval", "Mean"), lty=c(1, 1, 0), bty="n", pch=c(NA, NA,1), col = c("black","blue","red"), lwd=c(1,6,6), x.intersp=0.5, cex=0.9)

  }else {
    par(mfrow=c(1,1))
    plot(density(a_chain, n=GridLength), main = title, xlab = "Date", axes = F, ylim=c(0,max(density(a_chain, n=GridLength)$y) + step))
    # abscissa axis
    axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor(minValuex), floor(P1Valuex), floor(middleValuex), floor(P3Valuex), floor(maxValuex )))
    # ordinate axis
    axis(2, at=c(0, middleValuey , maxValuey),labels =c(0, round(middleValuey, 3), round(maxValuey, 3)) )

    segments(CredibleInterval(a_chain, level)[2], 0, CredibleInterval(a_chain, level)[3], 0, lwd=6, lty=1)
    points(mean(a_chain), 0 , lwd=6, pch=1)

    # legend
    legend(P3Valuex, maxValuey, c("Density", "Credible Interval", "Mean"), lty=c(1,1,0), pch=c(NA,NA,1), bty="n", lwd=c(1,6,6), x.intersp=0.5, cex=0.9)
  }

}

#####################################################
#          Hiatus between two dates                 #
#####################################################
#'  Test of the hiatus between two dates
#' Finds if it exists a gap between two dates that is the longest interval that satisfies : P(a_chain < IntervalInf < IntervalSup < b_chain | M) = level
#'
#' @param a_chain : numeric vector containing the output of the MCMC algorithm for the parameter a
#' @param b_chain : numeric vector containing the output of the MCMC algorithm for the parameter b
#' @param level probability corresponding to the level of confidence
#' @return The endpoints of the longest gap
#' @export
DatesHiatus <- function(a_chain, b_chain, level=0.95){

  if(length(a_chain) != length(b_chain)) {stop('Error : the parameters do not have the same length')} # test for the length of both chains
   
       gamma = mean((a_chain<b_chain))
       if (gamma < level) {print("No hiatus at this level")
       					   return(c(level=level, HiatusIntervalInf='NA',HiatusIntervalSup='NA')) } else # 
       {
    
      interval <- function(epsilon, P1End, P2Beginning, level)
      {
        q1 = quantile(P1End ,probs = 1-epsilon) ;
        indz = (P1End < q1)
        q2 = quantile(P2Beginning[indz],probs= (1-level-epsilon)/(1-epsilon))
        c(q1,q2)
      }
      hia = Vectorize(interval,"epsilon")

      indz = which(a_chain<b_chain)
      epsilon = seq(0,1-level,gamma)
      p = hia(epsilon, a_chain[indz], b_chain[indz], level/gamma)
      rownames(p)<- c("HiatusIntervalInf", "HiatusIntervalSup")

      D<- p[2,]-p[1,]
      DD = D[D>0]

      if (length(DD) > 0){
        I = which(D==max(DD))
        interval2 = round( p[,I], 0)
        if (p[2,I] != p[1,I]) {
          c(level=level, interval2[1], interval2[2])
        } else {
          c(level=level, HiatusIntervalInf='NA',HiatusIntervalSup='NA')
        }#end if (p[2,I] != p[1,I])

      } else {
        c(level=level, HiatusIntervalInf='NA',HiatusIntervalSup='NA')
        }#end if (length(DD) > 0)

      } # end if( sum(ifelse(PhaseBeginning < PhaseEnd, 1, 0) == length(PhaseBeginning) ) ) {  # test for Beginning < End


}




#####################################################
#                 Phase Time Range                 #
#####################################################

#' Phase Time Range
#'
#' Computes the shortest interval that satisfies : P(Phase1End < IntervalInf < IntervalSup < Phase2Beginning | M) = level
#'
#' @param PhaseBeginning_chain : numeric vector containing the output of the MCMC algorithm for the beginning of the phase
#' @param PhaseEnd_chain : numeric vector containing the output of the MCMC algorithm for the end of the phase
#' @param level probability corresponding to the desired level of confidence
#' @param max_decimal maximum number of decimal
#' @return The endpoints of the shortest time range associated with the desired level
#' @export
PhaseTimeRange <- function(PhaseBeginning_chain, PhaseEnd_chain, level=0.95, max_decimal=2){

  if(length(PhaseEnd_chain) != length(PhaseBeginning_chain)) { print('Error : the parameters do not have the same length')}   # test the length of both chains
  else{

    if( sum(ifelse(PhaseBeginning_chain <= PhaseEnd_chain, 1, 0)) != length(PhaseBeginning_chain) )  {  # test for Beginning < End
      print('Error : PhaseBeginning should be older than PhaseEnd')
    } else {

      periode <- function(epsilon, PBeginning, PEnd, level){
        q1 = quantile(PBeginning, probs = epsilon)    # Computes the 'level'th quantile of the beginning of the phase
        indz = (PBeginning > q1)
        q2 = quantile(PEnd[indz], probs= (level/(1-epsilon)))
        c(q1,q2)
      }   # end periode <- function(epsilon, PBeginning, PEnd, level){
      per = Vectorize(periode,"epsilon")

      epsilon = seq(0,1-level,.001)       # sequence of values used to compute
      p = per(epsilon,PhaseBeginning_chain, PhaseEnd_chain, level)
      rownames(p)<- c("TimeRangeInf", "TimeRangeSup")

      D<- p[2,]-p[1,]     # computes the length of all intervals
      I = which.min(D)    # finds the shortest interval
      range = round(p[,I], max_decimal)
      c(level=level, range[1], range[2]) # returns the endpoints of the shortest interval

    }
    # end if( sum(ifelse(PhaseBeginning < PhaseEnd, 1, 0) == length(PhaseBeginning) ) ) {  # test for Beginning < End

  }# end if(length(PhaseEnd) != length(PhaseBeginning)) {

} # end PhaseTimeRange <- function(PhaseBeginning, PhaseEnd, level, plot = F){





#####################################################
#             Statistics  for Phases               #
#####################################################

#' Summary statistics for phases
#'
#' Estimation of all usual statistics of the beginning and the end of a phase and the duration of the phase
#'
#' @param PhaseBeginning_chain numeric vector containing the output of the MCMC algorithm for the beginning of the phase
#' @param PhaseEnd_chain numeric vector containing the output of the MCMC algorithm for the end of the phase
#' @param level probability corresponding to the level of confidence used for the credible interval and the highest density region
#' @param max_decimal maximum number of decimal
#' @return A matrix of values corresponding to all the summary statistics
#' @export
PhaseStatistics <- function(PhaseBeginning_chain, PhaseEnd_chain, level=0.95, max_decimal=0){

  #Statistics according to PhaseBeginning_chain
  BeginningStat = MarginalStatistics(PhaseBeginning_chain, level, max_decimal)

  #Statistics according to PhaseEnd_chain parameter
  EndStat = MarginalStatistics(PhaseEnd_chain, level, max_decimal)

  #Statistics according to the duration of the phase
  DurationStat = MarginalStatistics(PhaseEnd_chain-PhaseBeginning_chain, level, max_decimal)

  # Resulted List
  
  if (length(BeginningStat) > length(EndStat)) {
    
    NbDiff = length(BeginningStat) - length(EndStat)
    Add = rep(NA,NbDiff) 
    EndStat = c(EndStat, Add)

  }else if (length(BeginningStat) < length(EndStat)) {
    
    NbDiff = length(EndStat) - length(BeginningStat)
    Add = rep(NA,NbDiff) 
    BeginningStat = c(BeginningStat, Add)
    
  }
  Mat1 = cbind(BeginningStat, EndStat)
  
  if (dim(Mat1)[1] > length(DurationStat)) {
    
    NbDiff = dim(Mat1)[1] - length(DurationStat)
    Add = rep(NA,NbDiff) 
    DurationStat = c(DurationStat, Add)
    
  }
  
  Mat = cbind(Mat1, DurationStat) 
  colnames(Mat) = c("Beginning", "End", "Duration")
  return(Mat)

}


#####################################################
#           Phase marginal density plot             #
#####################################################

#' Phase marginal density plot
#'
#' Plot of the density of Beginning and the End of a phase and summary statistics (mean, CI)
#'
#' @param PhaseBeginning_chain numeric vector containing the output of the MCMC algorithm for the beginning of the phase
#' @param PhaseEnd_chain numeric vector containing the output of the MCMC algorithm for the end of the phase
#' @param level probability corresponding to the level of confidence used for the credible interval and the time range
#' @param title The Title of the graph
#' @param colors if TRUE  -> use of colors in the graph
#' @param GridLength length of the grid used to estimate the density
#' @return A plot with the density of PhaseBeginning_chain + PhaseEnd_chain + additionnal summary statitsics
#' @export

PhasePlot <- function(PhaseBeginning_chain, PhaseEnd_chain, level=0.95, title = "Phase marginal posterior densities", colors = T, GridLength=1024){

  if(length(PhaseEnd_chain) != length(PhaseBeginning_chain)) { print('Error : the parameters do not have the same length')}   # test the length of both chains
  else{

  if( sum(ifelse(PhaseBeginning_chain <= PhaseEnd_chain, 1, 0)) == length(PhaseBeginning_chain) ) {

    minValuex <- min(density(PhaseBeginning_chain, n=GridLength)$x)
    maxValuex <- max(density(PhaseEnd_chain, n=GridLength)$x)
    middleValuex <- ( maxValuex + minValuex) / 2
    P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
    P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4


    maxValuey <- max ( max(density(PhaseBeginning_chain, n=GridLength)$y) , max(density(PhaseEnd_chain, n=GridLength)$y))
    middleValuey <- maxValuey /2
    step <- maxValuey  /20

    if (colors==T){
    # first graph
    par(las=1, mfrow=c(1,1), cex.axis=0.8)
    plot(density(PhaseEnd_chain, n=GridLength), main = title, xlab="Date", axes = F, ylim=c(0,maxValuey+step), xlim=c(minValuex, maxValuex), lty =1, lwd=2, col="steelblue4")
    lines(density(PhaseBeginning_chain, n=GridLength), lty =1, lwd=2, col ="steelblue1")

    # abscissa axis
    axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor(minValuex), floor(P1Valuex), floor(middleValuex), floor(P3Valuex), floor(maxValuex)))
    # ordinate axis
    axis(2, at=c(0, middleValuey, maxValuey), labels =c(0, round(middleValuey, 3), round(maxValuey, 3)) )

    # segment representing the CredibleInterval of the end of the phase
    CIEnd = CredibleInterval(PhaseEnd_chain, level)
    segments(CIEnd[2], 0, CIEnd[3], 0, lty = 1, lwd=6, col = "steelblue4")
    # point in red representing the mean
    points(mean(PhaseEnd_chain), 0, lwd=6, col = "steelblue4")
    # segment representing the Time range of the phase in green
    PTR = PhaseTimeRange(PhaseBeginning_chain, PhaseEnd_chain, level)
    segments(PTR[2], maxValuey+step, PTR[3], maxValuey+step, lwd=6, col = "violetred4")

    # segment representing the CredibleInterval in blue
    CIBeginning = CredibleInterval(PhaseBeginning_chain, level)
    segments(CIBeginning[2], step, CIBeginning[3], step, lty= 1, lwd=6, col = "steelblue1")
    # point in red representing the mean
    points(mean(PhaseBeginning_chain), step , lwd=6, col = "steelblue1")

    # legend
    legend(P3Valuex, maxValuey, c("Density of the Beginning", "Density of the End", "with Credible Interval" ," and Mean (o)", " Phase Time Range"), lty=c(1,1,0,0,1), bty="n",col = c("steelblue1","steelblue4","black","black","violetred4"), lwd=c(2,2,6,6,6), x.intersp=0.5, cex=0.9)

    } else {

      # first graph
      par(las=1, mfrow=c(1,1), cex.axis=0.8)
      plot(density(PhaseEnd_chain, n=GridLength), main = title, xlab="Date", axes = F, ylim=c(0,maxValuey+step), xlim=c(minValuex, maxValuex), lty =2, lwd=2)
      lines(density(PhaseBeginning_chain, n=GridLength), lty =3, lwd=2)

      # abscissa axis
      axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor(minValuex), floor(P1Valuex), floor(middleValuex), floor(P3Valuex), floor(maxValuex)))
      # ordinate axis
      axis(2, at=c(0, middleValuey, maxValuey), labels =c(0, round(middleValuey, 3), round(maxValuey, 3)) )

      # segment representing the CredibleInterval in blue
      CIEnd = CredibleInterval(PhaseEnd_chain, level)
      segments(CIEnd[2], 0, CIEnd[3], 0, lty = 2, lwd=6, col = 1)
      # point in red representing the mean
      points(mean(PhaseEnd_chain), 0, lwd=6, col = 1)
      # segment representing the Time range of the phase in green
      PTR = PhaseTimeRange(PhaseBeginning_chain, PhaseEnd_chain, level)
      segments(PTR[2], maxValuey+step, PTR[3], maxValuey+step, lwd=6, col = 1)

      # segment representing the CredibleInterval in blue
      CIBeginning = CredibleInterval(PhaseBeginning_chain, level)
      segments(CIBeginning[2], step, CIBeginning[3], step, lty= 3, lwd=6, col = 1)
      # point in red representing the mean
      points(mean(PhaseBeginning_chain), step , lwd=6, col = 1)

      # legend
      legend(P3Valuex, maxValuey, c("Density of the Beginning", "Density of the End", "with Credible Interval", "and Mean (o)", " Phase Time Range"), lty=c(3,2,0,0,1), bty="n",col = c(1,1,1,1,1), lwd=c(2,2,6,6,6), x.intersp=0.5, cex=0.9)

    }

  } else {
    print('Error : PhaseBeginning_chain should be older than PhaseEnd_chain')
  }

  } # end if(length(PhaseEnd_chain) != length(PhaseBeginning_chain))

}


#####################################################
#     Phase duration marginal density plot          #
#####################################################

#' Phase duration marginal density plot 
#'
#' Plot of the density of the duration of a phase and summary statistics (mean, CI)
#'
#' @param PhaseBeginning_chain numeric vector containing the output of the MCMC algorithm for the beginning of the phase
#' @param PhaseEnd_chain numeric vector containing the output of the MCMC algorithm for the end of the phase
#' @param level probability corresponding to the level of confidence used for the credible interval and the time range
#' @param title The Title of the graph
#' @param colors if TRUE  -> use of colors in the graph
#' @param GridLength length of the grid used to estimate the density
#' @return A plot with the density of the duration of the phase + additionnal summary statitsics
#' @export
PhaseDurationPlot <- function(PhaseBeginning_chain, PhaseEnd_chain, level=0.95, title = "Duration of the phase", colors = T, GridLength=1024){
  
  if(length(PhaseEnd_chain) != length(PhaseBeginning_chain)) { print('Error : the parameters do not have the same length')}   # test the length of both chains
  else{
    
    if( sum(ifelse(PhaseBeginning_chain < PhaseEnd_chain, 1, 0)) == length(PhaseBeginning_chain) ) {
      
      par(las=1, mfrow=c(1,1), cex.axis=0.8)
      MarginalPlot(PhaseEnd_chain -PhaseBeginning_chain, level, title = title, colors = colors, GridLength=GridLength)
      
      
    } else {
      print('Error : PhaseBeginning_chain should be older than PhaseEnd_chain')
    }
    
  } # end if(length(PhaseEnd_chain) != length(PhaseBeginning_chain))
  
}


#####################################################
#          Hiatus between two phases             #
#####################################################
#'  Gap/Hiatus between two successive phases (for phases in temporal order constraint)
#'
#' Finds if it exists a gap between two phases that is the longest interval that satisfies : P(Phase1End < IntervalInf < IntervalSup < Phase2Beginning | M) = level
#'
#' @param Phase1End_chain numeric vector containing the output of the MCMC algorithm for the end of the oldest phase
#' @param Phase2Beginning_chain numeric vector containing the output of the MCMC algorithm for the beginning of the following phase
#' @param level probability corresponding to the level of confidence
#' @param max_decimal maximum number of decimal
#' @return The endpoints of the longest gap
#' @export
PhasesGap <- function(Phase1End_chain, Phase2Beginning_chain, level=0.95, max_decimal=0){

  if(length(Phase1End_chain) != length(Phase2Beginning_chain)) { stop('Error : the parameters do not have the same length')} # test for the length of both chains
    else{

      if( sum(ifelse(Phase1End_chain <=Phase2Beginning_chain, 1, 0)) != length(Phase2Beginning_chain) )  {  # test for Beginning < End
        stop('Error : Phase1End_chain should be older than Phase2Beginning_chain')
      } else {

      interval <- function(epsilon, P1End, P2Beginning, level)
      {
        q1 = quantile(P1End ,probs = 1-epsilon) ;
        indz = (P1End < q1)
        q2 = quantile(P2Beginning[indz],probs= (1-level-epsilon)/(1-epsilon))
        c(q1,q2)
      }
      hia = Vectorize(interval,"epsilon")

      epsilon = seq(0,1-level,.001)
      p = hia(epsilon, Phase1End_chain, Phase2Beginning_chain, level)
      rownames(p)<- c("HiatusIntervalInf", "HiatusIntervalSup")

      D<- p[2,]-p[1,]
      DD = D[D>0]

      if (length(DD) > 0){
        I = which(D==max(DD))
        interval2 = round( p[,I], max_decimal)
        if (p[2,I] != p[1,I]) {
          c(level=level, interval2[1], interval2[2])
        } else {
          c(level=level, HiatusIntervalInf='NA',HiatusIntervalSup='NA')
        }#end if (p[2,I] != p[1,I])

      } else {
        c(level=level, HiatusIntervalInf='NA',HiatusIntervalSup='NA')
        }#end if (length(DD) > 0)

      } # end if( sum(ifelse(PhaseBeginning < PhaseEnd, 1, 0) == length(PhaseBeginning) ) ) {  # test for Beginning < End

    } # if(length(Phase1End_chain) != length(Phase2Beginning_chain)) {

}





#####################################################
#                Phases Transition                  #
#####################################################

#'  Transition range between two successive phases (for phases in temporal order constraint)
#'
#' Finds if it exists the shortest interval that satisfies : P(TransitionRangeInf < Phase1End < Phase2Beginning < TransitionRangeSup  | M) = level
#'
#' @param Phase1End_chain numeric vector containing the output of the MCMC algorithm for the end of the oldest phase
#' @param Phase2Beginning_chain numeric vector containing the output of the MCMC algorithm for the beginning of the following phase
#' @param level probability corresponding to the level of confidence
#' @param max_decimal maximum number of decimal
#' @return the endpoints of the transition interval
#' @export
PhasesTransition <- function(Phase1End_chain, Phase2Beginning_chain, level=0.95, max_decimal=0){

  result = as.matrix( PhaseTimeRange(Phase1End_chain, Phase2Beginning_chain, level=level, max_decimal=max_decimal))
  rownames(result)<- c(level=level, "TransitionRangeInf", "TransitionRangeSup")
  result <- t(result)

  return(result[1,])
}




#####################################################
#             Succession Plot                  #
#####################################################

#' Density Plots of two successive phases (for phases in temporal order constraint)
#'
#' Plot of the densities of two successive phases + statistics (mean, CI, HPDR)
#'
#' @param Phase1Beginning_chain numeric vector containing the output of the MCMC algorithm for the beginning of the oldest phase
#' @param Phase1End_chain numeric vector containing the output of the MCMC algorithm for the end of the oldest phase
#' @param Phase2Beginning_chain numeric vector containing the output of the MCMC algorithm for the beginning of the following phase
#' @param Phase2End_chain numeric vector containing the output of the MCMC algorithm for the end of the following phase
#' @param level probability corresponding to the level of confidence
#' @param title title of the graph
#' @param GridLength length of the grid used to estimate the density
#' @return a plot of all densities + CI + mean + HDR
#' @export

SuccessionPlot <- function(Phase1Beginning_chain, Phase1End_chain, Phase2Beginning_chain, Phase2End_chain, level=0.95,  title = "Phases marginal posterior densities", GridLength=1024){


  if(length(Phase1End_chain) != length(Phase2Beginning_chain)) { stop('Error : the parameters do not have the same length')} # test for the length of both chains
  else{

  if( sum(ifelse(Phase1Beginning_chain <= Phase1End_chain, 1, 0)) != length(Phase1Beginning_chain) ||  sum(ifelse(Phase2Beginning_chain <= Phase2End_chain, 1, 0)) != length(Phase1Beginning_chain) || sum(ifelse( Phase1End_chain <= Phase2Beginning_chain, 1, 0)) != length(Phase1Beginning_chain) ) {
    # test for Beginning < End and Phase1 < Phase2
    stop('Error : PhaseBeginning_chain should be older than PhaseEnd_chain')
  } else {

    minValuex <- min(density(Phase1Beginning_chain, n=GridLength)$x, density(Phase2Beginning_chain, n=GridLength)$x)
    maxValuex <- max(density(Phase1End_chain, n=GridLength)$x, density(Phase2End_chain, n=GridLength)$x)
    middleValuex <- ( maxValuex + minValuex) / 2
    P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
    P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4

    maxValuey <- max ( max(density(Phase1Beginning_chain, n=GridLength)$y) , max(density(Phase1End_chain, n=GridLength)$y), max(density(Phase2Beginning_chain, n=GridLength)$y) , max(density(Phase2End_chain, n=GridLength)$y))
    middleValuey <- maxValuey /2
    minValuey <- min ( min(density(Phase1Beginning_chain, n=GridLength)$y) , min(density(Phase1End_chain, n=GridLength)$y), min(density(Phase2Beginning_chain, n=GridLength)$y) , min(density(Phase2End_chain, n=GridLength)$y))

    haut = seq(minValuey,maxValuey,length.out=5)
    middleA <- maxValuey+ (haut[1] + haut[2]) / 2

    plot(density(Phase1End_chain, n=GridLength), main = title, ylab="Density", xlab = "Date", ylim=c(0,maxValuey+maxValuey), xlim=c(minValuex, maxValuex), bty='n',lty =1, lwd=2, axes=F, col = "steelblue")
    lines(density(Phase1Beginning_chain, n=GridLength), lty =1, lwd=2, col = "steelblue")

    # abscissa axis
    axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor( minValuex), floor( P1Valuex), floor( middleValuex), floor( P3Valuex), floor( maxValuex)))
    # ordinate axis
    axis(2, at=c(0, middleValuey, maxValuey), labels =c(0, round(middleValuey, 3), round(maxValuey, 3)) )

    ## Phase2
    lines(density(Phase2Beginning_chain, n=GridLength), lty =1, lwd=2, col ="violet")
    lines(density(Phase2End_chain, n=GridLength), lty =1, lwd=2, col ="violet")

    ## Phase Time Range
    PTR1 = PhaseTimeRange(Phase1Beginning_chain, Phase1End_chain, level=level)
    PTR2 = PhaseTimeRange(Phase2Beginning_chain, Phase2End_chain, level=level)
    segments(PTR1[2],maxValuey+haut[2],PTR1[3],maxValuey+haut[2],lwd=6,col="steelblue")
    segments(PTR2[2],maxValuey+haut[3],PTR2[3],maxValuey+haut[3],lwd=6, col ="violet")
    text(minValuex, middleA,"Time range",srt =90)

    ## Phase Transition
    PTrans = PhasesTransition(Phase1End_chain, Phase2Beginning_chain, level=level)
    segments(PTrans[2],maxValuey+haut[4],PTrans[3],maxValuey+haut[4],lwd=6, col = "steelblue")
    segments(PTrans[2],maxValuey+haut[4],PTrans[3],maxValuey+haut[4],lwd=6, col = "violet", lty=4)

    PGap = PhasesGap(Phase1End_chain, Phase2Beginning_chain, level=level)
    if (PGap[2] == "NA" || PGap[3] == "NA") {
      points( (PTrans[3]+PTrans[2])/2, maxValuey+haut[5], lwd=2, col = "steelblue", pch=4)
    } else {
      segments(PGap[2],maxValuey+haut[5],PGap[3],maxValuey+haut[5],lwd=6, col = "steelblue")
      segments(PGap[2],maxValuey+haut[5],PGap[3],maxValuey+haut[5],lwd=6, col = "violet", lty=4)
    }

    text(minValuex, maxValuey+haut[4],"Transition",srt =90)
    text(minValuex, maxValuey+haut[5],"Gap",srt =90)

  }

  }

}







#####################################################
#     Estimation of Credible interval        #
#####################################################

#' Bayesian credible interval for a series of MCMC chains
#'
#' Estimation of the shorest credible interval of the output of the MCMC algorithm for the parameter a
#'
#' @details A 100*level % credible interval is an interval that keeps N*(1-level) elements of the sample outside the interval
#' The 100*level % credible interval is the shortest of all those intervals.
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position numeric vector containing the position of the column corresponding to the MCMC chains of interest
#' @param level probability corresponding to the level of confidence used for the credible interval
#' @return The endpoints of the shortest credible interval
#' @export
#'
MultiCredibleInterval <- function(data, position, level=0.95){

  # number of chains
  L = length(position)

  # matrix of results for each pair of phases
  result = matrix(nrow=L, ncol=3)

  colnames(result) <- c("Level","CredibleIntervalInf", "CredibleIntervalSup")
  
  # names
  rownames(result) <- names(data)[position]

  for (i in 1:L) {

    sorted_sample <- sort(data[,position[i]])     # ordering the sample
    N = length(sorted_sample)                     # calculation of the sample size of the chain
    OutSample = N * (1-level)          # calculation of the number of data to be outside the interval

    I =  cbind(sorted_sample[1:(OutSample+1)] , sorted_sample[(N-OutSample):N])    #   combinasion of all credible intervals

    l = I[,2]-I[,1]   # length of intervals
    j <- which.min(l) # look for the shortest interval

    result[i,] =   c(level, round(I[j,1], digits = 0) , round(I[j,2],digits = 0) )   # returns the level and the endpoints

  }
  return(result)
}





#####################################################
#                   MultiHPD                        #
#####################################################

#' Bayesian HPD regions for a series of MCMC chains
#'
#' Estimation of the HPD region of the output of the MCMC algorithm for the parameter a
#'
#' @details Highest posterior density function region using the function 'hdr' from the package 'hdrcde'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position numeric vector containing the position of the column corresponding to the MCMC chains of interest
#' @param level probability corresponding to the level of confidence used for the credible interval
#' @return The endpoints of the shortest credible interval
#' @export
#'
MultiHPD <- function(data, position, level=0.95){
  
  # matrix of results for each pair of phases
    hdr = hdr(data[,position[1]], prob = c(level * 100))$hdr
    HPDR = round(hdr, digits = 0)  
    result = matrix( c(level, HPDR), nrow=1)
    dim = dim(result)[2]
    
    if(length(position) > 1){
      
      for (i in 2:length(position)) {
        
        hdr = hdr(data[,position[i]], prob = c(level * 100))$hdr
        HPDR = round(hdr, digits = 0) 
        res = c(level, HPDR)
        
        if (length(res) > dim) {
          NbCol = length(res) - dim
          AddColum = rep(NA, i-1) 
          Ajout = NULL 
          for (j in 1:NbCol){
            Ajout = cbind(Ajout, AddColum)
          }
          resultTemp = cbind(result, Ajout)
          result =  rbind(resultTemp, res)
          
        }else if (length(res) < dim) {
          NbCol = dim - length(res) 
          Add = rep(NA, NbCol) 
          Ajout = c(res,Add)

          result = rbind(result, Ajout)
          
        }else{
          result =  rbind(result, res)  
        }
        dim = dim(result)[2] 
      }

    }

    nom=c()
    for( k in (1:((dim-1)/2)) ) {
      nom=c(nom,paste("HPDRInf",k))
      nom=c(nom,paste("HPRDSup",k))
    }
    colnames(result) <- c("Level", nom)
    rownames(result) <- names(data)[position]
    
  return(result) # returns a matrix with the level and the endpoints
}

###############################################
#             MultiDatesPlot                  #
###############################################

#' Plot of credible intervals or HPD regions of a series of dates
#'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position numeric vector containing the position of the column corresponding to the MCMC chains of interest
#' @param level probability corresponding to the level of confidence
#' @param title title of the graph
#' @return a plot of the endpoints of the credible intervals of a series of dates
#' @export

MultiDatesPlot <- function(data, position, level=0.95, intervals = c("CI", "HPD"), title = "Plot of intervals"){
  
  if(intervals =="CI"){
  Bornes = MultiCredibleInterval(data, position, level=level) 
  }else if(intervals =="HPD") {
  Bornes = MultiHPD(data, position, level=level) 
  }
  Ordered = Bornes[order(Bornes[,2]),]
  nbCol = dim(Ordered)[2]
  
  NbDates = length(position)
  DatesNames <- rownames(Ordered)
  
  minValuex <- floor( min(Ordered[,2], na.rm = TRUE)/100) * 100
  maxValuex <- ceiling( max(Ordered[,-c(1,2)],na.rm = TRUE)/100) * 100
  middleValuex <- ( maxValuex + minValuex) / 2
  P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
  P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4
  seq = seq(minValuex, maxValuex)
  
  par(mar=c(5,6,4,2))
  plot(0, main = title, ylab="", xlab = "Time", ylim=c(0, NbDates), xlim=c(minValuex, maxValuex), type="n", axes=F)

  # abscissa axis
  axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) ) 
  # ordinate axis
  axis(2, at=1:NbDates, labels =DatesNames, las =2)
  
  ## Phase Time Range
  for (i in 1:NbDates ) { for (j in seq(2,(nbCol-1), by = 2)) { segments(Ordered[i,j], i, Ordered[i,j+1], i, lwd=6) } }
 
}



#####################################################
#         Multiple Phase Time Range                 #
#####################################################

#' Phase Time Range for multiple phases
#'
#' Computes the shortest interval that satisfies : P(Phase1End < IntervalInf < IntervalSup < Phase2Beginning | M) = level
#'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position_beginning numeric vector containing the column number corresponding to the beginning of the phases
#' @param position_end numeric vector containing the column number corresponding to the end of the phases set in the same order as in position_beginning
#' @param level probability corresponding to the desired level of confidence
#' @param max_decimal maximum number of decimal
#' @return The endpoints of the shortest time range associated with the desired level
#' @export


MultiPhaseTimeRange <- function(data, position_beginning, position_end=position_beginning+1, level=0.95, max_decimal=0){
  
  if (length(position_beginning)!= length(position_end)) {
    print('Error : the position vectors do not have the same length')
    } else {

  # number of phases
  L = length(position_beginning)
  
  # names
  names_beginning <- names(data)[position_beginning]
  names_end <- names(data)[position_end]
  
  # Construction of a new dataset containing the columns corresponding to the phases of interest
  phase = matrix(ncol = L*2, nrow=nrow(data))
  for (i in 1:L) {
    phase[,2*i-1] = data[,position_beginning[i]]
    phase[,2*i] = data[,position_end[i]]
  }
  
  # matrix of results
  result = matrix(nrow=L, ncol=3)
  colnames(result)<- c("Level","TimeRangeInf", "TimeRangeSup")
  
  phasenames <- vector(length = L)
  for (i in 1:L) { phasenames[i] = paste(names_beginning[i], names_end[i] ) }
  rownames(result)<- phasenames
  
  for (i in 1:L){
    result[i,] = PhaseTimeRange(phase[,2*i-1], phase[,2*i], level=level, max_decimal=max_decimal)
  }
  
  return(result)
  
  
  }
}

#####################################################
#       Hiatus between a succession of  phases      #
#####################################################

#'  Gap/Hiatus between a succession of phases (for phases in temporal order constraint)
#'
#' Finds if it exists a gap between two phases that is the longest interval that satisfies : P(Phase1End < IntervalInf < IntervalSup < Phase2Beginning | M) = level
#'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position_beginning numeric vector containing the column number corresponding to the beginning of the phases
#' @param position_end numeric vector containing the column number corresponding to the end of the phases set in the same order as in position_beginning
#' @param level probability corresponding to the level of confidence
#' @param max_decimal maximum number of decimal
#' @return The endpoints of the longest gap
#' @export

MultiPhasesGap <- function(data, position_beginning, position_end = position_beginning+1, level=0.95, max_decimal=0){
  
  if (length(position_beginning)!= length(position_end)) {
    print('Error : the position vectors do not have the same length')
  } else {
    
  # number of phases
  L = length(position_beginning)
  
  #names
  names_beginning <- names(data)[position_beginning]
  names_end <- names(data)[position_end]
  
  # Construction of a new dataset containing the columns corresponding to the phases of interest
  phase = matrix(ncol = L*2, nrow=nrow(data))
  for (i in 1:L) {
    phase[,2*i-1] = data[,position_beginning[i]]
    phase[,2*i] = data[,position_end[i]]
  }
  
  # matrix of results
  result = matrix(nrow=L-1, ncol=3)
  colnames(result)<- c("Level","HiatusIntervalInf", "HiatusIntervalSup")
  
  phasenames <- vector(length = (L-1))
  for (i in 1:L-1) { phasenames[i] = paste(names_end[i], "&", names_beginning[i+1]) }
  rownames(result)<- phasenames
  
  for (i in 1:(L-1)){
    result[i,] = PhasesGap(phase[,2*i], phase[,2*i+1], level=level, max_decimal=max_decimal)
  }
  
  return(result)
  
  }
}



#####################################################
#         Multiple Phases Transition               #
#####################################################

#'  Transition range for a succession of phases (for phases in temporal order constraint)
#'
#' Finds if it exists the shortest interval that satisfies : P(TransitionRangeInf < Phase1End < Phase2Beginning < TransitionRangeSup  | M) = level
#'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position_beginning numeric vector containing the column number corresponding to the beginning of the phases
#' @param position_end numeric vector containing the column number corresponding to the end of the phases set in the same order as in position_beginning
#' @param level probability corresponding to the level of confidence
#' @param max_decimal maximum number of decimal
#' @return the endpoints of the transition interval
#' @export


MultiPhasesTransition <- function(data, position_beginning, position_end = position_beginning+1, level=0.95, max_decimal=0){

  if (length(position_beginning)!= length(position_end)) {
    print('Error : the position vectors do not have the same length')
  } else {
    
  # number of phases
  L = length(position_beginning)
  
  #names
  names_beginning <- names(data)[position_beginning]
  names_end <- names(data)[position_end]
  
  # Construction of a new dataset containing the columns corresponding to the phases of interest
  phase = matrix(ncol = L*2, nrow=nrow(data))
  for (i in 1:L) {
    phase[,2*i-1] = data[,position_beginning[i]]
    phase[,2*i] = data[,position_end[i]]
  }

  # matrix of results
  result = matrix(nrow=L-1, ncol=3)
  colnames(result)<- c(level, "TransitionRangeInf", "TransitionRangeSup")
  
  phasenames <- vector(length = L-1)
  for (i in 1:L-1) { phasenames[i] = paste(names_end[i], "&", names_beginning[i+1]) }
  rownames(result)<- phasenames

  for (i in 1:(L-1)){
    result[i,] = PhaseTimeRange(phase[,2*i], phase[,2*i+1], level=level, max_decimal=max_decimal)
  }

  return(result)

  }
}





#####################################################
#        Multiple Phases Density Plots              #
#####################################################

#' Successive Phases Density Plots (for phases in temporal order constraint)
#'
#' Plot of the densities of several successive phases + statistics (mean, CI, HPDR)
#'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position_beginning numeric vector containing the column number corresponding to the beginning of the phases
#' @param position_end numeric vector containing the column number corresponding to the end of the phases set in the same order as in position_beginning
#' @param level probability corresponding to the level of confidence
#' @param title title of the graph
#' @return a plot of all densities + CI + mean + HDR
#' @export

MultiSuccessionPlot <- function(data, position_beginning, position_end = position_beginning+1, level=0.95, title = "Phases marginal posterior densities"){
  
  if (length(position_beginning)!= length(position_end)) {
    print('Error : the position vectors do not have the same length')
  } else {
    
    # construction of the new dataset
    L = length(position_beginning)
    
    GridLength=1024
    phase = matrix(ncol = L*2, nrow=nrow(data))
    densityX = matrix(ncol = L*2, nrow=GridLength)
    densityY = matrix(ncol = L*2, nrow=GridLength)
    
    for (i in 1:L) {
      phase[,2*i-1] = data[,position_beginning[i]]
      phase[,2*i] = data[,position_end[i]]
      
      densityX[,2*i-1] = density(data[,position_beginning[i]])$x
      densityX[,2*i] = density(data[,position_end[i]])$x
      
      densityY[,2*i-1] = density(data[,position_beginning[i]])$y
      densityY[,2*i] = density(data[,position_end[i]])$y
    }
    
    minValuex <- min (apply(densityX,2,min))
    maxValuex <- max( apply(densityX,2,max) )
    middleValuex <- ( maxValuex + minValuex) / 2
    P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
    P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4
    
    maxValuey <- max( apply(densityY,2,max))
    middleValuey <- maxValuey /2
    minValuey <- min( apply(densityY,2,min))
    
    haut = seq(minValuey,maxValuey,length.out=(3*L+1) )
    
    pal = rainbow(L)
    colors = pal[sample(x=1:L,L)]
    
    plot(density(phase[,1], n=GridLength), main = title, ylab="Density", xlab = "Date", ylim=c(0,2*maxValuey), xlim=c(minValuex, maxValuex), bty='n',lty =1, lwd=2, axes=F, col = colors[1])
    lines(density(phase[,2], n=GridLength), lty =1, lwd=2, col = colors[1])
    
    # abscissa axis
    axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor( minValuex), floor( P1Valuex), floor( middleValuex), floor( P3Valuex), floor( maxValuex)))
    # ordinate axis
    axis(2, at=c(0, middleValuey, maxValuey), labels =c(0, round(middleValuey, 5), round(maxValuey, 5)) )
    
    
    ## Following phases
    for(i in 2:L) {
      lines(density(phase[,2*i-1]), col=colors[i], lwd=2, lty=1)
      lines(density(phase[,2*i]), col=colors[i], lwd=2, lty=1)
    }
    
    ## Phase Time Range
    MPTR = MultiPhaseTimeRange(data=data, position_beginning=position_beginning, position_end=position_end, level=level)
    for (i in 1:(L) ) { segments(MPTR[i,2], maxValuey + haut[i+1], MPTR[i,3], maxValuey + haut[i+1], lwd=6,col=colors[i]) }
    text(minValuex, maxValuey + haut[2] , "Time range",srt =90)
    
    ## Phase Transition / Gap
    segments(minValuex, maxValuey+haut[L+2],maxValuex, maxValuey+haut[L+2], lwd=.2)
    segments(minValuex, maxValuey+haut[2*L+2],maxValuex, maxValuey+haut[2*L+2], lwd=.2)
    text(minValuex, maxValuey+haut[L+3 + trunc((L-1)/2)],"Transition ",srt =90)
    text(minValuex, maxValuey+haut[2*L+3 + trunc((L-1)/2)]," Gap",srt =90)
    #mtext(outer=TRUE, side = 1, "Transition", srt =90)
    
    PTrans = MultiPhasesTransition(data=data, position_beginning=position_beginning, position_end=position_end, level=level)
    PGap = MultiPhasesGap(data=data, position_beginning=position_beginning, position_end=position_end, level=level)
    
    for (i in 1:(L-1) ) {
      segments(PTrans[i,2], maxValuey+haut[L+2+i],PTrans[i,3], maxValuey+haut[L+2+i],lwd=6, col = colors[i])
      segments(PTrans[i,2], maxValuey+haut[L+2+i],PTrans[i,3], maxValuey+haut[L+2+i],lwd=6, col = colors[i+1], lty=4)
      
      if (PGap[i,2] == "NA" || PGap[i,3] == "NA") {
        points( (PTrans[i,3]+PTrans[i,2])/2, maxValuey+haut[2*L+2+i], lwd=2, col = colors[i], pch=4)
        
      } else {
        segments(as.numeric(PGap[i,2]), maxValuey+haut[2*L+2+i], as.numeric(PGap[i,3]), maxValuey+haut[2*L+2+i], lwd=6, col = colors[i])
        segments(as.numeric(PGap[i,2]), maxValuey+haut[2*L+2+i], as.numeric(PGap[i,3]), maxValuey+haut[2*L+2+i], lwd=6, col = colors[i+1], lty=4)
      }
      
    }
    
  }
}




#####################################################
#           Multiple Phases     Plots            #
#####################################################

#' Several Phases Density Plots
#'
#' Plot of the densities of several phases + statistics (mean, CI, HPDR)
#'
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position_beginning numeric vector containing the column number corresponding to the beginning of the phases
#' @param position_end numeric vector containing the column number corresponding to the end of the phases set in the same order as in position_beginning
#' @param level probability corresponding to the level of confidence
#' @param title title of the graph
#' @return a plot of all densities + CI + mean + HDR
#' @export

MultiPhasePlot <- function(data, position_beginning, position_end = position_beginning+1, level=0.95, title = "Phases marginal posterior densities"){

  if (length(position_beginning)!= length(position_end)) {
    print('Error : the position vectors do not have the same length')
  } else {
    
  # construction of the new dataset
  L = length(position_beginning) 
  
  phase = matrix(ncol = L*2, nrow=nrow(data))
  GridLength = 1024
  densityX = matrix(ncol = L*2, nrow=GridLength)
  densityY = matrix(ncol = L*2, nrow=GridLength)

  for (i in 1:L) {
    phase[,2*i-1] = data[,position_beginning[i]]
    phase[,2*i] = data[,position_end[i]]
    
    densityX[,2*i-1] = density(data[,position_beginning[i]])$x
    densityX[,2*i] = density(data[,position_end[i]])$x
    
    densityY[,2*i-1] = density(data[,position_beginning[i]])$y
    densityY[,2*i] = density(data[,position_end[i]])$y
  }

  minValuex <- min (apply(densityX,2,min))
  maxValuex <- max( apply(densityX,2,max) )
  middleValuex <- ( maxValuex + minValuex) / 2
  P1Valuex <- minValuex + ( maxValuex - minValuex ) / 4
  P3Valuex <- middleValuex + ( maxValuex - minValuex ) / 4

  maxValuey <- max( apply(densityY,2,max))
  middleValuey <- maxValuey /2
  minValuey <- min( apply(densityY,2,min))

  haut = seq(minValuey,middleValuey,length.out=(L+1) )

  pal = rainbow(L)
  colors = pal[sample(x=1:L,L)]

  plot(density(phase[,1], n=GridLength), main = title, ylab="Density", xlab = "Date", ylim=c(0,maxValuey+middleValuey), xlim=c(minValuex, maxValuex), bty='n',lty =1, lwd=2, axes=F, col = colors[1])
  lines(density(phase[,2], n=GridLength), lty =1, lwd=2, col = colors[1])

  # abscissa axis
  axis(1, at=c(minValuex, P1Valuex, middleValuex, P3Valuex, maxValuex) , labels =c(floor( minValuex), floor( P1Valuex), floor( middleValuex), floor( P3Valuex), floor( maxValuex)))
  # ordinate axis
  axis(2, at=c(0, middleValuey, maxValuey), labels =c(0, round(middleValuey, 5), round(maxValuey, 5)) )


  ## Following phases
  for(i in 2:L) {
    lines(density(phase[,2*i-1], n = GridLength), col=colors[i], lwd=2, lty=1)
    lines(density(phase[,2*i], n = GridLength), col=colors[i], lwd=2, lty=1)
  }

  ## Phase Time Range
  MPTR = MultiPhaseTimeRange(data=data, position_beginning=position_beginning, position_end=position_end, level=level)
  for (i in 1:L ) { segments(MPTR[i,2], maxValuey+haut[i+1],MPTR[i,3], maxValuey+haut[i+1],lwd=6,col=colors[i]) }
  text(minValuex, maxValuey+haut[2], "Time range", srt =90)
  }
}




####################################
###   Tempo plot   NEW 2016/09   ###

# The tempo plot introduced by T. S. Dye 
# A statistical graphic designed for the archaeological study of rhythms of the long term that embodies a theory of archaeological evidence for the occurrence of events
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position numeric vector containing the position of the column corresponding to the MCMC chains of interest
#' @param level probability corresponding to the level of confidence
#' @param count, if TRUE the counting process is a number, otherwise it is a probability
#' @param Gauss if TRUE, the Gaussian approximation of the CI is used
#' @param title title of the graph
#' @return a plot 
#' @export
TempoPlot <- function(data, position, level = 0.95 , count = TRUE, Gauss = FALSE, title = "Tempo plot") {

  # Construction of a new dataset containing the columns corresponding to the phases of interest
  L = length(position)
  
  groupOfDates = matrix(ncol = L, nrow=nrow(data))
  for (i in 1:L) {
    groupOfDates[,i] = data[,position[i]]
  }
  
   min = min(apply(groupOfDates,2, min))
   max = max(apply(groupOfDates,2, max))
   t = seq( min, max, length.out = 50*ncol(groupOfDates))
  
  f= function(x){
    g=ecdf(x)   
    y=g(t) 
    if (count)  y = y * ncol(groupOfDates)
    y
  } 
  F = t( apply( groupOfDates,1,f ) )
  moy = apply(F,2,mean)
  ec = apply(F,2,sd)
  qu = cbind( apply(F,2, quantile, probs  = (1-level)/2 , type = 8) ,  apply(F,2, quantile, probs  = 1-((1-level)/2) , type = 8)  )
  
  quG=cbind(moy+qnorm(1-(1-level)/2)*ec,moy-qnorm(1-(1-level)/2)*ec)
  
  if (Gauss==TRUE){

    matplot(t,cbind(moy,qu,quG) , lty = 1 ,xlab = "time " , ylab = "counting process " , type="l" , lwd = c(4,1,1,1,1), col=2  )
    polygon( c(t,rev(t)) , c(quG[,1] , rev(quG[,2] )),col = "lightpink", density=4)
    polygon( c(t,rev(t)) , c(qu[,1] , rev(qu[,2] )),col = "lightseagreen", density=4 ,angle= -45)
    lines(t,moy, col=2, lwd = 4)
    legend("bottomright" , legend=c("Bayes estimate", "Credible interval CI", "Gaussian Approx. of CI"), lty=1, col=c(2,"lightseagreen","lightpink"))
    title(title)
    
  }else{

    matplot(t,cbind(moy,qu) , lty = 1 ,xlab = "time " , ylab = "counting process " , type="l" , lwd = c(4,1,1), col=c(2)  )
    polygon( c(t,rev(t)) , c(qu[,1] , rev(qu[,2] )),col = "lightseagreen",density=4 ,angle= -45)
    legend("bottomright" , legend=c("Bayes estimate", "Credible interval CI"), lty=1, col=c(2,"lightseagreen"))
    title(title)
  }
  
}

##############################################
###   Tempo Activity plot   NEW 2016/09   ###

# A statistical graphic designed for the archaeological study of rhythms of the long term that embodies a theory of archaeological evidence for the occurrence of events
#' @param data dataframe containing the output of the MCMC algorithm 
#' @param position numeric vector containing the position of the column corresponding to the MCMC chains of interest
#' @param level probability corresponding to the level of confidence
#' @param count, if TRUE the counting process is a number, otherwise it is a probability
#' @param title title of the graph
#' @return a plot 
#' @export
TempoActivityPlot <- function(data, position, level = 0.95, count = TRUE, title = "Activity plot") {
  
  # Construction of a new dataset containing the columns corresponding to the phases of interest
  L = length(position)
  
  groupOfDates = matrix(ncol = L, nrow=nrow(data))
  for (i in 1:L) {
    groupOfDates[,i] = data[,position[i]]
  }
  
  min = min(apply(groupOfDates,2, min))
  max = max(apply(groupOfDates,2, max))
  t = seq( min, max, length.out = 50*ncol(groupOfDates))
  
  f= function(x){
    g=ecdf(x)   
    y=g(t) 
    if (count)  y = y * ncol(groupOfDates)
    y
  } 
  F = t( apply( groupOfDates,1,f ) )
  moy = apply(F,2,mean)
  ec = apply(F,2,sd)

  plot(x<-t[-1], y<-diff(moy)/diff(t), type = "l", xlab="Time", ylab="Incidence", main =title)  

}


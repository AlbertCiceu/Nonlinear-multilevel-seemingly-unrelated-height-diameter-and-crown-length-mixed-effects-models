function_Z<-function(X_HT, X_CL, dfplot,nSSP) {
  Za <- rbind(X_HT, X_CL)
  if (nSSP == 1) {
    Zb <- rbind(X_HT, X_CL)
  } else{
    nplot1 <- table(dfplot$SSP)[1] #no. obs ssp 1
    nplot2 <- table(dfplot$SSP)[2] #no. obs ssp 2
    total_nPlot<-nplot1+nplot2
    
    block1 <- X_HT[1:nplot1,]
    block1[] <- 0
    block2 <-
      X_HT[(nplot1 + 1):(total_nPlot),]
    block2[] <- 0
    
    HTZp1 <-
      X_HT[1:nplot1,] # Plot 1 partial derivates with the respect to a1, a2
    CLZp1 <-
      X_CL[1:nplot1,] # Plot 1 partial derivates with the respect to b1, b2
    
    Zb1 <- rbind(HTZp1, block2, CLZp1, block2)
    
    HTZp2 <-
      X_HT[(nplot1 + 1):(total_nPlot),] # Plot 2 partial derivates with the respect to a1, a2
    CLZp2 <-
      X_CL[(nplot1 + 1):(total_nPlot),] # Plot 2 partial derivates with the respect to alog, bLog
    Zb2 <- rbind(block1, HTZp2, block1, CLZp2)
    
    Zb <- cbind(Zb1, Zb2)
  }
  Z <- cbind(Za, Zb)
  return(Z)
}

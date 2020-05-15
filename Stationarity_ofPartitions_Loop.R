# Stationarity Tests on combinations of periods through Chi-square test of independence


setwd("C:/Users/hilbert/OneDrive - University of California, Davis/Analytics Software/R_wd/Stationarity")
openWikiPLAY<- read.csv('C:/Users/hilbert/OneDrive - University of California, Davis/Analytics Software/R_wd/Stationarity/openWikiPLAY.csv')
mydata<-c(openWikiPLAY$List.of.WWE.personnel_A7)    #Import dataset
A <- 37  # Alphabet size, adjust to dabase
s <- 5  # the average sample size for each transition count: on average, each transition should have the possibility to have s counts per possible transition

  g <- round(length(mydata)/(A*A*s), digits = 0) # set number of equal sized groups with s samples for division of total time series
  data<- c(tail(mydata, floor(length(mydata)/g)*g))    #eliminate as many early entires as to assure equal group size
  dataG <- as.data.frame(split(data, ceiling(seq_along(data)/(length(data)/g))))  # g equal length groups: dataG$X1, dataG$X2, ...

#run these functions to calculate test metrics before testing data
  
  X2value <- function(Fst, Snd) {
    x<-matrix(nrow=A,ncol=A,0)                     #create matrix for first half for the case of alphabet size A
    for (t in 1:(length(Fst)-1)) x[Fst[t],Fst[t+1]]<-x[Fst[t],Fst[t+1]]+1
    xv<-as.vector(t(x))
    
    z<-matrix(nrow=A,ncol=A,0)                     #create matrix for second half  for the of alphabet size A
    for (t in 1:(length(Snd)-1)) z[Snd[t],Snd[t+1]]<-z[Snd[t],Snd[t+1]]+1
    zv<-as.vector(t(z))
    
    SM <- xv+zv       # calculate side marginals
    Mx <- sum(xv)     # calculate bottom marginal for xv
    Mz <- sum(zv)     # calculate bottom marginal for zv
    TOT <- sum(SM)    # calculate total
    
    Ex <- SM*Mx/TOT   # calculate Expected counts for x assuming independence
    Ez <- SM*Mz/TOT   # calculate Expected counts for z assuming independence
    E <- c(Ex,Ez)     # concatenate both expectated counts under the independence assumption
    O <- c(xv,zv)     # concatenate both observed counts
    
    x2 <- ((O-E)^2)/E  # chi-square raw values
    x2[is.nan(x2)] <- 0  # replace NaN values with 0
    X2<-sum(x2)       # sum up raw values in order to get chi-squared value
    df<-length(xv)-1   # careful: for df it assumes that only two halves are compared: df = (length of vector - 1)*(number of vectors - 1) = in this case (49-1)*(2-1) = 48*1, so the "*1" is omitted
    
    X2  # chi-square value found =
    
  }
  
  
  X2sig99 <- function(Fst, Snd) {
    x<-matrix(nrow=A,ncol=A,0)                     #create matrix for first half for the case of alphabet size A
    for (t in 1:(length(Fst)-1)) x[Fst[t],Fst[t+1]]<-x[Fst[t],Fst[t+1]]+1
    xv<-as.vector(t(x))
    
    z<-matrix(nrow=A,ncol=A,0)                     #create matrix for second half  for the case of alphabet size A
    for (t in 1:(length(Snd)-1)) z[Snd[t],Snd[t+1]]<-z[Snd[t],Snd[t+1]]+1
    zv<-as.vector(t(z))
    
    SM <- xv+zv       # calculate side marginals
    Mx <- sum(xv)     # calculate bottom marginal for xv
    Mz <- sum(zv)     # calculate bottom marginal for zv
    TOT <- sum(SM)    # calculate total
    
    Ex <- SM*Mx/TOT   # calculate Expected counts for x assuming independence
    Ez <- SM*Mz/TOT   # calculate Expected counts for z assuming independence
    E <- c(Ex,Ez)     # concatenate both expectated counts under the independence assumption
    O <- c(xv,zv)     # concatenate both observed counts
    
    x2 <- ((O-E)^2)/E  # chi-square raw values
    x2[is.nan(x2)] <- 0  # replace NaN values with 0
    X2<-sum(x2)       # sum up raw values in order to get chi-squared value
    df<-length(xv)-1   # careful: for df it assumes that only two halves are compared: df = (length of vector - 1)*(number of vectors - 1) = in this case (49-1)*(2-1) = 48*1, so the "*1" is omitted
    
    qchisq(.99, df=df)  # low chi-squared value at .05 level =
    
  }
  
  
  X2df <- function(Fst, Snd) {
    x<-matrix(nrow=A,ncol=A,0)                     #create matrix for first half for the case of alphabet size A
    for (t in 1:(length(Fst)-1)) x[Fst[t],Fst[t+1]]<-x[Fst[t],Fst[t+1]]+1
    xv<-as.vector(t(x))
    
    z<-matrix(nrow=A,ncol=A,0)                     #create matrix for second half  for the case of alphabet size A
    for (t in 1:(length(Snd)-1)) z[Snd[t],Snd[t+1]]<-z[Snd[t],Snd[t+1]]+1
    zv<-as.vector(t(z))
    
    SM <- xv+zv       # calculate side marginals
    Mx <- sum(xv)     # calculate bottom marginal for xv
    Mz <- sum(zv)     # calculate bottom marginal for zv
    TOT <- sum(SM)    # calculate total
    
    Ex <- SM*Mx/TOT   # calculate Expected counts for x assuming independence
    Ez <- SM*Mz/TOT   # calculate Expected counts for z assuming independence
    E <- c(Ex,Ez)     # concatenate both expectated counts under the independence assumption
    O <- c(xv,zv)     # concatenate both observed counts
    
    x2 <- ((O-E)^2)/E  # chi-square raw values
    x2[is.nan(x2)] <- 0  # replace NaN values with 0
    X2<-sum(x2)       # sum up raw values in order to get chi-squared value
    df<-length(xv)-1   # careful: for df it assumes that only two halves are compared: df = (length of vector - 1)*(number of vectors - 1) = in this case (49-1)*(2-1) = 48*1, so the "*1" is omitted
    
    df  # degrees of freedom
    
  }
  
  
  X2p <- function(Fst, Snd) {
    x<-matrix(nrow=A,ncol=A,0)                     #create matrix for first half for the case of alphabet size A
    for (t in 1:(length(Fst)-1)) x[Fst[t],Fst[t+1]]<-x[Fst[t],Fst[t+1]]+1
    xv<-as.vector(t(x))
    
    z<-matrix(nrow=A,ncol=A,0)                     #create matrix for second half  for the case of alphabet size A
    for (t in 1:(length(Snd)-1)) z[Snd[t],Snd[t+1]]<-z[Snd[t],Snd[t+1]]+1
    zv<-as.vector(t(z))
    
    SM <- xv+zv       # calculate side marginals
    Mx <- sum(xv)     # calculate bottom marginal for xv
    Mz <- sum(zv)     # calculate bottom marginal for zv
    TOT <- sum(SM)    # calculate total
    
    Ex <- SM*Mx/TOT   # calculate Expected counts for x assuming independence
    Ez <- SM*Mz/TOT   # calculate Expected counts for z assuming independence
    E <- c(Ex,Ez)     # concatenate both expectated counts under the independence assumption
    O <- c(xv,zv)     # concatenate both observed counts
    
    x2 <- ((O-E)^2)/E  # chi-square raw values
    x2[is.nan(x2)] <- 0  # replace NaN values with 0
    X2<-sum(x2)       # sum up raw values in order to get chi-squared value
    df<-length(xv)-1   # careful: for df it assumes that only two halves are compared: df = (length of vector - 1)*(number of vectors - 1) = in this case (49-1)*(2-1) = 48*1, so the "*1" is omitted
    
    pchisq(X2, df, lower.tail=FALSE)  # probability that both halves are independent = are "the same" = are "stationary" = p =
    
  }
  

# NOW test the data HERE
  
  X2 <- apply(combn(ncol(dataG), 2), 2, function(x) X2value(dataG[,x[1]] , dataG[,x[2]])) # http://stackoverflow.com/questions/22446825/perform-pairwise-comparison-of-matrix 
  sigX2 <- apply(combn(ncol(dataG), 2), 2, function(x) X2sig99(dataG[,x[1]] , dataG[,x[2]]))
  df <- apply(combn(ncol(dataG), 2), 2, function(x) X2df(dataG[,x[1]] , dataG[,x[2]]))
  p <- apply(combn(ncol(dataG), 2), 2, function(x) X2p(dataG[,x[1]] , dataG[,x[2]]))
 
  c <- apply(combn(ncol(dataG), 2), 2, function(m) paste(names(dataG)[m], collapse='&'))
  
X2table <-  rbind(X2,sigX2,df,p)
colnames(X2table) <- c
X2table
"Number of roups:" 
g
"Average number of groups stationary:" 
mean(p)
write.csv(t(X2table), file = "C:/Users/hilbert/OneDrive - University of California, Davis/Analytics Software/R_wd/Stationarity/X2table.csv")




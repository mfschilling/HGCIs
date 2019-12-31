#########################################################################################
#                                                                                       #
#  The function CI.N() generates the confidence intervals for X = 0, 1, ..., n  #
#  for a specified n,M, and confidence level.                                           #
#                                                                                       #
#  Example: To generate all confidence intervals at n=10, M=25, and level=.90, use  #
#                                                                                       #
#                                                                                       #
#  > CI.N(10, 25, .90)                                                              #
#                                                                                       #
#########################################################################################


CI.N <- function(n,M,level)
{ 
  # TODO: Confirm that n and M can be switched, and then make a good comment here about when and why
  # This is necessary because the code breaks if number of successes is less than sample size
  if (M < n){
    swap <- M
    M <- n
    n <- swap
  }
  
  
  acurve<-function(N){phyper(0,M,N-M,n)}
  ii=max(M,n)
  start<-acurve(ii)
  
  while(start<level){
    start=acurve(ii)
    ii=ii+1
  }
  
  max<-ii-1
  
  acurve<-function(N){phyper(1,M,N-M,n)}
  ii=max(M,n)
  start<-acurve(ii)
  
  while(start<level){
    start=acurve(ii)
    ii=ii+1
  }
  
  cutoff.1<-ii-2
  
  if((cutoff.1+1-M)<1){
    cutoff.1<-max
  }
  
  
  
  
  tot.AC <- (n+1)*(n+2)/2   # Total number of acceptance curves is the triangle number of (n+1)
  
  AC.matrix <- matrix(NA,ncol=tot.AC,nrow=(cutoff.1+1-M)) 
  
  N <- seq(M,cutoff.1,1)
  
  #################################################################################################
  # ID.matrix identifies the column number of AC.matrix with its 
  # corresponding [l,u] values. [l,u] values will be included in
  # increasing order of span = u-l.
  
  ID.matrix <- matrix(NA,ncol=4,nrow=tot.AC) 
  colnames(ID.matrix)<-c("col","low","upp","span")
  
  column <- 0
  for (span in 0:n)
  {
    for (k in 0:(n-span))
    {
      column<-column+1
      l <- k
      u <- k+span
      ID.matrix[column,]<-c(column,l,u,(u-l))
    }
  }
  
  #################################################################################################
  # Generation of all acceptance curves
  
  for (cols in 1:tot.AC)
  {
    lower<-ID.matrix[(ID.matrix[,1]==cols),2]
    upper<-ID.matrix[(ID.matrix[,1]==cols),3]
    if (lower==upper){
      AC.matrix[,cols]<- dhyper(lower,M,N-M,n)
    }
    else if (lower==0){
      AC.matrix[,cols] <- phyper(upper,M,N-M,n)
    }
    else {
      AC.matrix[,cols] <- phyper(upper,M,N-M,n)-phyper((lower-1),M,N-M,n)
    }
  }
  
  #################################################################################################
  # Create initial cpf by choosing coverage probbility from highest
  # acceptance curve with minimal span.
  
  cpf.matrix <- matrix(NA,ncol=4,nrow=cutoff.1+1-M) 
  colnames(cpf.matrix)<-c("N","CovProb","low","upp")
  
  for (spot in 1:(cutoff.1+1-M))
  {
    # Determines the column number corresponding to coverage >= level
    # yet having minimum span. Recall that columns are shown in order
    # of INCREASING span. So the min position points to min span.
    min.col <- min(which(AC.matrix[spot,]>=level))
    
    # Looks up the corresponding minimum span of min.col
    min.span <- ID.matrix[(ID.matrix[,1]==min.col),4]
    
    # Pick out all column numbers with that minimum span
    span.cols<-ID.matrix[(ID.matrix[,4]==min.span),1]
    
    # Pick out maximum coverage among ACs with minimum span
    cov.prob <- max(AC.matrix[spot,span.cols])
    
    # Combine column numbers with corresponding cov. prob.
    col.covp <- cbind(span.cols,AC.matrix[spot,span.cols])
    
    # Pick out the column number corresponding to maximum coverage.
    # If multiple columns yield same maximum coverage, choose the
    # smaller column. 
    # Example: n=20, level=.90
    # AC[6,13] and AC[7,14] yield same max value at p = 0.5,
    # so we choose AC[6,13]
    min.span.col <- min(col.covp[(col.covp[,2]==cov.prob),1])
    
    x.lower <- ID.matrix[(ID.matrix[,1]==min.span.col),2]
    x.upper <- ID.matrix[(ID.matrix[,1]==min.span.col),3]
    
    cpf.matrix[spot,]<-c((spot-1+M),cov.prob,x.lower,x.upper)
    
  }
  
  HG.cpf.matrix <- cpf.matrix
  
  #################################################################################################
  # Gap Fix
  # If the previous step yields any violations in monotonicity in [l,u],
  # this will cause a gap. This will check if any violations occur.
  
  seg.cpf.matrix <- matrix(NA,ncol=6,nrow=cutoff.1+1-M) 
  colnames(seg.cpf.matrix)<-c("N","CovProb","low","upp","seg","break")
  
  seg.cpf.matrix[,1:4] <- cpf.matrix[(cutoff.1+1-M):1,1:4]
  
  # Apply segment numbers to p having the same [l,u] acceptance curve
  for (c in 1:nrow(seg.cpf.matrix))
  {
    if (c == 1)
    { 
      seg <- 1
      old.l <- seg.cpf.matrix[c,3]
      old.u <- seg.cpf.matrix[c,4]
      seg.cpf.matrix[c,5] <- seg
    }
    else {
      new.l <- seg.cpf.matrix[c,3]
      new.u <- seg.cpf.matrix[c,4]
      
      if (old.l==new.l & old.u==new.u) seg.cpf.matrix[c,5] <- seg
      else {
        seg <- seg + 1
        old.l <- new.l
        old.u <- new.u
        seg.cpf.matrix[c,5] <- seg
      }
    }
  }
  
  
  # Identify any breaks in [l,u] monotonicity rule
  b.flag <- 0
  for (c in 1:nrow(seg.cpf.matrix))
  {
    if (c == 1)
    { 
      status <- 0
      old.l <- seg.cpf.matrix[c,3]
      old.u <- seg.cpf.matrix[c,4]
      seg.cpf.matrix[c,6] <- status
    }
    else {
      new.l <- seg.cpf.matrix[c,3]
      new.u <- seg.cpf.matrix[c,4]
      
      if (old.l==new.l & old.u==new.u) seg.cpf.matrix[c,6] <- status
      else {
        if (old.l<=new.l & old.u<=new.u) {
          status <- 0
          seg.cpf.matrix[c,6] <- status
        }
        else {
          status <- 1
          seg.cpf.matrix[c,6] <- status
          b.flag <- 1
        }
        old.l <- new.l
        old.u <- new.u
      }
    }
  }
  
  #### If no breaks at all, then jump straight to CI construction
  
  #---------------------------------------------------------------#
  # Start condition on b.flag
  if (b.flag==1){
    
    # Identify the monotonicity break components
    broken<-seg.cpf.matrix[seg.cpf.matrix[,6]==1,]
    broken<-matrix(broken,ncol=6)
    
    # Identify the unique segment numbers of stretches on monotonicity break
    bad.seg <- unique(broken[,5])
    
    #.................................#
    
    # Perform cpf alteration
    for (b in bad.seg)
    {
      target.sub <- matrix(seg.cpf.matrix[(seg.cpf.matrix[,5]==(b+1)),],ncol=6,byrow=F)
      
      target.lower <- min(target.sub[,3])
      target.upper <- min(target.sub[,4])
      
      current.sub <- matrix(seg.cpf.matrix[(seg.cpf.matrix[,5]==(b)),],ncol=6,byrow=F)
      
      for (g in 1:nrow(current.sub))
      {
        point <- current.sub[g,1]
        
        # Update with the new [l,u] endpoints we'll be using
        current.sub[g,3] <- target.lower
        current.sub[g,4] <- target.upper
        
        if (target.lower==target.upper){
          current.sub[g,2]<- dhyper(target.lower,M,point-M,n)
        }
        else if (target.lower==0){
          current.sub[g,2] <- phyper(target.upper,M,point-M,n)
        }
        else {
          current.sub[g,2] <- phyper(target.upper,M,point-M,n)-phyper((target.lower-1),M,point-M,n)
        }
      }
      
      # Replace all rows of half.cpf corresponding to range of p's from 
      # current.sub with the entire current.sub. This will update with
      # the new coverage probabilities and corresponding [l,u] labels
      
      seg.cpf.matrix[(seg.cpf.matrix[,1]>=min(current.sub[,1]) & seg.cpf.matrix[,1]<=max(current.sub[,1])),]<-current.sub
      
    }
    
    #.................................#
    
    HG.cpf.matrix <- matrix(NA,ncol=4,nrow=cutoff.1+1-M) 
    colnames(HG.cpf.matrix)<-c("N","CovProb","low","upp")
    
    HG.cpf.matrix <- seg.cpf.matrix[(cutoff.1+1-M):1,1:4]
    
    
    #---------------------------------------------------------------#
  }
  
  # End condition on b.flag
  
  #################################################################################################
  # CI Generation
  
  ci.matrix <-  matrix(NA,ncol=3,nrow=n+1)
  rownames(ci.matrix) <- c(rep("",nrow(ci.matrix)))
  colnames(ci.matrix) <- c("x","lower","upper")
  
  for (x in 0:n)
  {
    ref.sub <- HG.cpf.matrix[(HG.cpf.matrix[,3]<=x & x<=HG.cpf.matrix[,4]),]
    ref.sub <- matrix(ref.sub,ncol=4)
    low.lim <- ref.sub[1,1]
    upp.lim <- ref.sub[nrow(ref.sub),1]
    ci.matrix[x+1,]<-c(x,low.lim,upp.lim)
  }
  
  ci.matrix[1,3] = Inf
  ci.matrix[2,3] = max-1
  
  heading <- matrix(NA,ncol=1,nrow=1)    
  
  heading[1,1] <- paste("Confidence Intervals for n = ",n,", M = ",M," and level = ",level,sep="")
  
  rownames(heading) <- c("")
  colnames(heading) <- c("")
  
  print(heading,quote=FALSE)
  
  # print(round(ci.matrix,digits=dp))
  print(ci.matrix)
  
  #################################################################################################
}


#########################################################################################
#                                                                                       #
#  Example: To generate all confidence intervals at n=10, M=25, and level=.90, use  #
#                                                                                       #
#                                                                                       #
#  > CI.N(10, 25, .90)                                                              #
#                                                                                       #
#########################################################################################
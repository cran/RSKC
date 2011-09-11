RSKC.trimkmeans <- function(data,k,trim=0.1, scaling=FALSE,
                            runs=100, points=NULL,
                       countmode=runs+1,
                       maxit=2*nrow(as.matrix(data))){
  # the code of trimkmeans is modified so that
  # for alpha = 0, no case is trimmed. i.e. k-means clustering is performed                     	
  data <- as.matrix(data)
  n <- nrow(data); nc<-ncol(data)
  nin <- ceiling((1-trim)*n)
  if (scaling) data <- scale(data)
  crit <- Inf
  oldclass <- iclass <- optclass <- rep(0,n)
  disttom <- rep(0,n)
#  optmeans <- data[sample(n,k),,drop=FALSE]
  for (i in 1:runs){
    if ((i/countmode)==round(i/countmode)) cat("Iteration ",i,"\n")
    if (is.null(points))
      means <- data[sample(n,k),,drop=FALSE]
    else
      means <- points
    wend <- FALSE
    itcounter <- 0
    while(!wend){
      itcounter <- itcounter+1
      reF<-.Fortran(
                    "disttom_iclass",
                    as.double(data),as.integer(n),as.integer(nc),
                    as.double(means),as.integer(k),
                    iclass=as.integer(iclass),
                    disttom=as.double(disttom),
                    PACKAGE="RSKC"
                    )            
      disttom<-reF$disttom
      iclass<-reF$iclass
    
      # define outliers **it could be that some clusters do not have obs**
      iclass0<-iclass;iclass0[order(disttom)[(nin+1):n]] <- 0
      if (trim!=0) iclass<-iclass0
#    newcrit <- sum(disttom[iclass>0])
#    cat("Iteration ",i," criterion value ",newcrit,"\n")

      # stopping criteria if the class labels are the same then stop the iteration..
      if (itcounter>=maxit | identical(oldclass,iclass)) wend <- TRUE
      else{
        for (l in 1:k){
        # if a cluster is empty then cluster center is the first outlier obs ??
          if (sum(iclass==l)==0) means[l,] <- data[iclass0==0,,drop=FALSE][1,]
          else{
            if (sum(iclass==l)>1){
              if (dim(means)[2]==1)
                means[l,] <- mean(data[iclass==l,])
              else
                means[l,] <- colMeans(data[iclass==l,])
            }
            else means[l,] <- data[iclass==l,]
          }
        }
        oldclass <- iclass
      }
    }
    newcrit <- sum(disttom[iclass>0])
    if (newcrit<=crit){
      optclass <- iclass
      crit <- newcrit
      optmeans <- means
    }
  }
  optclass[optclass==0] <- k+1
  out <- list(classification=optclass,means=optmeans,
              criterion=crit/nin,disttom=disttom,ropt=sort(disttom)[nin],
              k=k,trim=trim,runs=runs,scaling=scaling)
  #class(out) <- "tkm"
  out
}

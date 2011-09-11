RSKC.trimkmeans.missing <- function(data,k,w,
sumW=sum(w),trim=0.1,scaling=FALSE, runs=100, points=NULL,
                       countmode=runs+1, printcrit=FALSE,
                       maxit=2*nrow(as.matrix(data))){
   # w must be; length(w)==ncol(data)
  data <- as.matrix(data);
  n <- nrow(data);p <-ncol(data)
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
      reF <- wrapper.fort(data,nrow=n,ncol=p,mu=means,k=k,W=w,sumW=sumW)
      disttom <- reF$disttom
      iclass <- reF$iclass

      iclass0<-iclass;iclass0[order(disttom)[(nin+1):n]] <- 0
      if (trim!=0) iclass<-iclass0
    
      # stopping criteria if the class labels are the same
      # then stop the iteration..
      if (itcounter>=maxit | identical(oldclass,iclass)) wend <- TRUE
      else{
        for (l in 1:k){
        # if a cluster is empty then cluster center is the farthest case from 
        # its cluster center
        # if this is true, then cluster center can contains missing value..
          if (sum(iclass==l)==0) means[l,] <- data[iclass0==0,,drop=FALSE][1,]
          else{
            if (sum(iclass==l)>1){
              if (dim(means)[2]==1) # if only one feature in data
                means[l,] <- mean(data[iclass==l,],na.rm=TRUE)
              else
                means[l,] <- colMeans(data[iclass==l,],na.rm=TRUE)
            }
            # if the l^th cluster cotains only one element,
            # cluster center can contains missing values..
            else means[l,] <- data[iclass==l,] 
          }
        }
       oldclass <- iclass
      }
    }
    newcrit <- sum(disttom[iclass>0])
    if (printcrit) cat("Iteration ",i," criterion value ",
                       newcrit/nin,"\n")
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



wrapper.fort <- function(data,nrow,ncol,mu,k,W,sumW)
            {
              miss_data <- matrix(1,nrow,ncol)
              miss_data[is.na(data)] <- 0
              data[miss_data==0] <- 0
              
              miss_mu <- matrix(1,k,ncol)
              miss_mu[is.na(mu)] <- 0
              mu[miss_mu==0] <- 0
              # outputs
              iclass <- rep(0,nrow)
              disttom <- rep(0,nrow)
              
              result<-.Fortran("disttom_iclass_missing",
                            as.double(data),
                            as.integer(miss_data),
                            as.integer(nrow),as.integer(ncol),
                            as.double(mu),
                            as.integer(miss_mu),
                            as.integer(k),
                            iclass=as.integer(iclass),
                            disttom=as.double(disttom),
                            as.double(W),
                            as.double(sumW),
                            PACKAGE="RSKC")
              if (sum(is.na(result$disttom))!=0)
                stop("L1 is too small (sparse) for a dataset with missing values!")
              return(list(disttom=result$disttom,iclass=result$iclass))
            }


WDISC <-
function(D,mu,ncl,n,w,sumW){
 # D, mu and w must be in reduced dimention	
 # ncol(D)==length(w)
  dist<-matrix(NA,nrow=n,ncol=ncl);  
  for ( i in 1:ncl){
          s <- sc <- scale(D, center=mu[i,], scale=FALSE)^2
          s[ is.na(sc) ] <- 0;vec<-rowSums(s)

          adjust<-sumW/((!is.na(sc))%*%w) # n by 1 adjusted sclaers 
          dist[,i]<-vec*adjust 
           }      
  return(dist)
  }
  

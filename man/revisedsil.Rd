\name{revisedsil}
\alias{revisedsil}
\title{
The revised silhouette
}

\description{
This function returns a revised silhouette plot, 
cluster centers in weighted squared Euclidean distances and  
a matrix containing the weighted squared Euclidean distances between cases and each cluster center.
Missing values are adjusted.
}



\usage{
revisedsil(d,reRSKC=NULL,CASEofINT=NULL,col1="black",
	CASEofINT2 = NULL, col2="red", print.plot=TRUE, 
	W=NULL,C=NULL,out=NULL)
}

\arguments{


\item{d}{
A numerical data matrix, \code{N} by \code{p}, where \code{N} is the number of cases and \code{p} is the number of features.
}


\item{reRSKC}{
 A list output from RSKC function.
}



\item{CASEofINT}{
Necessary if print.plot=TRUE.
A vector of the case indices that appear in the revised silhouette plot. 
The revised silhouette widths of these indices are colored in \code{col1} if \code{CASEofINT != NULL}.
The average silhouette of each cluster printed in the plot is computed EXCLUDING these cases.
}

\item{col1}{
See \code{CASEofINT}.
}

\item{CASEofINT2}{
A vector of the case indices that appear in the revised silhouette plot.
The indices are colored in \code{col2}.
}

\item{col2}{
See \code{CASEofINT2}
}

\item{print.plot}{
If \code{TRUE}, the revised silhouette is plotted.
}

\item{W}{
Necessary if \code{reRSKC = NULL}. A positive real vector of weights of length \code{p}.
}

\item{C}{
Necessary if \code{reRSKC = NULL}.  An integer vector of class labels of length \code{N}.
}

\item{out}{
Necessary if \code{reRSKC = NULL}. 
Vector of the case indices that should be excluded in the calculation of cluster centers.
In \code{RSKC}, cluster centers are calculated without the cases that have the furthest 100*\code{alpha} \% Weighted squared Euclidean distances to their closest cluster centers.
If one wants to obtain the cluster centers from \code{RSKC} output, set \code{out} = \code{<RSKCoutput>$oW}. 
}


}


\value{

\item{trans.mu}{
Cluster centers in reduced weighted dimension.
See example for more detail.
}

\item{WdisC}{
\code{N} by \code{ncl} matrix, where \code{ncl} is the prespecified number of clusters. It contains the weighted distance between each case and all cluster centers.
See example for more detail.
}

\item{sil.order}{
Silhouette values of each case in the order of the case index.
}
\item{sil.i}{
Silhouette values of cases ranked by decreasing order within clusters.
The corresponding case index are in \code{obs.i}
}

}


%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\author{
Yumi Kondo <y.kondo@stat.ubc.ca>
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\references{
Yumi Kondo (2011), Robustificaiton of the sparse K-means clustering algorithm, MSc. Thesis, University of British Columbia
\url{http://hdl.handle.net/2429/37093}
}
\examples{

# little simulation function 
sim <-
function(mu,f){
   D<-matrix(rnorm(60*f),60,f)
   D[1:20,1:50]<-D[1:20,1:50]+mu
   D[21:40,1:50]<-D[21:40,1:50]-mu  
   return(D)
   }


### output trans.mu ###

p<-200;ncl<-3
# simulate a 60 by p data matrix with 3 classes 
d<-sim(2,p)
# run RSKC
re<-RSKC(d,ncl,L1=2,alpha=0.05)
# cluster centers in weighted squared Euclidean distances by function sil
sil.mu<-revisedsil(d,W=re$weights,C=re$labels,out=re$oW,print.plot=FALSE)$trans.mu
# calculation 
trans.d<-sweep(d[,re$weights!=0],2,sqrt(re$weights[re$weights!=0]),FUN="*") 
class<-re$labels;class[re$oW]<-ncl+1
MEANs<-matrix(NA,ncl,ncol(trans.d))
for ( i in 1 : 3) MEANs[i,]<-colMeans(trans.d[class==i,,drop=FALSE])
sil.mu==MEANs
# coincides 

### output WdisC ###

p<-200;ncl<-3;N<-60
# generate 60 by p data matrix with 3 classes 
d<-sim(2,p)
# run RSKC
re<-RSKC(d,ncl,L1=2,alpha=0.05)
si<-revisedsil(d,W=re$weights,C=re$labels,out=re$oW,print.plot=FALSE)
si.mu<-si$trans.mu
si.wdisc<-si$WdisC
trans.d<-sweep(d[,re$weights!=0],2,sqrt(re$weights[re$weights!=0]),FUN="*") 
WdisC<-matrix(NA,N,ncl)
for ( i in 1 : ncl) WdisC[,i]<-rowSums(scale(trans.d,center=si.mu[i,],scale=FALSE)^2)
# WdisC and si.wdisc coincides

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


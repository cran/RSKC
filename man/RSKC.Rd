\name{RSKC}
\alias{RSKC}
\alias{RSKC.trimkmeans}
\alias{RSKC.trimkmeans.missing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robust Sparse K-means}
\description{
The robust sparse \emph{K}-means clustering method by Kondo (2011). In this algorithm, sparse \emph{K}-means (Witten and Tibshirani (2010)) is robustified by iteratively trimming the prespecified proportion of cases in the weighted squared Euclidean distances and the squared Euclidean distances.
}
\usage{
RSKC(d, ncl, alpha, L1 = 12, nstart = 200, 
silent=TRUE, scaling = FALSE, correlation = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

\item{d}{
A numeric matrix of data, \code{N} by \code{p}, where \code{N} is the number of cases and \code{p} is the number of features. 
Cases are partitioned into \code{ncl} clusters. Missing values are accepted. 
}

\item{ncl}{
The prespecified number of clusters.
}

\item{alpha}{
0 <= \code{alpha} <= 1, the proportion of the cases to be trimmed in robust sparse \emph{K}-means. 

If \code{alpha} > 0 and \code{L1} >= 1 then \code{RSKC} performs robust sparse \emph{K}-means. 

If \code{alpha} > 0 and \code{L1} = \code{NULL} then \code{RSKC} performs trimmed \emph{K}-means. 

If \code{alpha} = 0 and \code{L1} >=1 then \code{RSKC} performs sparse \emph{K}-means (with the algorithm of Lloyd (1982)).
  
If \code{alpha} = 0 and \code{L1} = \code{NULL} then \code{RSKC} performs \emph{K}-means (with the algorithm of Lloyd). 

For more details on trimmed \emph{K}-means, see Gordaliza (1991a), Gordaliza (1991b).
}


 \item{L1}{
A single L1 bound on weights (the feature weights). If \code{L1} is small, then few features will
have non-zero weights. If \code{L1} is large then all features
will have non-zero weights.
If \code{L1} = \code{NULL} then \code{RSKC} performs nonsparse clustering (see \code{alpha}).
}
  \item{nstart}{
The number of random initial sets of cluster centers in every step (a) which performs \emph{K}-means or trimmed \emph{K}-means.
}
\item{silent}{
If \code{TRUE}, then the processing step is not printed.
}

\item{scaling}{
If \code{TRUE}, \code{RSKC} subtracts the each entry of data matrix by the corresponding column mean and divide it by the corresponding column SD: see \code{\link{scale}}
}

\item{correlation}{
If \code{TRUE}, \code{RSKC} centers and
scales the rows of \code{data} before the clustering is performed. i.e., \code{trans.d = t(scale(t(d)))}
The squared Euclidean distance between cases in the transformed dataset \code{trans.d} is proportional to the dissimilality measure based on the correlation between the cases in the dataset \code{d}
}

}

\value{
\item{N}{The number of cases.}
\item{p}{The number of features.}
\item{ncl}{See \code{ncl} above.}
\item{L1}{See \code{L1} above.}
\item{nstart}{See \code{nstart} above.}
\item{alpha}{See \code{alpha} above.}
\item{scaling}{See \code{scaling} above.}
\item{correlation}{See \code{correlation} above.}
\item{missing}{It is \code{TRUE} if at least one point is missing in the data matrix, \code{d}.}
          

\item{labels}{
An integer vector of length \code{N}, set of cluster labels for each case.
Note that trimmed cases also receive the cluster labels.
}

\item{weights}{
A positive real vector of length \code{p}, containing weights on each feature.}

\item{WBSS}{
A real vector containing the weighted between sum of squares at each Step (b). The weighted between sum of squares is the objective function to maximize, excluding the prespecified proportions of cases.
The length of this vector is the number of times that the algorithm iterates the process steps (a),(a-2) and (b) before the stopping criterion is satisfied.
This is returned only if \code{L1} is numeric and > 1.
}

\item{WWSS}{
A real number, the within cluster sum of squares at a local minimum.
This is the objective function to minimize in nonsparse methods.
For robust clustering methods, this quantity is calculated without the prespecified proportions of cases.
This is returned only if \code{L1}=\code{NULL}, 
}

\item{oE}{
Indices of the cases trimmed in squared Euclidean distances.
}

\item{oW}{
Indices of the cases trimmed in weighted squared Euclidean distances. 
If \code{L1} =\code{NULL}, then \code{oW} are the cases trimmed in the Euclidean distance, because all the features have the same weights, i.e., \code{1}'s.}
}

\details{
Robust sparse \emph{K}-means is a clustering method that extends the sparse \emph{K}-means clustering of Witten and Tibshirani to make it resistant to oultiers by trimming a fixed proportion of observations in each iteration. 
%
These outliers are flagged both in terms of their weighted and unweighted distances to eliminate the effects of outliers in the selection of feature weights and the selection of a partition.
%
In Step (a) of sparse \emph{K}-means, given fixed weights, the algorithm aims to maximize the objective function over a partition i.e. it performs \emph{K}-means on a weighted dataset.
Robust sparse \emph{K}-means robustifies Step (a) of sparse \emph{K}-means by performing trimmed \emph{K}-means on a weighted dataset: it trims cases in weighted squared Euclidean distances.  
%
Before Step (b), where, given a partition, the algorithm aims to maximize objective function over weights, the robust sparse 
\emph{K}-means has an intermediate robustifying step, Step (a-2).
At this step, it trims cases in squared Euclidean distances.
%
Given a partition and trimmed cases from Step (a) and Step (a-2), the objective function is maximized over weights at Step(b).
The objective function is calculated without the trimmed cases in Step (a) and Step(a-2).
%
The robust sparse \emph{K}-means algorithm repeat Step (a), Step (a-2) and Step (b) until a stopping criterion is satisfied.
%
For the calculation of cluster centers in the weighted distances, see \code{revisedsil}.
}

\references{
Y. Kondo, M. Salibian-Barrera, R.H. Zamar. RSKC: An R Package for a Robust and Sparse K-Means Clustering Algorithm.,Journal of Statistical Software, 72(5), 1-26, 2016.
	
A. Gordaliza. Best approximations to random variables based on trimming procedures. Journal of Approximation Theory, 64, 1991a.
 
A. Gordaliza. On the breakdown point of multivariate location estimators based on trimming procedures. Statistics & Probability Letters, 11, 1991b. 

Y. Kondo (2011), Robustificaiton of the sparse K-means clustering algorithm, MSc. Thesis, University of British Columbia
\url{http://hdl.handle.net/2429/37093}

D. M. Witten and R. Tibshirani. A framework for feature selection in 
clustering. Journal of the American Statistical Association, 105(490)
713-726, 2010. 

S.P. Least Squares quantization in PCM. IEEE Transactions on information theory, 28(2): 129-136, 1982.
}

\author{
Yumi Kondo <y.kondo@stat.ubc.ca>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
# little simulation function 
sim <-
function(mu,f){
   D<-matrix(rnorm(60*f),60,f)
   D[1:20,1:50]<-D[1:20,1:50]+mu
   D[21:40,1:50]<-D[21:40,1:50]-mu  
   return(D)
   }

set.seed(1);d0<-sim(1,500)# generate a dataset
true<-rep(1:3,each=20) # vector of true cluster labels
d<-d0
ncl<-3
for ( i in 1 : 10){
   d[sample(1:60,1),sample(1:500,1)]<-rnorm(1,mean=0,sd=15)
}

# The generated dataset looks like this...
pairs(
      d[,c(1,2,3,200)],col=true, 
      labels=c("clustering feature 1",
      "clustering feature 2","clustering feature 3",
      "noise feature1"),
      main="The sampling distribution of 60 cases colored by true cluster labels", 
      lower.panel=NULL) 


# Compare the performance of four algorithms
###3-means
r0<-kmeans(d,ncl,nstart=100)
CER(r0$cluster,true)

###Sparse 3-means
#This example requires sparcl package
#library(sparcl)
#r1<-KMeansSparseCluster(d,ncl,wbounds=6)
# Partition result
#CER(r1$Cs,true)
# The number of nonzero weights
#sum(!r1$ws<1e-3)

###Trimmed 3-means

r2<-RSKC(d,ncl,alpha=10/60,L1=NULL,nstart=200)
CER(r2$labels,true)

###Robust Sparse 3-means
r3<-RSKC(d,ncl,alpha=10/60,L1=6,nstart=200)
# Partition result
CER(r3$labels,true)
r3

### RSKC works with datasets containing missing values...
# add missing values to the dataset
set.seed(1)
for ( i in 1 : 100)
{   
d[sample(1:60,1),sample(1,500,1)]<-NA
}
r4 <- RSKC(d,ncl,alpha=10/60,L1=6,nstart=200)

}
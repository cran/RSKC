\name{Clest}
\alias{Clest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An implementation of Clest with robust sparse K-means. 
CER is used as a similarity measure.
}
\description{
The function \code{Clest} performs Clest ( Dudoit and Fridlyand (2002)) with CER as the measure of the agreement between two partitions (in each training set).
The following clustering algorithm can be used: \emph{K}-means, trimmed \emph{K}-means, sparse \emph{K}-means and robust sparse \emph{K}-means. 
}

\usage{
Clest(d, maxK, alpha, B = 10, B0 = 20, nstart = 100, 

      L1 = 12, beta = 0.1, pca = TRUE, silent=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
A numerical data matrix (\code{N} by \code{p}) where \code{N} is the number of cases and \code{p} is the number of features. The cases are clustered.
}
  \item{maxK}{
The maximum number of clusters that you suspect.
}
  \item{alpha}{
See \code{\link{RSKC}}.
}
  \item{B}{
The number of times that an observed dataset \code{d} is randomly partitioned into a learning set and a training set.
Note that each generated reference dataset is partitioned into a learning and a testing set only once to ease the computational cost.
}

  \item{B0}{
The number of times that the reference dataset is generated.
}
  \item{nstart}{
The number of random initial sets of cluster centers at Step(a) of robust sparse \emph{K}-means clustering.
}
  \item{L1}{
See \code{\link{RSKC}}.
}

  \item{beta}{
0 <= \code{beta} <= 1: significance level. 
Clest chooses the number of clusters that returns the strongest significant evidence against the hypothesis H0 : K = 1. 
}
  \item{pca}{
Logical, if \code{TRUE}, then reference datasets are generated from a PCA reference distribution. 
If \code{FALSE}, then the reference data  set is generated from a simple reference distribution.
}

   \item{silent}{
Logical, if \code{TRUE}, then the number of iteration on progress is not printed.
}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{

\item{K}{
The solution of Clest; the estimated number of clusters.
}

\item{result.table}{
A real matrix (\code{maxK-1} by 4). 
Each row represents \code{K=2},...,\code{maxK} and columns represent the test statistics (=observed CER-reference CER), observed CER, reference CER and \emph{P}-value. 
}

\item{referenceCERs}{
A matrix (\code{B0} by \code{maxK-1}), containing CERs of testing datasets from generated datasets for each \code{K=2,...,maxK}.
}

\item{observedCERs}{
A matrix (\code{B} by \code{maxK-1}), containing CERs of \code{B} testing sets for each \code{K=2,...,maxK}.
} 

\item{call}{The matched call.}

}

\references{
Yumi Kondo (2011), Robustificaiton of the sparse K-means clustering algorithm, MSc. Thesis, University of British Columbia
\url{http://hdl.handle.net/2429/37093}

S. Dudoit and J. Fridlyand. A prediction-based resampling method for estimating the number of clusters in a dataset. Genome Biology, 3(7), 2002. 
}


\author{
Yumi Kondo  <y.kondo@stat.ubc.ca>
}

%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
\dontrun{
# little simulation function 
sim <-
function(mu,f){
   D<-matrix(rnorm(60*f),60,f)
   D[1:20,1:50]<-D[1:20,1:50]+mu
   D[21:40,1:50]<-D[21:40,1:50]-mu  
   return(D)
   }
 
 set.seed(1)
 d<-sim(1.5,100); # non contaminated dataset with noise variables
 
# Clest with robust sparse K-means
rsk<-Clest(d,5,alpha=1/20,B=3,B0=10, beta = 0.05, nstart=100,pca=TRUE,L1=3,silent=TRUE);
# Clest with K-means
k<-Clest(d,5,alpha=0,B=3,B0=10, beta = 0.05, nstart=100,pca=TRUE,L1=NULL,silent=TRUE);
}
}

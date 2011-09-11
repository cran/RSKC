\name{CER}
\alias{CER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classification Error Rate (CER)
}
\description{
Compute the classification error rate of two partitions. 
}
\usage{
CER(ind, true.ind,nob=length(ind))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ind}{
Vector, containing the cluster labels of each case of a partition 1. 
}
  \item{true.ind}{
Vector, containing the cluster labels of each case of a partition 2. 
}

\item{nob}{
The number of cases (the length of the vector ind and true ind)
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Return a CER value.
CER = 0 means perfect agreement between two partitions and CER = 1 means complete disagreement of two partitions.
Note: 0 <= \code{CER} <= 1
}
\references{
H. Chipman and R. Tibshirani. Hybrid hierarchical clustering with 
applications to microarray data. Biostatistics, 7(2):286-301, 2005.
}
\author{
Yumi Kondo <y.kondo@stat.ubc.ca>
}

\note{
This function uses \code{comb}, which generates all combinations of the elements in the vector \code{ind}.
For this reason, the function \code{CER} is not suitable for vector in a large dimension.  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
vec1<-c(1,1,1,2,3,3,3,2,2)
vec2<-c(3,3,3,1,1,2,2,1,1)
CER(vec1,vec2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

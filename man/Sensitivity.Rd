\name{Sensitivity}
\alias{Sensitivity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the sensitivities (probability of true positive) of each cluster
}
\description{
The sensitivity of cluster \emph{k} is calculated as follows:
%
First, the proportions of observations whose true cluster label is \emph{k}
are computed for each classified clusters.
%
Then the largest proportion is selected as sensitivity.
%
Since this calculation can return 1 for sensitivities of all clusters if all observations
belong to one cluster, we also report the observed cluster labels
returned by the algorithms.
}

\usage{
Sensitivity(L,trueL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{
The cluster label vector from any clustering algorithms.
}
 
\item{trueL}{
The true cluster labels.
}
}

\author{
Yumi Kondo <y.kondo@stat.ubc.ca>
}



%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
vec1<-c(1,1,1,2,3,3,3,2,2)
vec2<-c(3,3,3,1,1,2,2,1,1)
Sensitivity(vec1,vec2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

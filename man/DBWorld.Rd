\name{DBWorld}
\alias{DBWorld}
\alias{rawDBWorld}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
E-mails from DBWorld mailing list
}

\description{
The dataset contains n= 64 bodies of e-mails in binary bag-of-words
representation which Filannino manually collected
from DBWorld mailing list.
%
DBWorld mailing list announces conferences, jobs, books, software and
grants.
%
Filannino applied supervised learning algorithm to classify e-mails
between ``announces of conferences'' and ``everything else''.
%
Out of 64 e-mails, 29 are about conference announcements and 35 are not.


Every e-mail is represented as a vector containing p binary values,
where p is the size of the vocabulary extracted from the entire
corpus with some constraints:
the common words such as ``the'', ``is'' or ``which'', so-called stop words, 
and words that have less than 3 characters or more than 30 chracters
are removed from the dataset. 
%
The entry of the vector is 1 if the corresponding word belongs to the
e-mail and 0 otherwise.
%
%
The number of unique words in the dataset is p=4702.
%
The dataset is originally from the UCI Machine Learning
Repository DBWorldData.


\code{rawDBWorld} is a list of 64 objects containing the original E-mails.
}

\usage{
data(DBWorld)
data(rawDBWorld)
}

%\arguments{
%}

\details{
See Bache K, Lichman M (2013). for details of the data descriptions. 
The original dataset is freely available from USIMachine Learning Repository website 
\url{http://archive.ics.uci.edu/ml/datasets/DBWorld+e-mails}
}

\references{
Bache K, Lichman M (2013). UCI Machine Learning Repository." \url{http://archive.ics.uci.edu/ml/datasets}

Filannino, M., (2011). 'DBWorld e-mail classification using a very small corpus', Project of Machine Learning course, University of Manchester. 
}
\author{
Yumi Kondo <y.kondo@stat.ubc.ca>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

\examples{
\dontrun{
data(DBWorld)
data(rawDBWorld)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

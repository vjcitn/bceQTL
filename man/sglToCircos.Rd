\name{sglToCircos}
\alias{sglToCircos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
use RCircos to display SNP-gene associations
}
\description{
use RCircos to display SNP-gene associations
}
\usage{
sglToCircos(sgl, chr.exclude = NULL, snpgr = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sgl}{
list with vectors of gene symbols, names are SNP ids
}
  \item{chr.exclude}{
for RCircos
}
  \item{snpgr}{
if supplied, start() will be used to provide SNP locations, otherwise
BSgenome snpsById will be used
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%%  ~Describe the value returned
%%%  If it is a LIST, use
%%%  \item{comp1 }{Description of 'comp1'}
%%%  \item{comp2 }{Description of 'comp2'}
%%% ...
%}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(normTab1)
data(bcGR)
sglToCircos(normTab1, snpgr=bcGR)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")

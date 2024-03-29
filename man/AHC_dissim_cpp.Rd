% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{AHC_dissim_cpp}
\alias{AHC_dissim_cpp}
\title{Dissimilarity method using C++}
\usage{
AHC_dissim_cpp(X, method = "average", dissim_matrix = FALSE)
}
\arguments{
\item{X}{Data matrix by default (if dissim_matrix = True, then this is the dissimilarity matrix)}

\item{method}{type of linkage for the Lance-Williams distance ("average" by default)}

\item{dissim_matrix}{false by default, put true if you give the dissimilarity matrix directly}
}
\value{
P clustering matrix (each row corresponds to each iteration)
}
\description{
Dissimilarity method using C++
}
\examples{
X = as.matrix(iris[-5])
method = "single"
P = AHC_dissim_cpp(X, method)

}

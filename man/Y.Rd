\name{Y}
\alias{Y}
\docType{data}
\title{An example single cell count expression matrix (Yan)}
\source{
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36552}
}
\description{
Y is a count expression matrix which belongs to "matrix" class. The data includes 124 cells about human preimplantation embryos and embryonic stem cells. It contains 19304 genes after removing genes with extreme high dropout rate.
}
\references{
Yan, Liying, et al. "Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells." Nature structural & molecular biology 20.9 (2013): 1131.
}
\format{An object of "matrix" class contains the count expressions}
\usage{data("Yan")}
\examples{
data("Yan")
Y[1:10, 1:4]
}
\keyword{datasets}

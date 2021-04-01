\name{trueclass}
\alias{trueclass}
\docType{data}
\title{An example single cell dataset for the cell label information (Yan)}
\description{
  The true cell type labels for Yan dataset. It includes 8 different cell types.
}
\source{
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36552}
}
\references{
Yan, Liying, et al. "Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells." Nature structural & molecular biology 20.9 (2013): 1131.
}
\format{A character vector contains the cell type label}
\usage{data("Yan")}
\examples{
data("Yan")
table(trueclass)
}
\keyword{datasets}

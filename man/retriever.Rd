% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/retriever.R
\name{retriever}
\alias{retriever}
\title{Generate robust disease-specific drug-response profiles}
\usage{
retriever(cellLines, corThreshold = 0.6)
}
\arguments{
\item{cellLines}{A vector of cell lines from where the disease-specific drug response profiles will be generated.}

\item{corThreshold}{A numeric value between 0 and 1 that represent the minimum Spearman correlation coefficient value required by a profile to remain in the analysis.}
}
\value{
Robust disease-specific drug-response profiles that represent the unique transcriptional changes induced by a compound in a given disease.
}
\description{
Generate disease-specific drug-response profiles that are independent of time, concentration, and cell-line. Based on the cell lines used as surrogates, the returned profiles represent the unique transcriptional changes induced by a compound in a given disease.
}
\details{
In the first step, we take the response profile of a given cell line to the same compound under the same concentration at different time points and averaged them. Then, the descriptive power of the generated profile to represent the drug response at a given concentration in the cell line, independently of time, is evaluated using Spearman's correlation coefficient. If the correlation with the generated profile is larger than 0.6, then the averaged profile is returned. The original drug response profiles that do not reach the threshold are removed, and the averaged profile is recomputed using the ones above the threshold, this procedure ensures the removal of aberrant or insufficient cellular responses. Only averaged signatures of at least two profiles are used in the second step.

In the second step, we take the stable time-independent signature profiles of the compounds at a particular concentration in the same cell line. To remove the concentration dependency, we applied the same procedure described in the first step over the averaged profiles. The profiles returned by the second step are the stable ones independent of the time and concentration under the same cell line.

In the last step, to generate disease-specific drug response profiles we again applied the procedure described in the step 1 to the stable response profiles to the same compounds in the three cell lines used as surrogate of the triple-negative breast cancer in the LINCS-L1000 project. The profiles returned by the third step are robust disease-specific transcriptional signatures representing the changes that are unique for a compound in triple-negative breast cancer.
}
\examples{
\donttest{

# Generate a robust profiles across different breast cancer cell lines.
BRCA <- retriever(cellLines = c('MDAMB231', 'MCF7', 'SKBR3', 'HS578T', 'BT20'))
}
}
\author{
Daniel Osorio <daniecos@uio.no>
}

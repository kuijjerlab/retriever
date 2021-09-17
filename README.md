# <img src="https://raw.githubusercontent.com/kuijjerlab/retriever/main/inst/retriever.png" width="30" title="retriever logo"> retriever
An R package to generate robust disease-specific response signatures from the LINCS-L1000 data that are independent of time, concentration, and cell-line. Based on the cell lines used as surrogates, the returned profiles represent the unique transcriptional changes induced by a compound in a given disease.

## Method
The procedure to generate disease-specific drug-response profiles was divided into three steps. The first is to remove the time dependency, the second is to remove the concentration dependency, and the third is to remove the cell line dependency and generate disease-specific robust drug response profiles. 

1 - (Panel A) We took the response profile of a given cell line to the same compound under the same concentration at different time points and averaged them. Then, the descriptive power of the generated profile to represent the response to the drug at a given concentration in the cell line independently of the time was tested using Spearman’s correlation coefficient. If the correlation with the generated profile is larger than 0.6, then the averaged profile is returned. The original drug response profiles that do not reach the threshold are removed, and the averaged profile is recomputed using the ones above the threshold, this procedure ensures the removal of aberrant or insufficient cellular responses. Only averaged signatures of at least two profiles are used in the second step. 

2 - (Panel B) We took the stable time-independent signature profiles of the compounds at a particular concentration in the same cell line. To remove the concentration dependency, we applied the same procedure described in the first step over the averaged profiles. The profiles returned by the second step are the stable ones independent of the time and concentration under the same cell line. 

3 - (Panel C) To generate disease specific drug-response profiles we again applied the procedure described in the step 1 to the stable response profiles to the same compounds in this case, in three cell lines used as surrogate of the triple-negative breast cancer in the LINCS-L1000 project. Thus, the profiles returned by the third step are robust disease-specific transcriptional signatures representing the changes that are unique for a compound in triple-negative breast cancer. 

![method](https://raw.githubusercontent.com/kuijjerlab/retriever/main/inst/method.png)

## Usage
**retriever** is under active development, you can install it, using the following command:
```{r}
library(remotes)
install_github('kuijjerlab/retriever')
library(retriever)
```

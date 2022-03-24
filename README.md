
# msiScores

Statistical & multivariate scores implementation as described Erich et al. (reference below). The scores can be used for systematic evaluation of workflows and comparison of sample preperation methods to check for reproducibility and reliability of the generated data from different preparation protocols.

### Installation

The [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package could be used to install _msiScores_.

```r
install.packages("devtools")
devtools::install_github("CeMOS-Mannheim/msiScores")
```
Note that _msiScores_ was created with `renv`to help manage R package dependencies and computational reproducibility. for more info, check [the renv guide page](https://rstudio.github.io/renv/articles/renv.html). 

### Example

``` r
library(msiScores)


data("exampleScores")
test <- getScores(mat, sample(c("1", "2"), nrow(mat), replace = TRUE), compute.r2 = FALSE, compute.fc = FALSE)

```
### Citing _msiScores_

Erich, Katrin, et al. "Scores for standardization of on-tissue digestion of formalin-fixed paraffin-embedded tissue in
MALDI-MS imaging." Biochimica et Biophysica Acta (BBA)-Proteins and Proteomics 1865.7 (2017): 907-915.

### Contact

You are welcome to submit suggestions and bug-reports at: <https://github.com/CeMOS-Mannheim/msiScores/issues>.

### License

See [license document](LICENSE.md).

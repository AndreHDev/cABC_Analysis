cABCanalysis
=============

R and python implementation of the computed ABC Analysis for rational selection of the most informative variables. Features caclulation using monotonic splines, edge-case handling (see details [here](https://github.com/AndreHDev/cABC_Analysis/blob/main/R/cABC_special_cases.R)), base R plots as well as GGPLots.

<img width="1072" height="602" alt="cABCExample" src="https://github.com/user-attachments/assets/9612fb9b-6ba9-487d-8235-f6dfb4c6791f" />

Install R package
-----------------
From the cloned repository root, install the package from the current directory:

```r
install.packages("path/to/repo", repos = NULL, type = "source")
library(cABCanalysis)
```

Reference
-----------------
If you use this tool in your research, please cite the original publication:

Ultsch A, Lötsch J (2015) "Computed ABC Analysis for Rational Selection of Most Informative Variables in Multivariate Data". PLoS ONE 10(6): e0129767.

[doi:10.1371/journal.pone.0129767](https://doi.org/10.1371/journal.pone.0129767)

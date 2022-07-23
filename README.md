# Global Adaptive Generative Adjustment (GAGA)

GAGA is a parameter estimation method, and it has efficient model selection ability. Its hyperparameters are easy to set. In general, alpha can be set to 1, 2 or 3, but when the collinearity of the load matrix is serious, the hyperparameters can be selected larger, such as 5. In the comparative experiment of dealing with linear models, the GAGA algorithm is competitive with adaptive lasso, SCAD, MCP. In the comparative experiment of dealing with generalized linear models, the GAGA algorithm is competitive with glmnet.

At present, GAGA can deal with "gaussian", "binomial", "poisson", "multinomial" and "cox" problems in linear models.

For the first edition of the paper, see: <https://arxiv.org/abs/1911.00658v1> The first edition of the paper only proves the theoretical properties of Gaga linear model when the load matrix is orthogonal. At present, Xiaofei Wang has completed the proof of the theoretical properties when the load matrix is non-orthogonal. We have submitted the manuscript and hope to be admitted.

## Setup
### If OS is Windows
Install from Github:

1. Install Rtools from https://cran.r-project.org/bin/windows/Rtools/

2. Install devtools package **install.packages("devtools")**

3. Install GAGA package  **devtools::install_github("WangBin0431/GAGA")**

4. Use **? GAGA** to get help doc 

Install from package file:

1. Download GAGA package from: https://drive.google.com/file/d/17MbbfWuBXhK1EUuElf8WyJ_wQoaGjmRS/view?usp=sharing

2. Install it.

3. Use **? GAGA** to get help doc

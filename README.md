# Global Adaptive Generative Adjustment (GAGA)

GAGA is a parameter estimation method, and it has efficient model selection ability. Its hyperparameters are easy to set. In general, alpha can be set to 1, 2 or 3, but when the collinearity of the load matrix is serious, the hyperparameters can be selected larger, such as 5. In the comparative experiment of dealing with linear models, the GAGA algorithm is competitive with adaptive lasso, SCAD, MCP. In the comparative experiment of dealing with generalized linear models, the GAGA algorithm is competitive with glmnet.

At present, GAGA can deal with "gaussian", "binomial", "poisson", "multinomial" and "cox" problems in linear models.

For the first edition of the paper, see: <https://arxiv.org/abs/1911.00658v1> The first edition of the paper only proves the theoretical properties of Gaga linear model when the load matrix is orthogonal. At present, Xiaofei Wang has completed the proof of the theoretical properties when the load matrix is non-orthogonal. We have submitted the manuscript and hope to be admitted.

## Setup

Install from Github:

1.  Install Rtools from <https://cran.r-project.org/bin/windows/Rtools/>

2.  Install devtools package: `install.packages("devtools")`

3.  Install GAGA package: `devtools::install_github("WangBin0431/GAGA")`

4.  Use `?GAGA` to get help doc

Install from package file:

1.  Download GAGA package from: <https://drive.google.com/drive/folders/1IGaYlAIedbjBLExdf8VdooUlhZ9jLTAH?usp=sharing>

2.  Install it.

3.  Use `?GAGA` to get help doc

## Tip

Improve R performance by installing optimized BLAS libraries.

**If you like this program, please give me a star. Thank you.**

# Participants

Bin Wang, eatingbeen@hotmail.com, <https://scholar.google.com/citations?user=7HrcFqEAAAAJ&hl=en>

Xiaofei Wang, <wangxf341@nenu.edu.cn>, <https://scholar.google.com/citations?user=WBnzo_QAAAAJ&hl=en>

Jianhua Guo, jhguo@nenu.edu.cn

Huaian Diao, diao@jlu.edu.cn, <https://scholar.google.com/citations?user=p9M8-YgAAAAJ&hl=en>


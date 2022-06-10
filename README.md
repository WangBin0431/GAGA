# GAGA
Global Adaptive Generative Adjustment  (GAGA) algorithm
GAGA is a parameter estimation method, and it has efficient model selection ability. Its hyperparameters are easy to set. In general, alpha can be set to 2 or 3, but when the collinearity of the load matrix is serious, the hyperparameters can be selected larger, such as 5. In our comparative experiments, the GAGA algorithm is competitive with adaptive lasso, SCAD, MCP.
It can be applied to linear and non-linear model problems. We only upload the R version of the linear model for now. Please wait for other versions.

For the first edition of the paper, see: https://arxiv.org/abs/1911.00658v1
The first edition of the paper only proves the theoretical properties of the case where the load matrix is orthogonal. The proof for the case where the load matrix is not orthogonal has been completed and has been submitted, hoping to be accepted.

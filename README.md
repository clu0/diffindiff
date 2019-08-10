# diffindiffdev: Robust Non-parametric Difference-in-differences estimation

This package implements methods for non-parametric difference-in-differences estimation, as proposed by Hirshberg and Wager (2017). We consider a
setup where we observe data `(X, Y, S, T)` generated according
to the following model:
```
X ~ P(X)
S ~ P(S|X), T ~ P(T|X), where S, T takes values {0,1}
Y = b(X) + T*f(X) + S*h(X) + T*S*tau(X) + noise
```
and we want to estimate the average effect `TAU = E[tau(X)]`. The general framework is flexible: auxilary parameters such as `b(X)`, `f(X)` and `h(x)` can be estimated with any preferred machine learning method. Our methods are implemented with the lasso (glmnet), and `cv.glmnet` is used to perform cross-fitting.

### Installation

To install this package in R, run the following commands:
```R
library(devtools) 
install_github("clu0/diffindiff")
```
### Example usage

```R
library(diffindiff)
library(glmnet)
library(rlearner)
library(EQL)

n = 100; p=12
X = matrix(rnorm(n * p), n, p)
tau = 1
f = X[,3]
h = X[,4]
S.treat.prob = 0.4
T.treat.prob = 0.4
Si = rbinom(n, 1, S.treat.prob)
Ti = rbinom(n, 1, T.treat.prob)
b = pmax(X[,1] + X[,2], 0)
Y = b + Ti*f + Si*h + T*S*tau + rnorm(n)

# generate a basis
make_matrix = function(x) stats::model.matrix(~.-1, x)
X = data.frame(X) %>% make_matrix
order = 3
Basis = generate.basis(X,order)

# If we assume that the treatment effect tau is constant, which is true
# for the setup above, then we can run the transformed regression
TR_fit = DiD(X = Basis,
  Y = Y,
  T = T,
  S = S,
  constant_eff = "constant")
  
# If we do not assume tau is constant, we can run DiD-AMLE
gamma.minimax = minimax(Basis, T, S)
DiD_AMLE_fit = DiD(X = Basis,
  Y = Y,
  T = T,
  S = S,
  constant_eff = "non_constant",
  gamma = gamma.minimax)

```

#### References
Chen Lu, Xinkun Nie, Stefan Wager.
<b>Robust Nonparametric Difference-in-Differences
Estimation.</b>
2019.
[<a href="https://arxiv.org/pdf/1905.11622.pdf">arxiv</a>]

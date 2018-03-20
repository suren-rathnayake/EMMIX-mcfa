## MCFA: Mixtures of Common Factor Analyzers

(NOTE: This pacakge is now incorporated into EMMIXmfa.)

Fits a mixture of factor analyzers with a common component matrix for the
factor loadings before the transformation of the latent factors to be
white noise. It is designed specifically for the task of displaying the
observed data points in a lower (_q_-dimensional) space, 
where _q_ is the number of factors adopted in the
factor-analytic representation of the observed vector.

It also provides a greater reduction in the number of parameters in the model.
Component distributions can either be from the family of multivariate normals
or from the family of multivariate _t_-distributions.
Maximum likelihood estimators of model parameters are obtained using the Expectation-Maximization algorithm.

Fitting a MCFA model with there components using two factors for the Iris data available in
R can be done using,  
```
fit <- mcfa(Y = iris[, -5], g = 3, q = 2)
```

The groupings can be visualized in the _q_-dimensional factor space.
```
plot_factors(fit)
```
MCFA fits multivariate normals to the data, fitting _t_-distributions can be achieved
using `mctfa` function. Further, there are functions to generate data from a `emmix-mcfa`
models (`rmix`), estimate factor scores (`factor_scores`), estimate adjusted Rand Index (`ari`),
find the number of misallocations (`err`), among others.

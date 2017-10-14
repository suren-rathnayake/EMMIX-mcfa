## MCFA: Mixtures of Common Factor Analyzers

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

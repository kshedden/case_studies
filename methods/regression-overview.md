# Overview of regression methods

## Introduction

Regression analysis is arguably the most widely-used tool in applied
statistics, and has also inspired many important developments in statistical
theory.  Here we discuss some concepts that can be used to understand
several of the most widely-used approaches to regression.  Then we review
some specific regression methods along with their key properties.

Before proceeding, note that regression itself is somewhat difficult to
define in a way that differentiates it from the rest of statistics.  In most
cases, regression focuses on a conditional distribution, e.g. the conditional
distribution of a variable $y$ given another variable $x$.  The symbols $y$
and $x$ may denote either scalars or vectors.  Any analysis focusing on such
a conditional distribution can be seen as a form of regression analysis.

The variables $y$ and $x$ have various essentially equivalent names.
The variable $y$ is known as the "outcome", the "response variable", and
the "dependent variable", while "x" may be referred to as a "covariate", a
"predictor", a "regressor", or an "explanatory variable".

## Major concepts

* __Prediction__: A major use case for regression is to predict (or
"forecast") unobserved values of the dependent variable $y$ when the
explanatory variables $x$ are observed.  However not all applications
of regression analysis focus on prediction.  In basic research we may use
regression analysis to understand whether and how two quantities are related,
even when there is no practical need to make predictions.  It is also often
of interest to understand relationships between variables that are too
weak to be useful for prediction.  This can give insight into mechanisms
relating the phenomena measured by $x$ and $y$.

* __Mean regression__: this term refers to any regression analysis where the
population target is the conditional mean function $E[y | x]$; also known as
"mean structure modeling".

* __Variance regression__: this is a class of approaches that model the
conditional variance ${\rm var}[y|x]$ alone, or that model the conditional
variance along with the conditional mean.  A common technique is to model
$\log {\rm var}[y|x]$ since a variance is always positive.

* __Quantile regression__: this refers to any method that models a conditional
quantile in terms of covariates.  A common form of quantile reqression is
"median regression", which is sometimes favored over more conventional
approaches targeting the conditional mean because it is more "robust".
Quantile regression using probability points other than 0.5 can be used to
understand the tails of the conditional distribution of $y$ given $x$.

* __Single index models__: a single index model is any regression model that
is expressed in terms of one "linear predictor" involving the covariates, i.e.
$b_1x_1 + \cdots + b_px_p$, where the $x_j$ are observed covariates (data)
and the $b_j$ are unknown coefficients (parameters).  This means that all
the information in $x$ that is related to $y$ is contained in the given
linear combination of the components of $x$.

* __Linear model__: Depending on the context, this term can mean any of the
following: (i) some property of the conditional distribution, most commonly the
conditional mean, is modeled as a linear function of the covariates, (ii) some
property of the conditional distribution is linear in the parameters, or (iii)
the fitted values and/or parameter estimates are linear in the observed values
of the response variable.  In common use, "linear model" refers to a linear
model for the mean structure, e.g. $E[y|x] = \alpha + \beta x$, that is fit to
the data using linear least squares.  This model always has properties (ii)
and (iii) (in reference to the conditional mean).  It may not have property
(i), since $E[y|x] = \alpha + \beta x + \gamma x^2$ can also be considered
to be a linear model with properties (i) and (iii) but does not have property
(ii).  It is incorrect to assert that linear models must have property (ii),
and as a result cannot represent regression relationships between $y$ and $x$
that are nonlinear.  An example of a model that does not have properties (i),
(ii), or (iii) is the mean structure model $E[y|x] = \exp(\alpha + \beta x)$.

* __Regression for independent observations__: Most basic regression methods
are directly applicable for samples of independent observations. These basic
methods may or may not give meaningful results when the observations are
statistically dependent.  More advanced regression methods have been devised
for use when the observations are statistically dependent.

* __Heteroscedasticity__: If the conditional variance ${\rm Var}[y|x]$ is
constant (i.e. does not depend on $x$), then the conditional distribution
of $y$ given $x$ is homoscedastic, otherwise it is heteroscedastic
(alternative terminology is "constant variance" and "nonconstant variance").
Heteroscedasticity can be accommodated by some regression procedures
(e.g. Poisson regression works best if the mean and variance are equal,
which is a strong form of heteroscedasticity).  Some regression procedures,
like ordinary least squares, work best when the population is homoscedastic,
but can still give meaningful results (with some loss of power) if the
population is heteroscedastic.

* __Mean/variance relationship__: If the variance of a distribution is a
function of the mean, there is a mean/variance relationship.  For example,
in the Poisson distribution, the variance is equal to the mean, and
in the negative binomial distribution the variance has the form $\mu +
\alpha\mu^2$, where $\alpha \ge 0$ is a "shape parameter".  However in the
Gaussian distribution, the variance is unrelated to the mean (i.e. it is a
constant function of the mean) so there is no mean/variance relationship.

* __Overdispersion/underdispersion__: these terms refer to misspecification of
a working regression model relative to the population.  If the true conditional
variance is greater than the conditional variance of the model being fit to
the data, then there is overdispersion.  If the true conditional variance
is less than the conditional variance of the model being fit to the data
there is underdispersion.

* __Repeated measures__: This is one reason that data may be non-independent.
Repeated measures (or "clustering") refers to any setting in which the data
fall into groups, and the observations in any one group are more similar to
each other than they are to observations in other groups.  This can happen,
for example, due to unobserved covariates that are stable (constant) within
a group.

* __Marginal regression__: This is a form of regression analysis where the
estimation target is the marginal regression function $E[y|x]$, even though
the data may be clustered or otherwise dependent.  Some methods for marginal
regression also give insight into the marginal variance function ${\rm
Var}[y|x]$ and marginal covariances ${\rm Cov}[y_1, y_2|x_1, x_2]$.  If the
data are dependent, a "complete" model should target the joint distribution
$P(y_1, \ldots, y_n | x_1, \ldots, x_n)$.  A marginal model averages over
most of the data, allowing us to focus on simpler "marginal" relationships.

* __Multilevel regression__: This is an alternative term for "random effects
modeling".  It emphasizes the fact that in many data sets, there are complex
inter-relationships between the observations that are not explained by
the covariates.  These inter-relationships allow us to speak in terms of
unobserved "random effects" that are incorporated into the linear predictors
of one or more observations.  The random effects give rise to dependence,
and also, in nonlinear models, they give rise to different ways of defining
a "regression effect".  Multilevel models can also be viewed as a way to
model variances and covariances, but modeling them through random effects,
rather than directly.

* __Conditional/marginal effect__ (in multilevel regression): In a
multilevel model, a "marginal effect" is usually defined as the change in
$E[y|x]$ corresponding to a one unit change of a specific covariate $x_k$.
A "conditional effect" is usually defined as the change in $E[y|x,u]$ for a one
unit change in $x_k$, where $u$ is an unobserved random effect.  For linear
models, conditional and marginal effects are the same, but in nonlinear
models the two types of effects differ.  Methods for nonlinear regression
target either the marginal effects, or the conditional effects, but usually
not both.  In most cases the conditional effect will be numerically larger
than the marginal effect.  Note that the word "effect", while widely used,
conveys causality that may not be warranted.

* __Conditional/marginal effect__ (in single-level regression): Another
use of the term "marginal effect" arises in single-level regression models.
In this case the marginal effect is the change in $E[y|x_k]$ corresponding
to a one unit change of $x_k$, while the conditional effect is the change
in $E[y| x_1, \ldots, x_p]$ corresponding to a single unit change in $x_k$,
with the other variables $x_j$ for $j\ne k$ held fixed.  When referring to
this type of marginal effect, the marginal and conditional effects differ
even in a linear model.

* __Nonparametric regression__: There isn't always a bright line between
"parametric" and "nonparametric" regression methods, but in general parametric
methods are less "flexible" and may only work well if certain strong conditions
on the population hold.  Nonparametric methods may work well in a broader
class of population settings, but often with lower power and precision.

## Models, fitting procedures, and algorithms

It is important to distinguish between the various regression model
structures (e.g. different model parameterizations), and different ways for
fitting a regression model structure to data.  For example, the "linear mean
structure" model is one prominent structural model for regression, in which
the conditional mean function $E[y|x]$ is expressed as a linear function of
the predictors in $x$ (i.e. in terms of a linear predictor).  There are many
"fitting procedures" (algorithms) that enable one to fit this linear model to
data, including least squares, penalized least squares, and many variations
of robust regression, maximum likelihood regression, and Bayesian regression.
However all of these fitting procedures are fitting the same model structure
to the data.

As an example, least squares is a fitting procedure that can be used to fit a
model to data.  The least squares fitting procedure has statistical properties
(i.e. it is known to be efficient, consistent, etc. in some settings).
A different (e.g. Bayesian or penalized) procedure for fitting the same class
of models will have its own, potentially different properties (e.g. it may
be consistent in some settings where least squares is not and vice-versa).

Algorithms are specific numerical procedures used to fit a model to data.
For example, we may use iteratively reweighted least squares to fit a
generalized linear model to a dataset.  In many cases, the algorithm exactly
minimizes an explicit loss-function, and therefore the algorithm itself does
not impact the statistical properties of the analysis (e.g. we can use the
QR or SVD approaches to solving the linear least squares problem and always
get the same result unless the model is degenerate). In a few settings,
e.g. nearest neighbor regression, regression trees, or deep neural networks,
people may say that "the algorithm is the model".  In these settings,
there is a mean structure model, but the model is extremely flexible and
the properties of the procedure result mostly from the algorithm rather than
from the structural form of the model.

## Some specific regression analysis methods

* __Least squares__: ordinary least squares (OLS) is the most basic type of
curve fitting.  It is most effective when the conditional mean function is
linear in the covariates, and the conditional variance is constant (i.e. we
have "homoscedasticity").  Both of these restrictions can be worked around,
however.  Nonlinearity of the mean function can be accommodated using basis
functions, and heteroscedasticity can be accommodated using inverse variance
weights (in which case were are doing "weighted least squares", or WLS).
Also, heteroscedasticity only impacts statistical efficiency, which may not
be a major concern when fitting simple models to large datasets.

* __Generalized Linear Models (GLM)__: GLM's are an extension of linear
models that introduce _link functions_ and _mean/variance relationships_.
The link function allows the expected value of the response variable to be
expressed as a known transformation of the linear predictor.  The mean/variance
relationship expresses how the conditional variance of the response given
the predictors relates to the conditional mean of the response given the
predictors.  GLM's are often (but not always) a better alternative to using
linear least squares with a transformed predictor (e.g. instead of regressing
$\log y$ on $x$ using a linear model, regress $y$ on $x$ using a GLM with
a log link function).  GLMs can be used as a likelihood-based approach,
but there is also a quasi-likelihood theory for GLMs that justifies their
use in much broader settings, e.g. we can use Poisson regression if the data
are not integer-valued, or do not follow a Poisson distribution, or even if
the mean/variance relationship is violated.

* __Generalized Estimating Equations (GEE)__: GEE is an extension of GLM that
allows for certain types of statistical dependencies between the observations.
A GEE is determined by specifying the GLM that it is derived from, and a
"working model" for the correlation structure.  The fitting and inference in
a GEE is robust in that the working dependence model can be misspecified,
and the estimates and inferences will still be valid (this can be stated
in more precise terms but we will not do that here).  GEE estimates the
"marginal mean structure".  In the linear case, GEE is closely related to
the more basic technique of "generalized least squares" (GLS).

* __Multilevel linear models__: multilevel (or mixed) linear models are an
extension of the basic linear model in which there are (usually) one or more
covariates, and also "random effects" which describe how the observations
are correlated with each other.  These unobserved random effects can be
viewed as missing information that reflects additional structure in the
population not captured through the covariates.  There is essentially
a 1-1 correspondence between mixed linear models and GLS/GEE models,
in that both estimate the same population target (the conditional mean
function), but using different estimators.  The mixed linear model will
in most cases give better estimates of variance parameters than GLS/GEE,
but may be less robust to misspecification of the dependence structure.
It is a very rich framework that can be used to account for a variety of
structures in the population that are difficult to model in other ways,
including clustering, multilevel (nested) clustering, crossed clustering,
and heterogeneous partial associations (e.g. the coefficient for a covariate
differs across many known subpopulations).

* __Multilevel GLM's__: these are one of the most challenging classes of
regression models, especially from a computational perspective.  Structurally,
they are very similar to linear mixed models, and in practice, can be
interpreted in a similar way, except for the important distinction that
in a multilevel GLM, the marginal and conditional mean structures differ
(which is not the case for a multilevel linear model).

### Other forms of regression:

* __Survival regression__ -- this is a large set of techniques used
   for handling censored data

* __Conditional regression__ -- this is a useful but narrowly applicable
"trick" in which by conditioning on certain statistics, a multilevel model is
essentially converted into a single-level model.  The most familiar forms of
this technique are single-level conditional logistic and Poisson regression.
In both cases, we can have clustered data (which would more often be handled
using mixed effects or GEE), but by conditioning on the observed total of
the outcome values within each group, the observations become conditionally
independent, and can be rigorously fit using a single-level likelihood
approach.

* __Local regression__ -- this is a very flexible approach to capturing
nonlinear regression relationships.  It is an example of a regression method
that is not fitting a single-index model; it is generally seen as being
limited by the "curse of dimensionality", so that it cannot be applied with
more than a handful of covariates, unless the sample size is very large.

* __Additive regression__ -- this is a way to restrict the general kernel
regression technique to avoid the curse of dimensionality.  The conditional
mean function $E[y|x]$ is modeled as $g_1(x_1) + \cdots + g_p(x_p)$, where
the $g_j()$ are unknown univariate functions.  The model is additive over
the covariates, which is a strong restriction, but generalizes classical
linear models by allowing each covariate to be transformed in an arbitrary way.

* __Dimension reduction regression__ -- this is a very unique and distinct
class of regression approaches that posit a multi-index structure and an
unknown link function.  Specifically, $E[y|x]$ may be modeled as having the
form $g(b_1^\prime x, \ldots, b_k^\prime x)$, where the $b_j$ are vectors
of regression coefficients, and $g$ is an unknown link function.  The focus
is on estimating the regression "directions" $b_j$, not on the link function.

* __Generalized method of moments (GMM)__ -- this is a technique for
efficiently estimating the parameters of nonlinear models using only the
moments.  It is mainly used when it is important to estimate regression
effects without requiring a model for the full conditional distribution to
be specified.

* __Multivariate regression__ -- these are techniques for regressing a vector
of dependent variables on a vector of independent variables.  Information is
shared across the regressions to allow them to be fit more accurately as a
collection, compared to performing separate regressions to the components
of the vector dependent variable.

* __Machine learning/algorithmic regression__ -- this is a broad and
loosely-defined collection of methods for regression analysis in which complex
representations of regression functions, e.g. trees, ensembles of trees, or
neural networks, are fit to data, usually with substantial regularization.
The distinction between "machine learning" and "statistical regression" is
generally artificial and not useful to make, but some techniques, especially
neural networks, are often viewed as being part of machine learning.

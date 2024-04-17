# Multilevel regression

Regression analysis aims to understand the conditional distribution of
a scalar response $y$ in relation to explanatory variables
$x \in {\mathbb R}^p$.  The conditional mean $E[y | x]$ is an important
feature of this conditional distribution, but in many cases we may
also be interested in the conditional variances ${\rm Var}[y | x]$ and the
conditional covariances ${\rm Cov}[y | x]$.  Multilevel regression is a
framework for regression analysis that is especially useful if we
have covariances between observations.

Conditional covariances often arise due to the way that the data were
collected.  Formally, the data are statistically dependent if the joint
pdf does not factor as the product of marginal pdf's, i.e.

$$
P(y_1, \ldots, y_n | x_1, \ldots, x_n) \ne \prod_i P(y_i | x_i).
$$

Multilevel regression is a means to understand the
conditional mean $E[y|x]$, conditional variance ${\rm Var}(y | x)$,
conditional covariances ${\rm Cov}(y_i, y_j | x_i, x_j)$, and other
forms of dependence among the observations.

A common manner in which correlated data arise in practice is when
the data are collected as *repeated measures*.  For example,
we may be studying a characteristic of individual people such as income,
and we collect this data from each person every year for multiple years.
Two repeated measures of one person's income are likely to be more similar
than incomes for two different people, so there is an
[intraclass correlation](https://en.wikipedia.org/wiki/intraclass_correlation)
between two income observations made on the same person.
Data collected in this way are called *longitudinal data*.

Another typical setting in which correlated data arise is when we have
a *clustered sample*.  For example, suppose that we randomly sample
communities, and then within each community we randomly sample 100
people and obtain their incomes.  Since two people in a community tend
to have similar incomes, any two observations made in the same
community will be correlated.

Generically, we can refer to a collection of repeated measures that
may be statistically dependent as a
*block*.  If we have longitudinal data, then each person is a block.
If we have a cluster sample, then each cluster is a block.

## A single level of blocking

Let $y_{ij}$ denote the $j^{\rm th}$ repeated measure for the $i^{\rm
th}$ block.  One way to account for the correlations in the data is to
introduce a *random effect*, with the most basic random effect being a
*random intercept*.  A random intercept is a random variable
$\theta_i$ that arises in the following model:

$$
y_{ij} = \beta^\prime x_{ij} + \theta_i + \epsilon_{ij}.
$$

In this model, $i=1, \ldots n$ indexes blocks and $j=1, \ldots, n_i$
indexes observations within block $i$ ($n_i$ is the size of block
$i$).

The *mean structure* is parameterized through the linear predictor

$$
\beta^\prime x_{ij} = \beta_0 + \beta_1x_{ij1} + \cdots \beta_p x_{ijp}.
$$

This linear predictor is exactly the same as would be present in a
conventional linear model.

The random intercepts $\theta_i$ are a collection of independent and
identically distributed (IID) random variables assumed to follow a
distribution with mean zero and variance $\tau^2$. Finally, we have
*unexplained variation* specific to each observation that is represented
through the random variables
$\epsilon_{ij}$, which are IID random variables with mean
zero and variance $\sigma^2$ that are independent of the $\theta_i$.

We can now study the marginal moments of the multilevel model.  The
marginal mean is

$$
E[y_{ij} | x_{ij}] = \beta^\prime x_{ij},
$$

since the random effects $\theta_i$ and the unexplained "errors"
$\epsilon_{ij}$ all have mean zero.

The marginal variance is

$$
{\rm var}[y_{ij} | x_{ij}] = {\rm var}(\theta_i) + {\rm var}(\epsilon_{ij}) = \tau^2 + \sigma^2.
$$

Note that

$$
{\rm var}(\theta_i + \epsilon_{ij}) = {\rm var}(\theta_i) + {\rm var}(\epsilon_{ij})
$$

since $\theta_i$ and $\epsilon_{ij}$ are independent.

The marginal covariance between two observations in the same block is

$$
{\rm cov}[y_{ij}, y_{ij^\prime} | x_{ij}, x_{ij^\prime}] = {\rm cov}(\theta_i+\epsilon_{ij}, \theta_i+\epsilon_{ij^\prime}) = \tau^2.
$$

Further, the marginal correlation between two observations in the same
block is $\tau^2/(\tau^2+\sigma^2)$, which is also known as the
*intra-class correlation*, or *ICC*.

Two observations in different blocks are independent so the covariance
and correlation between them is zero.

It is important to note that the random effects $\theta_i$ are neither
data (which must be observed) nor are they parameters which are
unknown values to be estimated based on the data (e.g. using maximum likelihood).
Instead, the random
effects are "latent variables" that are not observed and are integrated out
of the model before estimating the parameters using a procedure such
as maximum likelihood estimation.

The parameters of the model discussed above are $\beta$, $\tau^2$, and
$\sigma^2$.  All these parameters are estimated jointly, usually via
maximum likelihood estimation (MLE) or the closely-related restricted
maximum likelihood estimation (REML) which we will not define here.

### Alternatives to multilevel regression

Random effects models are not the only way to accommodate dependence
in a regression analysis.  Some alternatives are:

* Robust inference -- this involves fitting a regression model ignoring
any dependence induced by the sampling design (e.g. longitudinal or
cluster-wise dependence) and then adjusting the standard errors or
test statistics.  For example, suppose we use OLS when the data are
dependent, yielding an estimate $\check{\beta}$.  It can be shown that
if the mean structure is linear, $\check{\beta}$ and $\hat{\beta}$ are
both valid estimates of $\beta$.  However the standard errors and test
statistics obtained by OLS are not correct when dependence is present.
A procedure known as "robust inference", "Huber-White inference", or
"sandwich inference" can be used to obtain valid standard errors.

* Estimating equations regression -- this involves reweighting the
OLS estimator to obtain efficient estimates with valid standard
errors in the presence of dependence.  The most common implementation
of this approach is "generalized estimating equations" (GEE).

* Fixed effects regression -- this involves "one hot coding" the
group variable to produce a potentially large set of dummy variables
that are included in the regression along with the other covariates.
There is an extensive literature on the relationship between fixed
effects and random effects approaches to modeling clustered data.
An advantage of the fixed effects approach is that it does not
assume a model for the cluster effects.
The main drawback of the fixed effects approach is that when we have
a large number of small groups, there is not enough information to
estimate the fixed effects as parameters, and the estimates of
parameters for variables of interest may be inconsistent.  This
is known as the *Neyman-Scott problem*.

## Model formulas and long form data

Most modern software for fitting multilevel models expects the data to
be in long form.  That is, instead of having the data indexed by two
subscripts $y_{ij}$ with $i$ indexing blocks and $j$ indexing
observations within blocks, in long form we have only a single index
$i$ that runs over all observations in all blocks.  A second variable $b$ is used to
denote the block for each observation.  For example, we may have a
response variable $y$, a block variable $b$, and covariates $x_1$ and
$x_2$.  Observation $y_i$ belongs to block $b_i$ and has covariates
$x_{1i}$, $x_{2i}$, etc.

The "random intercepts" model discussed above can be expressed using
the formula *y ~ x1 + x2 + (1 | b)*, where the special syntax
*(1 | b)* indicates that the model has a random intercept, with distinct
values of $b$ defining groups that each have their own random
intercept.

## Random slopes

Now consider a certain explanatory variable, say $x_1$.  The "fixed"
slope for $x_1$ is $\beta_1$, meaning that the average response
differs by $\beta_1$ units for two observations that differ by one
unit in $x_1$ and are identical in all other covariates.  However it
may be that this variable operates differently in different blocks.
That is, there may be a slope $\gamma_i$ specifically for block $i$,
leading to the following model:

$$
y_{ij} = \beta^\prime x_{ij} + \theta_i + \gamma_i x_{1ij} + \epsilon_{ij}.
$$

Here we view $\gamma_i$ as a latent random variable with mean zero and
variance $\tau_1^2$, and as above $\theta_i$ is random with mean zero
and variance $\tau^2$.

The formula for a model with random slopes can be expressed with the
formula *y ~ x1 + x2 + (1 + x1 | b)*.

Random slopes are a way to capture *heterogeneous treatment effects*,
a situation where the association between a covariate $x$ and the outcome
$y$ differs based on the group $b$.

## Predicting random effects

There are three distinct types of values in a multilevel regression:
the data, the random effects, and the parameters.

The data are always observed, and above would be denoted $y$, $x$, and
$b$.  The data can be further partitioned into data that are modeled
as random (here $y$), and data that are conditioned upon so can be
viewed as either fixed (deterministic) or random (here $x$ and $b$).

Parameters are unknown deterministic (nonrandom) quantities to be
estimated.  In a multilevel analysis, the parameters can be
partitioned into the *mean structure parameters* $\beta$ and the
*variance structure parameters* $\sigma^2$ and $\tau^2$.

The random effects $\theta_i$ and $\gamma_i$ are random variables but
are not observed.

In frequentist statistics, we *estimate* parameters and *predict*
unobserved random variables.  When we predict the random effects we
usually use an approach called *Best Linear Unbiased Prediction*
(BLUP).

The difference between estimation and prediction is subtle.  One way
of understanding this issue is that parameters determine the
distribution of every data value, and thus as we gain more and more
data, the parameter estimates will become arbitrarily precise (this is
*statistical consistency*).

Random effects generally relate to only a small number of
observations, which does not necessarily increase as our sample size
increases.  For example, we may have a longitudinal study in which we
have $n$ people observed at birth and at age five years.  Suppose we
use the model *y ~ x + (1 | b)* to describe this data.  We can let $n$
grow to recover the parameters of the model, which are $\beta_0$ (the
intercept), $\beta_1$ (the slope for $x$), $\tau^2$ (the variance of
the random intercepts), and $\sigma^2$ (the unexplained variation).
As $n$ grows, we will recover these four parameters to arbitrary
precision.  However if we use the BLUP $\hat{b}_i$ to predict the
random interept $b_i$, we will only ever have two observations for
case $i$ and therefore the prediction cannot be consistent for the
truth.

# Generalized linear multilevel regression

We may wish to use a GLM-like method to accommodate a non-linear single index
mean structure, and also accommodate non-indepence in the data.  Generalized linear
multilevel models (also known as GLIMMIX models) have been developed
for this purpose.

For a setting with a single set of blocks (i.e. a random intercept model), the mean structure of a GLIMMIX model is

$$
g(E[y_{ij} | x_{ij}, \eta_i]) = \beta^\prime x_{ij} + \eta_i
$$

where $\eta_i$ is a random effect with mean zero and variance $\tau^2$.

As with conventional GLMs, there is a mean/variance relationship, so the conditional
variance is

$$
{\rm var}[y_{ij} | x_{ij}, \eta_i] = \phi \cdot v(E[y_{ij} | x_{ij}, \eta_i])
$$

for a variance function $v: {\mathbb R}\rightarrow {\mathbb R}^+$ and a scale
parameter $\phi$.

GLIMMIX models do not enjoy the distributional robustness of conventional GLMs.  For
example, if we model the data as Poisson (given random effects), the data need to
actually be Poisson for most theoretical guarantees to hold.  In contrast,
in a conventional GLM, obtaining valid results from a "Poisson model"
only requires that the conditional mean and conditional variance of the Poisson model
hold.

A major phenomenon that arises in the multilevel GLM setting with a non-linear
link function is that the marginal and conditional covariate effects differ.  That is,
in a GLIMMIX model $\beta_1$ is the change in expected response (on the scale of the link function)
for a unit chance in $x_1$, for a fixed value of the random effect $\eta$.  This
leads to an important distinction between "subject specific" and "population average" covariate
effects.

Another important aspect of GLIMMIX models is that in practice they
are quite difficult to fit using maximum likelihood methods.  While algorithms
and software exist for this purpose, fitting such models can be slow and
issues such as non-convergence can occur.

# Multilevel regression

Regression analysis aims to understand the conditional distribution of
a scalar response $y$ in relation to explanatory variables
$x \in {\mathbb R}^p$.  In many familiar forms of regression, we focus on the
marginal distribution $y\\,|\\,x$.  However it may be that we collect data
in such a way that different observed values of $y$ are not statistically
independent of each other, i.e.

$$
P(y_1, \ldots, y_n | x_1, \ldots, x_n) \ne \prod_i P(y_i | x_i).
$$

Multilevel regression is a means to understand the
conditional mean $E[y|x]$, conditional variance ${\rm var}(y | x)$,
conditional covariances ${\rm cov}(y_i, y_j | x_i, x_j)$, and other
forms of dependence among the observations.

The usual manner in which correlated data arise in practice is when
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
zero and variance $\sigma^2$.

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
unknown values to be estimated based on the data.  Instead, the random
effects are "latent" data that are not observed and are integrated out
of the model before estimating the parameters using a procedure such
as maximum likelihood estimation.

The parameters of the model discussed above are $\beta$, $\tau^2$, and
$\sigma^2$.  All these parameters are estimated jointly, usually via
maximum likelihood estimation (MLE) or the closely-related restricted
maximum likelihood estimation (REML) which we will not define here.

### Using OLS to model multilevel data

It would be possible to fit a model to correlated data using OLS,
ignoring the correlations. Let $\check{\beta}$ denote this estimator,
while $\hat{\beta}$ denotes the estimated mean parameters for the
mixed model.  Further, let $\check{\sigma}^2$ denote the MLE of
$\sigma^2$ for simple linear regression (essentially the sample
variance of the residuals).

It is a fact that $\check{\beta} \approx \hat{\beta}$.  That is, we
can still recover the mean parameters even when ignoring the
correlations present in the data.  Further, $\check{\sigma}^2$ will be
approximately equal to $\sigma^2 + \tau^2$, the total variance from
both random effects and unexplained variation.  Nevertheless, there
are at least two important reasons that OLS is not the best choice for
analyzing such data.  First, and most important, the OLS standard
errors for $\check{\beta}$ will be wrong, usually too small.  This
will lead to overconfidence about any findings.  Second, OLS will not
recover $\beta$ as efficiently as possible, where "efficiency" here
means that we get the most accurate estimate possible for a given
sample size.

## Model formulas and long form data

Most modern software for fitting multilevel models expects the data to
be in long form.  That is, instead of having the data indexed by two
subscripts $y_{ij}$ with $i$ indexing blocks and $j$ indexing
observations within blocks, in long form we have only a single index
for all observations in all blocks.  A second variable $b$ is used to
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
and variance $\tau^2$.  Further, we take $\gamma_i$ and $\theta_i$ to
be independent of each other.

The formula for a model with random slopes can be expressed as
*y ~ x1 + x2 + (1 + x1 | b)*.

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

In frequentist statistics, we can *estimate* parameters and *predict*
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

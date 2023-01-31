# Generalized Linear Models

Generalized Linear Models (GLMs) are a type of single-index regression
model that substantially extends the range of analyses that can be
meaningfully carried out compared to using conventional linear models.

A single-index model expresses the conditional mean function
$E[Y|X=x]$ through a single *linear predictor*, which is a linear
function of the covariates having the form:

$$
\beta_0 + \beta_1x_1 + \cdots + \beta_p x_p.
$$

A linear model is a specific type of GLM that models the conditional
mean function directly as the linear predictor:

$$
E[Y|X=x] = \beta_0 + \beta_1x_1 + \cdots + \beta_p x_p.
$$

A GLM relates the conditional mean to the linear predictor via a link
function $g: {\cal R}\rightarrow {\cal R}$:

$$
g(E[Y|X=x]) = \beta_0 + \beta_1x_1 + \cdots + \beta_p x_p.
$$

The link function should be invertible, so we can also write

$$
E[Y|X=x] = g^{-1}(\beta_0 + \beta_1x_1 + \cdots + \beta_px_p).
$$

__Relationship to linear models with a transformed response
variable:__ For a nonlinear function $g$, the value of $g(E[Y|X=x])$
is not equal to the value of $E[g(Y)|X=x]$.  While there is a
superficial similarity between fitting a GLM and fitting a linear
model with a transformed dependent variable $g(Y)$, these approaches
to regression can perform quite differently.

__Additive and multiplicative mean structures:__ By far the most
common link functions are the identify function $g(x)=x$ and the
logarithm function $g(x)=\log(x)$.  These give very different models
for how the covariates are related to the mean response.  Suppose we
have two covariates $x_1$ and $x_2$, leading to the linear predictor
$\beta_0 + \beta_1x_1 + \beta_2x_2$.  With an identity link, the
relationship between the covariates and the conditional mean is
additive.  If we compare two cases with the same value of $x_2$ and
whose values of $x_1$ differ by 1 unit, then the case with the greater
value of $x_1$ will be expected to have a response that is $\beta_1$
units greater compared to the case with the lower value of $x_1$.
Written concisely, we have:

$$
E[Y|X_1=x_1+1, X_2=x_2] - E[Y|X_1=x_1, X_2=x_2] = \beta_1
$$

If instead of a linear link we have a log link function, then the case
with the greater value of $x_1$ will be expected to have a response
that is $\exp(\beta_1)$ times greater than the case with the lower
value of $x_1$:

$$
E[Y|X_1=x_1+1, X_2=x_2] / E[Y|X_1=x_1, X_2=x_2] = \exp(\beta_1).
$$

__Variance functions:__ GLM's specify a condfitional variance function $V(\cdot)$
as well as a mean function.  For a basic GLM, the variance function
has the form

$$
{\rm Var}[Y|X=x] = \phi \cdot V(E[Y|X=x]),
$$

where $V(\cdot): {\cal R}\rightarrow R^+$ is a given function and
$\phi \ge 0$ is a *scale parameter*.  In words, the
variance is expressed through a "mean/variance relationship", i.e.
the variance can only depend on the covariates through the mean that
they imply.

__GLM's and parametric probability models:__ Many basic GLM's are
equivalent to using maximum likelihood analysis to fit a parametric
probability model to the data.  There is an alternative
"quasi-likelihood" approach to understanding GLMs that does not
emphasize likelihoods or probability models; we discuss this further
below.

Some examples of GLM's that are equivalent to parametric probability
models are:

* _Gaussian linear model_: In the standard (canonical) Gaussian GLM,
$E[Y|X=x]$ is the linear predictor, which means that $g$ is the
identity function, and the variance function is constant,
$V(\mu) = 1$.  The scale parameter $\phi$ is identical to the usual "error
variance" $\sigma^2$ in a linear model.  This is equivalent to using
maximum likelihood estimation to fit the model in which
$P(Y|X=x) = N(\beta_0 +\beta_1x_1 +\cdots \beta_px_p, \sigma^2)$,
where $N(\mu, \sigma^2)$
is the normal distribution with mean $\mu$ and variance $\sigma^2$.

* _Poisson log-linear model_: In the canonical Poisson GLM,
$\log(E[Y|X=x])$ is equal to the linear predictor, which means that
$g$ is the logarithm function.  The variance function is the identity
function $V(\mu) = \mu$.  The scale parameter is fixed at $\phi=1$.  This is
equivalent to using maximum likelihood estimation to fit the model in
which $P(Y|X=x)$ is parameterized as a Poisson PMF with mean
$\exp(\beta_0 + \beta_1x_1+\cdots+ \beta_px_p)$.

When viewed in terms of parameterized probability distributions, a GLM
can be seen as indexing an infinite family of distributions through
the covariate vector $x$.  That is, in a Poisson GLM, $P(Y|X=x)$ is
Poisson for every $x$, but with a different mean value (and variance)
in each case.

Note that a GLM describes a collection of conditional distributions,
$P(Y|X=x)$, for different values of $x$.  A GLM (like any regression
procedure) says very little about the marginal distribution $P(Y)$.
In a Poisson GLM, repeated observations of $Y$ at the same value of
$X$ will follow a Poisson distribution, but the marginal distribution
of $Y$, which describes the pooled set of all $Y$ values (taken at
different values of $X$) will not be Poisson.  Similarly, in a
Gaussian linear model, $Y$ values taken at the same $X$ are Gaussian,
but the marginal distribution of $Y$ is not Gaussian.

## Overview of different GLM families

As noted above, some GLMs are "inspired" by parameteric families of
probability distributions.  The Poisson and Gaussian GLMs are very
widely used, but there are many other useful GLMs that can be
specified through different choices of the family, link function, and
variance function.  In fact there are infinitely many possible GLMs.
We will discuss a few of the most prominent ones here.

Note that all of the approaches discussed below are suitable for
non-negative response variables.  The main GLM family that is used
with data that can take on both positive and negative values is the
Gaussian family.

Other than the Gaussian/linear model, all of the GLM's discussed
here most commonly use the logarithm as the
link function, so the mean structure model is
$E[Y|X=x] = \exp(\beta^\prime x)$, although alternative link functions are
possible, giving rise to different mean structures.

The _negative binomial_ GLM can be seen as an extension of the Poisson
GLM.  Its mean/variance relationship involves a shape parameter
$\alpha$, such that

$$
{\rm Var}[Y|X=x] = E[Y|X=x] + \alpha \cdot E[Y|X=x]^2.
$$

If $\alpha=0$, the mean/variance relationship is the same as in the
Poisson model.  If $\alpha > 0$, the variance increases faster with
the mean compared to the Poisson GLM.  The _quasi negative binomial_
model introduces a scale parameter to the mean/variance relationship,
which becomes

$$
{\rm Var}[Y|X=x] = \phi\cdot (E[Y|X=x] + \alpha \cdot E[Y|X=x]^2).
$$

The _Tweedie GLM_ interpolates the mean/variance relationship between
a Poisson and a negative binomial GLM.  It has the form

$$
{\rm Var}[Y|X=x] = \phi\cdot E[Y|X=x]^p
$$

where $1 \le p \le 2$ is a "power parameter".

A _Gamma GLM_ has a mean/variance relationship that is

$$
{\rm Var}[Y|X=x] = \phi \cdot E[Y|X=x]^2.
$$

Note that this implies that the coefficient of variation is

$$
{\rm CV}[Y|X=x] = {\rm SD}[Y|X=x]/E[Y|X=x] = \sqrt{\phi},
$$

which is constant.  By contrast, in a Poisson GLM, the CV is
$1/E[Y|X=x]^{1/2}$, meaning that the relative error is smaller when
the mean is larger.

An _inverse Gaussian GLM_ has

$$
{\rm Var}[Y|X=x] = \phi \cdot E[Y|X=x]^3
$$

as its mean/variance relationship.  The conditional variance in this
case grows very fast with the mean.

## Mean/variance relationships and "overdispersion"

For data sets that involve counting the number of times that a rare
event happens, it is natural to think in terms of a Poisson model.
For example, if we have the number of car accidents that occur at a
particular roadway intersection in a city during one year, we might
imagine that this count is Poisson-distributed.  Recall that a Poisson
distribution describes the number of times that an event occurs when
the event has a constant probability of happening in each small time
interval, and if these occurrences are independent.  If there is a
fixed, small probability of a car accident occurring at a particular
intersection on each day, then the number of accidents at that
intersection per year might follow a Poisson distribution, as long as
the accidents on different days are independent.  If all the roadway
intersections in a city have the same probability of having an
accident per day, then the set of counts for all such intersections
may behave like a sample from a Poisson distribution.

If either the independence or homogeneity (equal probability)
conditions for a Poisson process stated above do not hold, then the
resulting counts may not be Poisson-distributed.  For example, if
there are differences in risk for different days of the week, or
during different seasons, then the total number of events in one year
will generally not follow a Poisson distribution.  For example, even
if the count for each intersection is Poisson-distributed, the
collection of counts for all intersections will not be a sample from a
Poisson distribution unless the daily probabilities of an accident
occurring at the different roadway intersections are the same.

To apply Poisson regression, we do not need the data to be exactly
Poisson-distributed.  We only need the specified mean and variance
structures to hold.  In this case the model is being fit from a
"quasi-likelihood" perspective, and this does not make the findings
any less meaningful than if the underlying data were truly Poisson.

A key property of a Poisson distribution is that the mean is equal to
the variance.  "Overdispersion" occurs when the variance is greater
than the mean (underdispersion is also possible but we will not
consider that here).  Overdispersion typically results from
"heterogeneity" -- the total count being modeled is a sum of counts
over temporal or spatial sub-units with different event rates.

## Scale parameters

The variance function $V(\mu)$, where $\mu = \mu(x)$ is the mean
structure, describes how the mean and variance are related in a GLM.
The value of $V(\cdot)$ must be either explicitly or implicitly
specified when fitting a GLM.  However it turns out that $V(\cdot)$
only needs to be the correct variance function up to an unknown
constant multiplier.  That is, as long as the true conditional
variance can be represented in the form
$\textrm{Var}(Y|X=x) = \phi\cdot V(\mu)$, where $\phi \ge 0$
is a constant, then we can still
estimate $\beta$ with good theoretical guarantees.

When working with binary data, the variance is completely determined
by the mean and the scale parameter must be set to $\phi=1$.

### Scale parameter estimation

One challenge to understanding the variance structure is that, as
noted above, a GLM defines a probability model for each value of the
covariate $x$.  That is, we observe a value $y$ drawn from the
probability model defined by the GLM at $X=x$.

It is the conditional distribution of $Y$ at a single value of $x$,
not the marginal distribution of $Y$, that matches the "family"
(Poisson, Gaussian, etc.).  Unless there are ties (replicate
observations at one $x$), we only observe one $y$ for each $x$, so we
only observe one $y$ from each distribution.  With a sample of size 1
(at each $x$), it is not straightforward to assess the shape of the
distribution or even its variance.

Fortunately, there is a way to estimate the scale parameter $\phi$
even when there are no replicates.  The most basic way of doing this
is using the following estimator that pools information across the
whole sample, adjusting for both the mean structure and the unscaled
variance function $V(\cdot)$

$$
\hat{\phi} = (n-p)^{-1}\sum_i (y_i - \hat{\mu}_i)^2 / V(\hat{\mu}_i)
$$

where $\hat{\mu}_i = g^{-1}(\hat{\beta}^\prime x_i)$ is the estimated
mean value for one observation, $n$ is the sample size, and $p$ is the
number of covariates.

### Robust estimation of the scale parameter

In some settings, we know or suspect that our data are contaminated
with outliers, but we do not wish to come up with a rule for deciding
which observations are problematic.  In this, case we can use a
so-called _robust_ estimate of the scale parameter that is less
impacted by outliers than the standard estimate.  Note that this does
not address the impact of the outliers on estimation of $\beta$, it
only addresses their impact on estimation of $\phi$.  In some
situations, this is sufficient.

An effective way to robustly estimate the scale parameter is _Huber's
proposal 2_.  If $r_i = (y_i - \hat{\mu}_i)/V(\hat{\mu_i})^{1/2}$ is
the Pearson standardized residual for observation $i$, the standard
estimate of the scale parameter solves the equation

$$
\sum_i (r_i/\hat{\phi})^2 = n-p,
$$

where $n$ is the sample size and $p$ is the number of covariates.
Huber's proposal 2 modifies this estimating equation to

$$
\sum_i h_c(r_i/\hat{\phi})^2 = k\cdot (n-p),
$$

where
$h_c(x) = x\cdot I_{|x| \le c} + c\cdot \textrm{sign}(x)\cdot I_{|x| > c}$.
The value of $c$ is a tuning parameter, often set to
$c=1.345$.  The value of $k$ is set to a value so that if the $r_i$
are Gaussian, the estimate of $\hat{\phi}$ coincides with the
population variance.

### Variance diagnostics

If we have enough data, it is possible to assess the goodness of fit
of the working variance model directly.  Taking the Poisson case
as an example, we can fit a preliminary Poisson regression, then
stratify the data into bins based on the fitted values $\hat{y}$ from
this model.  We then estimate the scale parameter within individual
bins.  If the variance model is correct, these
estimated scale parameters should be approximately constant with
respect to the estimated means.

## GLMs via quasi-likelihood

Just as most of the favorable properties of linear regression do not
depend on Gaussianity of the data, most of the favorable properties of
GLMs do not require the data to follow the specific family used to
define the model.  For example, Poisson regression can often be used
even if the data do not follow a Poisson distribution.  The primary
requirement for attaining good performance with a GLM is that the mean
structure is (approximately) correct.  An important secondary
requirement is to have a reasonably accurate model for the conditional
variance.  Other properties of the distribution (other than the
conditional mean and variance) are less important. For example,
domain constraints can generally be ignored, so a Poisson GLM can be
fit to data that includes non-integer values, as long as the mean and
variance models hold.

One useful example of a GLM fit using quasi-likelihood is
"quasi-Poisson" regression, which results from using Poisson
regression, but allowing the scale parameter $\phi$ to take on values
other than 1. In a quasi-Poisson regression, the variance is equal to
$\phi$ times the mean, i.e.
${\rm Var}(Y|X=x) = \phi\cdot E(Y|X=x) = \phi\exp(\beta^\prime x)$.
Conducting Poisson regression in a setting
where we suspect that the variance structure
${\rm Var}(Y|X=x) = \phi\cdot E(Y|X=x)$
holds for $\phi \ne 1$ is an appropriate way to
estimate the mean structure parameters $\beta$, but it is not maximum
likelihood.  This procedure is justified under the quasi-likelihood
theory for GLMs.

If we believe that the specified mean structure, and the variance
model are approximately correct, then quasi-likelihood theory
justifies that we (i) use standard GLM fitting procedures to estimate
$\beta$, ignoring the scale parameter $\phi$ at the time of fitting,
and (ii) adjust the standard errors for the regression parameter
estimates by multiplying them by a factor of $\hat{\phi}^{1/2}$. Some
other tools from maximum likelihood analysis can be adapted to the
quasi-likelihood setting, such as using the quasi-AIC instead of the
standard AIC.  Note that some ideas from maximum likelihood theory
such as likelihood ratio testing are not directly applicable in the
quasi-likelihood setting.

### Estimation without valid inference

An important robustness property of quasi-GLMs is that the mean
parameters can be estimated consistently even when the variance is
mis-specified.  Incorrect specification of the variance leads to a
loss of efficiency, but does not completely invalidate the analysis.
For example, we can use quasi-Poisson regression with a _working
variance structure_ ${\rm Var}(Y|X=x) = \phi\exp(\beta^\prime x)$ to
estimate the mean structure even if we are not confident that this
working variance structure is correct, as long as we believe that the
mean structure (e.g. based on the log link function) is correct.  In
this setting, we can consistently estimate the mean structure
parameters $\beta$, but most tools for statistical inference
(e.g. standard errors) may be incorrect.  This approach can be useful
if we are willing to sacrifice a (usually small) amount of efficiency
for the sake of simplicity, and do not need standard errors (e.g. if
our research aims are primarily predictive).

### Non-standard variance functions

There are canonical choices of the variance function for basic GLMs
that are modeled on parametric probability models.  For example,
$V(\mu) = 1$ matches the Gaussian distribution, and $V(\mu) = \mu$
matches the Poisson distribution, where $\mu = E[Y|X=x]$.  However it
is entirely legitimate to specify a non-standard variance function if
that provides a better fit to the data.  For example, we can fit a
Poisson regression model with variance function $V(\mu) = \mu^p$, for
a given value of $p$ suggested by the data.  This is a quasi-GLM,
since the fitting process is not equivalent to maximum likelihood
estimation for a model with the specified mean and variance
structures.

### Robust inference for GLMs

Adjusting the standard errors by multiplying with a factor of
$\hat{\phi}^{1/2}$ is appropriate if the variance model
$\textrm{Var}(Y|X=x) = \phi V(E[Y|X=x])$ is correct.  If we use a
working variance structure that may not be correct, then we should use
a "robust" approach to obtaining the standard errors, sometimes called
the "Huber-White" or "heteroscedasticity robust" approach to
inference.  We will not present the details of this approach here, but
it is supported by most software.

### Model selection

AIC and other information criteria (BIC, etc.) are very popular tools
for model selection, but in their standard implementation can only be
used for models fit with maximum likelihood.  Since quasi-likelihood
estimates are generally not maximum likelihood estimates, AIC cannot
be used in this setting.  As an alternative, there is a model
selection statistic called _QAIC_, which is defined as

$$
\textrm{QAIC}_j = -2\cdot Q_j/\tilde{\phi} + 2p_j,
$$

where $j=1,\ldots, m$ indexes models that we wish to compared, $Q_j$
is the quasilikelihood for model $j$, calculated using a scale
parameter of $\phi=1$, $\tilde{\phi}$ is the estimated scale parameter
for the parent model of all $m$ models being compared, and $p_j$ is
the number of parameters for model $j$.

The need to have a common parent model makes the QAIC a bit more
difficult to work with than the AIC, which can be computed for one
model without reference to what it will be compared to.

### Estimation process

GLMs are fit by solving the following system of _estimating
equations_, where $\mu_i = \mu_i(\beta) = g^{-1}(x_i^\prime \beta)$
is the mean of observation $i$ implied by a given value of $\beta$:

$$
\sum_i \partial \mu_i/\partial \beta \cdot (y_i - \mu_i) / V(\mu_i) = 0.
$$

Note that $\partial \mu_i/\partial \beta$ is a $p$-dimensional vector,
where $p$ is the dimension of $\beta$, so this is a system of $p$
equations in $p$ unknowns.  Since this is a non-linear system, there
is no guarantee that a solution exists or is unique for finite
samples, but the quasi-likelihood theory of GLMs gives conditions
under which there exists a sequence of solutions that are consistent
for the true mean structure parameters $\beta$.

In some cases, this estimating equation corresponds to a _score
function_, meaning that it is the gradient of a proper log-likelihood.
In particular, this score equation can be defined by considering
an exponential family of distributions, and constructing the
score function as the derivative of the log-likelihood.  This corresponds
to a maximum-likelihood approach to estimation.

An alternative way of deriving the score equations is based on the
method of moments.  Suppose we wish to estimate the mean parameters
$\beta$ in a model with a given mean function $\mu_i(\beta) = g^{-1}(x_i^\prime \beta)$
and a given variance function $V(\mu_i)$.  The inverse variance weighted
sum of squared residuals is

$$
\sum_i (y_i - \mu_i)^2/V(\mu_i).
$$

If we ignore the fact that $V(\cdot)$ depends on $\beta$ (imagine
that when computing $V(\mu_i)$ we plug-in an estimate of $\beta$
instead of including it in the optimization),
then the derivative of the inverse variance weighted sum of squared
residuals is (up to a multiplicative constant) the score equation
given above.  In this way, we can derive the score equation without
ever specifying a log-likelihood.  In doing so, we are choosing to measure
goodness of fit with the squared residuals, which corresponds to
maximum likelihood in the Gaussian setting but not for other distributions.
Nevertheless, using squared residuals to measure fit can be shown
to produce good estimators for a wide variety of distributions.  In addition,
this score equation accounts for heteroscedasticity by inverse variance
weighting with respect to $V(\cdot)$, and also considers the curvature
of the mean function through the Jacobian $\partial\mu/\partial\beta$.

However in some other cases, there is no log-likelihood for which this
expression is the gradient.  This is the reason that quasi-likelihood
analysis is not the same as maximum likelihood analysis, in general.

An important contribution of RWM Wedderburn is that the function

$$
\int_0^\mu (y - u)/V(u) du
$$

is the antiderivative of
$\partial \mu/\partial \beta \cdot (y - \mu) / V(\mu)$.  Therefore,
this expression can be seen as providing a
concrete quasi-likelihood that exhibits the mean and variance
structures we have chosen to fit to our data.  This has been used to
produce quasi-likelihood counterparts to important quantities from
likelihood analysis, including AIC, score testing, and log-likelihood
ratio testing.  Note that the AIC derived from Wedderburn's
quasi-likelihood is called QIC, and is generally different from the
QAIC discussed above.

## Generalized estimating equations

Many datasets contain data that are statistically dependent.  In some
settings, this dependence can be ignored, but in most cases it is
important to account for it in a data analysis.  There are many
approaches to working with dependent data.  Here we focus on a
framework known as "Generalized Estimating Equations" (GEE), which
extends the GLM approach to accommodate dependent data.

First, we present two examples of data that would likely be dependent:

* Suppose that we measure a person's C-reactive protein (CRP) levels,
a biomarker for inflammation.  A longitudinal study design might have
each subject return for annual assessments.  The repeated measures
taken on one subject would likely be statistically dependent.  This
form of dependence is called _serial dependence_ because it results
from having repeated measures taken over time.

* Suppose that we consider the number of COVID-19 deaths per day in
each US county, over the span of several months.  These counts could
be statistically dependent over time (serial dependence), and could
also be dependent within states.  That is,
all counties within a state may move up or down together.

Note that the data being dependent is mostly a property of the way in
which the data were collected, rather than being intrinsic to one type
of measurement.  For example, CRP data would be independent if
collected in a cross-sectional study, but is dependent when collected
in a longitudinal study.

Basic regression analysis focuses on the conditional mean relationship
$E[Y|X=x]$.  Modeling
multidimensional probability distributions is a much harder task.  The
goal of GEE, like GLM, is to focus on the conditional mean
relationship, while accommodating the presence of statistical
dependence.

Above we used working variance models and quasi-likelihood analysis to
justify the use of GLM fitting procedures in settings where the GLM
probability model may not hold.  We can now extend this idea to
develop a GLM-like procedure to accommodate non-independent data.

The key idea here is that of a _working correlation structure_.  This
is a model for how observations in the dataset are related.  Most GEE
software arranges the data into _clusters_, or _groups_.  Two
observations in different groups are always independent.  Two
observations in the same group may be dependent.  The working
correlation structure is an attempt to specify how these dependencies
are structured.  Formally, $R_i(\alpha)$ is the $n_i\times n_i$ working
correlation matrix for cluster $i$.  The vector $\alpha$ contains parameters
that determine the correlation structure.

The mean structure parameters $\beta$ are estimated by solving the
following estimating equations

$$
\sum_i \partial \mu_i/\partial \beta \cdot V(\mu_i)^{-1/2}R_i(\alpha)^{-1}V(\mu_i)^{-1/2}(y_i - \mu_i) = 0.
$$

The working correlation structure $R_i(\alpha)$
does not need to be correct in order to obtain meaningful results.
The mean structure parameters $\beta$ can be estimated, and standard
errors can be obtained, even if the working correlation structure is
wrong.  However the estimate of $\beta$ will be more efficient if the
working correlation structure is approximately correct.

There are many possible working correlation structures, but these are
some of the most common:

* _Independence_: In the independence working correlation structure,
all observations are taken to be independent, both within and between
clusters, so $R_i(\alpha) = I_{n_i\times n_i}$.  The parameter
estimates $\beta$ will be the same as
obtained in a GLM fit, but the standard errors will be different.

* _Exchangeable_: In an exchangeable working correlation structure,
any two observations in the same cluster are correlated at level
$\rho$, which is a parameter to be estimated from the data.  The
correlation model is
$R_i(\alpha) = \alpha\cdot 1_{n_i\times n_i} + (1-\alpha)\cdot I_{n_i\times n_i}$.

* _Autoregressive_: In an autoregressive working correlation
structure, the observations within a cluster are ordered, and the
correlation between observations $j$ and $k$ within cluster $i$
is $R_i(j,k) = \alpha^{|j-k|}$, where
$\alpha$ is a parameter to be estimated from the data.

* _Stationary_: In a stationary working correlation structure, the
observations within a cluster are ordered, and the correlation between
cases $i$ and $j$ is $\rho_{|i-j|}$ if $|i-j|\le m$, and the two cases
are uncorrelated otherwise. The values of $\alpha=(\rho_1, \ldots, \rho_m)$ are
estimated from the data, and $m$ is a tuning parameter.

The algorithm for fitting a GEE alternates between updating the
estimate of $\beta$, and updating the estimates of any correlation
parameters $\alpha$.  It is a quasi-likelihood approach, so has
many robustness properties such as being robust to mis-specification
of the variance structure, and of the correlation structure.  Note
that it remains a requirement for observations in different groups to
be independent.

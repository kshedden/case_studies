# Generalized Linear Models

*Generalized Linear Models* (GLMs) are a type of single-index regression
model that substantially extends the range of analyses that can be
meaningfully carried out compared to using conventional linear models.

A single-index model expresses the conditional mean function
$E[Y|X=x]$ through a single *linear predictor*, which is a linear
function of the covariates having the form:

$$
\beta_0 + \beta_1x_1 + \cdots + \beta_p x_p.
$$

A *linear model* is a specific type of GLM that models the conditional
mean function directly as the linear predictor:

$$
E[Y|X=x] = \beta_0 + \beta_1x_1 + \cdots + \beta_p x_p.
$$

A *generalized linear model* relates the conditional mean to the linear predictor via a link
function $g: {\mathbb R}\rightarrow {\mathbb R}$:

$$
g(E[Y|X=x]) = \beta_0 + \beta_1x_1 + \cdots + \beta_p x_p.
$$

The link function should be invertible, so we can also write

$$
E[Y|X=x] = g^{-1}(\beta_0 + \beta_1x_1 + \cdots + \beta_px_p).
$$

__Relationship to linear models with a transformed response
variable:__ For a nonlinear function $g$, the value of $g(E[Y|X=x])$
is not equal to the value of $E[g(Y)|X=x]$.  Thus, while there is a
superficial similarity between fitting a GLM and fitting a linear
model with a transformed dependent variable $g(Y)$, these approaches
to regression can perform quite differently.

__Additive and multiplicative mean structures:__ By far the most
common link functions are the identify function $g(x)=x$ and the
logarithm function $g(x)=\log(x)$.  These give very different
*mean structure models*
for how the covariates are related to the mean response.  Suppose we
have two covariates $x_1$ and $x_2$, leading to the linear predictor
$\beta_0 + \beta_1x_1 + \beta_2x_2$.  With an identity link, the
relationship between the covariates and the conditional mean is
additive.  If we compare two cases with the same value of $x_2$ and
whose values of $x_1$ differ by 1 unit, then the case with the greater
value of $x_1$ is expected to have a response that is $\beta_1$
units greater compared to the case with the lower value of $x_1$,
that is

$$
E[Y|X_1=x_1+1, X_2=x_2] - E[Y|X_1=x_1, X_2=x_2] = \beta_1.
$$

The covariate effect in this case is *additive*.

If instead of a linear link we have a log link function, then the case
with the greater value of $x_1$ is expected to have a response
that is $\exp(\beta_1)$ times greater than the case with the lower
value of $x_1$:

$$
E[Y|X_1=x_1+1, X_2=x_2] / E[Y|X_1=x_1, X_2=x_2] = \exp(\beta_1).
$$

The covariate effect in this case is *multiplicative*.

__Variance functions:__ GLM's specify a condfitional variance function $V(\cdot)$
as well as a mean function.  The conditional variance of the response is determined
from the variance function via the relationship

$$
{\rm Var}[Y|X=x] = \phi \cdot V(E[Y|X=x]),
$$

where $V(\cdot): {\mathbb R}\rightarrow {\mathbb R}^+$ is a given function and
$\phi \ge 0$ is a *scale parameter*.  In words, the
variance is expressed through a *mean/variance relationship*, i.e.
the conditional variance can only depend on the covariates through the conditional mean that
they imply.

Mean/variance relationships are a powerful way to accommodate *heteroscedasticity* in
regression analyses.  They provide a simpler alternative to *variance regression*
in which the conditional variance function ${\rm Var}[Y|X=x]$ is an arbitrary function
of $x$.

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
variance" $\sigma^2$ in a linear model.  Fitting a Gaussian linear model
as a GLM is equivalent to using maximum likelihood estimation to fit the model in which
$P(Y|X=x) = N(\beta_0 +\beta_1x_1 +\cdots \beta_px_p, \sigma^2)$,
where $N(\mu, \sigma^2)$
is the normal distribution with mean $\mu$ and variance $\sigma^2$.

* _Poisson log-linear model_: In the canonical Poisson GLM,
$\log(E[Y|X=x])$ is equal to the linear predictor, which means that
$g$ is the logarithm function.  The variance function is the identity
function $V(\mu) = \mu$.  The scale parameter is fixed at $\phi=1$.
Fitting the Poisson log-linear model as a GLM is
equivalent to using maximum likelihood estimation to fit the model in
which $P(Y | X=x)$ is parameterized as a Poisson probability mass
function (PMF) with mean
$\exp(\beta_0 + \beta_1x_1+\cdots+ \beta_px_p)$.

A GLM can be seen as indexing a collection of distributions through
the covariate vector $x$.  If we have observations
$(y_i, x_i), i=1, \ldots, n$, the conditional distribution $y_i | x_i$
of each observation follows a specific probability
distribution in the collection, based on $x_i$.
For example, in a Poisson GLM, $P(Y_i | X=x_i)$ is a Poisson distribution,
and different observations in general will follow different
Poisson distributions.

Here we express GLMs through expressions for the conditional mean structure
$E[Y|X=x]$ and conditional variance structure ${\rm Var}[Y | X=x]$.
Except for the Gaussian case, there is no
meaningful way to write a GLM in the "generative form"
$y_i = \beta^\prime x_i + \epsilon_i$ for an *error term* $\epsilon_i$.

Note that a GLM describes the conditional distributions
$P(Y | X=x)$, for different values of $x$.  A GLM (like any regression
procedure) says very little about the marginal distribution $P(Y)$.
In a Poisson GLM, repeated observations of $Y$ at the same value of
$X$ will follow a Poisson distribution, but the *pooled* (or *marginal*,
or *unconditional*)
distribution of $Y$, which describes the all $Y$ values taken at
different values of $X$, will not be Poisson.  Similarly, in a
Gaussian linear model, $Y$ values taken at the same $X$ are Gaussian,
but the marginal distribution of $Y$ is not Gaussian.

## Overview of different GLM families

Some GLMs such as the Gaussian, Poisson, and binomial GLMs are derived
from parameteric families of probability distributions.  We discuss
these familiar GLMs here.

The Gaussian GLM typically uses the identity link funtion $g(x) = x$,
the Poisson GLM typically uses the log link function $g(x) = \log(x)$,
and the binomial GLM typically uses the logistic link function
$g(x) = {\rm log}(x / (1-x))$.
It is possible to use alternative link functions with any of these families,
e.g. we can use the log link function for a Gaussian GLM, but we do not consider
this possibility further here.

All of the examples introduced below use the logarithm as the
link function, so the mean structure model is
$E[Y | X=x] = \exp(\beta^\prime x)$.

The _quasi Poisson_ GLM introduces a non-negative scale parameter $\phi$
into the Poisson variance model (in which $\phi=1$ is fixed):

$$
{\rm Var}[Y | X=x] = \phi E[Y | X=x].
$$

The _negative binomial_ GLM introduces a shape parameter
$\alpha$ into the variance model, such that

$$
{\rm Var}[Y | X=x] = E[Y | X=x] + \alpha \cdot E[Y | X=x]^2.
$$

The _quasi negative binomial_ model has both scale and shape parameters, so that
the mean/variance relationship is

$$
{\rm Var}[Y|X=x] = \phi\cdot (E[Y|X=x] + \alpha \cdot E[Y|X=x]^2).
$$

In the negative binomial models, if $\alpha=0$, the mean/variance relationship
is the same as in the Poisson model or quasi-Poisson models.  If $\alpha > 0$,
the variance increases faster as a convex non-linear function of the mean.

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
$1 / E[Y|X=x]^{1/2}$, meaning that the conditional variance is
smaller in relative terms when the conditional mean is larger.

An _inverse Gaussian GLM_ has

$$
{\rm Var}[Y|X=x] = \phi \cdot E[Y|X=x]^3
$$

as its mean/variance relationship.  The conditional variance in this
case grows very fast with the mean.

## Mean/variance relationships and "overdispersion"

For data sets that involve counting the number of times that an
event happens, it is natural to think in terms of a Poisson model.
For example, if we have the number of car accidents that occur at a
particular roadway intersection in a city during one year, we might
imagine that this count is Poisson-distributed.  Recall that a Poisson
distribution describes the number of times that an event occurs over a
long period of time when
the event has a constant probability of happening in each small time
interval, and if these occurrences are independent.  If there is a
fixed, small probability of a car accident occurring at a particular
intersection on each day, then the number of accidents at that
intersection per year might follow a Poisson distribution, as long as
the accidents on different days are independent.  If all the roadway
intersections in a city have the same probability of having an
accident per day, then the set of counts for all such intersections
may behave like a sample from a (common) Poisson distribution.

If either the independence or homogeneity (equal probability)
conditions for a Poisson process stated above do not hold, then the
resulting counts may not be Poisson-distributed.  For example, if
there are differences in risk for different days of the week, or
during different seasons of the year, then the total number of events in one year
will generally not follow a Poisson distribution.  For example, even
if the count for each intersection is Poisson-distributed, the
collection of counts for all intersections will not be a sample from a
Poisson distribution unless the daily probabilities of an accident
occurring at the different roadway intersections are the same.

To make use of Poisson regression, we do not need the data to be exactly
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
Overdispersion can be accommodated by using the quasi-Poisson,
negative binomial, quasi-negative binomial, Tweedie, Gamma, or
inverse Gaussian models as discussed above.

## Scale parameters

The variance function $V(\mu)$, where $\mu = \mu(x)$ is the mean
structure, describes how the mean and variance are related in a GLM.
The function $V(\cdot)$ must be
specified when fitting a GLM.  However it turns out that $V(\cdot)$
only needs to be the correct variance function up to an unknown
constant multiplier.  That is, as long as the true conditional
variance can be represented in the form
$\textrm{Var}(Y|X=x) = \phi\cdot V(\mu)$, where $\phi \ge 0$
is a constant, then we can still
estimate $\beta$ with good theoretical guarantees.

A special case arises when working with binary data, where the variance
is completely determined by the mean (i.e. there is only one possible
conditional distribution for a binary random variable once its mean
is known).  In this case, the scale parameter must be set to $\phi=1$.

### Scale parameter estimation

One challenge to understanding the variance structure is that, as
noted above, a GLM defines a probability model for each value of the
covariate $x$.  That is, we observe one value $y$ drawn from the
probability model defined by the GLM at $X=x$.

As noted above, it is the conditional distribution of $Y$ at a single value of $x$,
not the marginal distribution of $Y$, that matches the "family"
(Poisson, Gaussian, etc.).  Unless there are ties (replicate
observations at one $x$), we only observe one $y$ for each $x$, so we
only observe one $y$ from each distribution.  With a sample of size 1
(at each $x$), it is not straightforward to assess the shape of these
conditional distributions or even their variances.

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

It is possible to assess the goodness of fit of the variance model directly
by considering the *Pearson residuals*, which are defined to be
$r_i = (y_i - \hat{y}_i) / V(\hat{y}_i)^{1/2}$.  These quantities
should be approximately independent with mean zero and variance $\phi$.

A powerful way to detect miss-specification of the variance model is
to consider the relationship between the Pearson residuals $r_i$ and
the fitted values $\hat{y}_i$.  This is a generalization of the
"residuals on fitted values plot" used in basic multiple regression
analysis.  We can plot $r_i$ on $\hat{y}_i$ and visually assess
non constant variance, which if present would imply that the heteroscedasticity
in the data is not correctly captured by the chosen variance function.
A more focused way to detect heteroscedasticity is to
plot $|r_i|$ against $\hat{y}_i$ and use a scatterplot smoothing technique
like Lowess to estimate the
conditional mean relationship between these two
quantities, which should be constant.

## GLMs via quasi-likelihood

Just as most of the favorable properties of linear regression do not
depend on Gaussianity of the data, most of the favorable properties of
GLMs do not require the data to follow the specific family used to
define the model.  For example, Poisson regression can often be used
even if the data do not follow a Poisson distribution.  The primary
requirement for attaining good performance with a GLM is that the
specified mean and variance structures are correct.
Other properties of the distribution such as higher moments do not
need to match the specified parametric family. For example,
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
justifies that we may (i) use standard GLM fitting procedures to estimate
$\beta$, ignoring the scale parameter $\phi$ at the time of fitting,
and (ii) adjust the standard errors for the regression parameter
estimates by multiplying them by a factor of $\hat{\phi}^{1/2}$. Some
other tools from maximum likelihood analysis can be adapted to the
quasi-likelihood setting, such as using the quasi-AIC instead of the
standard AIC.  Note that other ideas from maximum likelihood theory
such as likelihood ratio testing are not directly applicable in the
quasi-likelihood setting.

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
\sum_i (\partial \mu_i/\partial \beta) \cdot (y_i - \mu_i) / V(\mu_i) = 0.
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

We can use this as a loss function for estimating $\beta$.
If we ignore the fact that $V(\cdot)$ depends on $\beta$ (imagine
that when computing $V(\mu_i)$ we plug-in an estimate of $\beta$
instead of including it in the optimization),
then the derivative of the inverse variance weighted sum of squared
residuals is (up to a multiplicative constant) the score equation
given above.  In this way, we can derive the score equations without
ever specifying a log-likelihood.  In doing so, we are choosing to measure
goodness of fit with the (weighted) sum of squared residuals, which corresponds to
maximum likelihood in the Gaussian setting but not for other distributions.
Nevertheless, using squared residuals to measure goodness of fit can be shown
to produce good estimators for a wide variety of distributions.  In addition,
this score equation accounts for heteroscedasticity by inverse variance
weighting with respect to $V(\cdot)$, and also considers the curvature
of the mean function through the Jacobian $\partial\mu/\partial\beta$.

Depending on how $\mu(\cdot)$ and $V(\cdot)$ are specified, there may or
may not be a log-likelihood for which the score equations derived above
are gradient.  For this reason, we may refer to these as *quasi-score equations*.
This is the reason that quasi-likelihood analysis is not in general the same as maximum
likelihood analysis.

An important contribution of RWM Wedderburn is that the function

$$
\int_0^\mu (y - u)/V(u) du
$$

is the antiderivative of
$(\partial \mu/\partial \beta) \cdot (y - \mu) / V(\mu)$.  Therefore,
this expression can be seen as providing a
quasi-likelihood function corresponding to the quasi-score function.
This quasi-likelihood may not be an actual likelihood, but it can be used to
produce quasi-likelihood counterparts to important quantities from
likelihood analysis, including AIC, score testing, and log-likelihood
ratio testing.  Note that the AIC derived from Wedderburn's
quasi-likelihood is called QIC, and is generally different from the
QAIC discussed above.

## Generalized Estimating Equations

Many datasets contain data that are statistically dependent, and
it is
important to account for this in data analysis.  There are many
approaches to working with dependent data.  Here we focus on a
framework known as *Generalized Estimating Equations* (GEE), which
extends the GLM approach to accommodate dependent data.

First, we present two examples of data that would likely be dependent:

* Suppose we have a collection of subjects $i=1, \ldots, n$ and repeatedly
measure a biomarker in each person's blood.  Let $y_{ij}$ denote the
$j^{\rm th}$ repeated measure of the biomarker for subject $i$, with this
measurement taking place at time $t_{ij}$, which does not need to
be the same for all subjects.  The repeated measures
taken on one subject would likely be statistically dependent.  This
form of dependence is called _serial dependence_ because it results
from having repeated measures taken over time.  In most such cases
we can assume that observations taken on different subjects are independent,
so we can partition our data into blocks such that there may
be dependence within blocks but not between blocks.

* Suppose that we consider the number of accidental deaths per day in
each US county, considering days ranging over a span of several months.
These counts could be statistically dependent within counties
over time (which is serial dependence as discussed above),
and moreover could
also be dependent between different counties within the same state.
This may be called _clustered dependence_.

Note that the data being dependent is mostly a property of the way in
which the data were collected, rather than being intrinsic to one type
of measurement.  For example, biomarker data would generally be independent if
collected in a cross-sectional study, but becomes dependent when collected
in a longitudinal study.

Modeling multidimensional probability distributions is very challenging.  The
goal of GEE, like GLM, is to focus on the conditional mean, but while
accommodating the presence of statistical dependence.

Above we used working variance models and quasi-likelihood analysis to
justify the use of GLM fitting procedures in settings where the GLM
probability model may not hold.  We can now extend this idea to
develop a GLM-like procedure to accommodate non-independent data.

The key idea here is that of a _working correlation structure_.  This
is a model for how observations in the dataset are thought to be
related.  Most GEE
software arranges the data into _clusters_ (also called _groups_ or
_blocks_).  Two
observations in different clusters are always independent, but two
observations in the same cluster may be dependent.  The working
correlation structure is an attempt to specify how these dependencies
are structured.  Formally, $R_i(\alpha)$ is the $n_i\times n_i$ working
correlation matrix for cluster $i$, where $n_i$ is the number of observations
in cluster $i$.  The vector $\alpha$ contains parameters
that determine the correlation structure (generally $\alpha$ has low dimension).

We can estimate the parameter vector $\beta$ by minimizing the *Mahalanobis distance*:

$$
\sum_i (y_i - \mu)_i^\prime \Sigma_i^{-1} (y_i - \mu_i),
$$

where $\Sigma_i = V_i^{1/2} R_i(\alpha) V_i^{1/2}$ is the covariance
matrix among observations on cluster $i$, with $V_i$ being
the diagonal $n_i\times n_i$ matrix with diagonal values equal to
$V(\mu_1), \ldots, V(\mu_{n_i})$.  Suppose for the moment that
$\Sigma_i$ is known, rather than depending on $\beta$ through
the $\mu_i$.  Then the value of $\beta$ that minimizes the Mahalanobis
distance solves the following *estimating equations*

$$
\sum_i (\partial \mu_i/\partial \beta)^\prime \cdot \Sigma_i \cdot (y_i - \mu_i) = 0.
$$

The Jacobian $J_i = \partial \mu_i/\partial \beta$ is a $n_i\times p$
matrix of partial derivatives and $y_i - \mu_i$ is a $n_i$-dimensional
vector of residuals.

The estimating equations can be solved using *Gauss-Seidel* iterations.
Essentially this involves using $\hat{\beta}$ from the previous
iteration to define $\Sigma_i$ and $J_i$, simplifying the
iterative update of $\beta$.  The resulting update is

$$
\hat{\beta} \leftarrow \hat{\beta} + \left(\sum_i J_i^\prime\Sigma_i^{-1}J_i\right)^{-1} \sum_i J_i^\prime \Sigma_i^{-1}(y_i - \hat{\mu}_i).
$$

The *robust* covariance matrix of $\hat{\beta}$ has the form $m^{-1}B^{-1}MB^{-1}$,
where $m$ is the number of clusters,

$$
B = m^{-1}\sum_i J_i^\prime\Sigma_i^{-1}J_i
$$

and

$$
M = m^{-1} \sum_i J_i^\prime \Sigma_i^{-1} r_ir_i^\prime \Sigma_i^{-1} J_i
$$

where $r_i = y_i - \hat{\mu}_i \in {\cal R}^{n_i}$ are working residuals.
This covariance matrix is robust to miss-specification of the working
correlation model $R_i(\alpha)$, however the estimate of $\beta$ will be
more efficient if the
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

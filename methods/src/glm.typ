#import "@preview/ilm:2.1.1": *

#show link: set text(fill: blue)

#set text(lang: "en")

#set math.equation(numbering: "1", supplement: none)

#show: ilm.with(
  title: [Generalized Linear Models],
  authors: "Kerby Shedden",
  figure-index: (enabled: true),
  table-index: (enabled: true),
  listing-index: (enabled: true),
  chapter-pagebreak: false
)

= Introduction

_Generalized Linear Models_ (GLMs) are a type of single-index regression model that substantially extends the range of analyses that can be meaningfully carried out compared to using conventional linear models (specifically, linear mean structure models fit using ordinary least squares).

A single-index model expresses the conditional mean function $E[Y|X=x]$ through a single _linear predictor_, which is a linear function of the covariates having the form:

$
beta_0 + beta_1x_1 + dots.h.c + beta_p x_p equiv eta.
$

A _linear model_ (or _linear mean structure model_) is a specific type of GLM that models the conditional mean function directly as the linear predictor:

$
E[Y|X=x] = beta_0 + beta_1x_1 + dots.h.c + beta_p x_p.
$

A _generalized linear model_ relates the conditional mean to the linear predictor via a link function $g: RR -> RR$:

$
g(E[Y|X=x]) = beta_0 + beta_1x_1 + dots.h.c + beta_p x_p.
$

The link function is invertible, so we can also write

$
E[Y|X=x] = g^(-1)(beta_0 + beta_1x_1 + dots.h.c + beta_p x_p) = g^(-1)(eta).
$

*Relationship to linear models with a transformed response variable*

For a nonlinear function $g$, the value of $g(E[Y|X=x])$ differs from the value of $E[g(Y)|X=x]$. Thus, while there is a superficial similarity between fitting a GLM with link function $g$, and fitting a linear model with a transformed dependent variable $g(Y)$, these approaches to regression can perform quite differently.

*Additive and multiplicative mean structures*

The most common link functions are the identify function $g(x)=x$ and the logarithm function $g(x) = log(x)$. These give very different _mean structure models_ for how the covariates are related to the mean response. Suppose we have two covariates $x_1$ and $x_2$, leading to the linear predictor $beta_0 + beta_1x_1 + beta_2x_2$. With an identity link, the relationship between the covariates and the conditional mean is additive. If we compare two cases with the same value of $x_2$ and whose values of $x_1$ differ by 1 unit, then the case with the greater value of $x_1$ is expected to have a response that is $beta_1$ units greater compared to the case with the lower value of
$x_1$, that is

$
E[Y|X_1=x_1+1, X_2=x_2] - E[Y|X_1=x_1, X_2=x_2] = beta_1.
$

If instead of a linear link we have a log link function, then the case with the greater value of $x_1$ is expected to have a response that is $exp(beta_1)$ times greater than the case with the lower value of $x_1$:

$
E[Y|X_1=x_1+1, X_2=x_2] / E[Y|X_1=x_1, X_2=x_2] = exp(beta_1).
$

The covariate effect in this case is _multiplicative_.

*Variance functions*

GLMs specify a conditional variance function $V(dot.c)$, such that the conditional variance of the response is determined from the variance function via the relationship

$
"Var"[Y|X=x] = phi dot.c V(E[Y|X=x]),
$

where $V(dot.c): RR -> RR^+$ is a given function and $phi >= 0$ is a _scale parameter_. In words, the variance is expressed through a _mean/variance relationship_, i.e. the conditional variance can only depend on the covariates through the conditional mean that they imply.

Mean/variance relationships are a powerful way to accommodate _heteroscedasticity_ in regression analyses. They provide a simpler alternative to _variance regression_ in which the conditional variance function $"Var"[Y|X=x]$ is an arbitrary function of $x$.

*GLMs and parametric probability models*

Many basic GLMs are equivalent to using maximum likelihood analysis to fit a parametric probability model to the data. This approach to GLMs can be developed through the use of exponential families.  However there is an alternative "quasi-likelihood" approach to understanding GLMs that does not emphasize likelihoods or exponential families; that is the approach taken here and developed further below.

Some examples of GLMs that are equivalent to parametric probability models
are:

- _Gaussian linear model_: In the Gaussian GLM, $E[Y|X=x]$ is the linear predictor, which means that $g$ is the identity function, and the variance function is constant, $V(mu) = 1$. The scale parameter $phi$ is identical to the usual "error variance" $sigma^2$ in a linear model. Fitting a Gaussian linear model as a GLM is equivalent to using maximum likelihood estimation to fit the model in which $P(Y|X=x) = N(beta_0 +beta_1 x_1 + dots.h.c beta_p x_p, sigma^2)$, where $N(mu, sigma^2)$ is the normal distribution with mean $mu$ and variance $sigma^2$.

- _Poisson log-linear model_: In the Poisson GLM, $log(E[Y|X=x])$ is equal to the linear predictor, which means that $g$ is the logarithm function. The variance function is the identity function $V(mu) = mu$. The scale parameter is fixed at $phi=1$. Fitting the Poisson log-linear model as a GLM is equivalent to using maximum likelihood estimation to fit the model in which $P(Y | X=x)$ is parameterized as a Poisson probability mass function (PMF) with mean $exp(beta_0 + beta_1 x_1 + dots.h.c + beta_p x_p)$.

A GLM can be seen as indexing a collection of distributions through the covariate vector $x$. If we have observations $(y_i, x_i), i=1, ..., n$, the conditional distribution $y_i | x_i$ of each observation follows a specific probability distribution in the collection, based on $x_i$. For example, in a Poisson GLM, $P(Y_i | X=x_i)$ is a Poisson distribution, and
different observations in general will follow different Poisson distributions.

Here we express GLMs through expressions for the conditional mean structure $E[Y|X=x]$ and conditional variance structure $"Var"[Y | X=x]$. Except for the Gaussian case, there is no meaningful way to write a GLM in the "generative form" $y_i = beta' x_i + epsilon_i$ for an _error term_ $epsilon_i$.

Note that a GLM describes the conditional distributions $P(Y | X=x)$, for different values of $x$. A GLM (like any regression procedure) says very little about the marginal distribution $P(Y)$. In a Poisson GLM, repeated observations of $Y$ at the same value of $X$ will follow a Poisson distribution, but the _pooled_ (or _marginal_, or _unconditional_) distribution of $Y$, which describes the all $Y$ values taken at different values of $X$, will not be Poisson. Similarly, in a Gaussian linear model, $Y$ values taken at the same $X$ are Gaussian, but the marginal distribution of $Y$ is not Gaussian.

= Overview of different GLM families

Some GLMs such as the Gaussian, Poisson, and binomial GLMs are derived from parameteric families of probability distributions. We discuss these familiar GLMs here.

The Gaussian GLM typically uses the identity link funtion $g(x) = x$, the Poisson GLM typically uses the log link function $g(x) = log(x)$, and the binomial GLM typically uses the logit link function $g(x) = log(frac(x, (1-x), style: "horizontal"))$. It is possible to use alternative link functions with any of these families, e.g. we can use the log link function for a Gaussian GLM or an identity link for a Poisson GLM, and there are sometimes good reasons for doing so.

All of the examples introduced below use the logarithm as the link function, so the mean structure model is $E[Y | X=x] = exp(beta' x)$.

The _Poisson_ GLM specifies $V(mu) = mu$, consistent with the Poisson distribution.  This gives the variance model

$
"Var"[Y | X=x] = E[Y | X=x],
$

where $phi = 1$.

The _quasi-Poisson_ GLM introduces a non-negative scale parameter $phi$ into the Poisson variance model, so that

$
"Var"[Y | X=x] = phi dot.c E[Y | X=x].
$

The _negative binomial_ GLM has both linear and quadratic terms in the variance model, and their combination is governed by a _shape parameter_ $alpha$, such that

$
"Var"[Y | X=x] = E[Y | X=x] + alpha dot.c E[Y | X=x]^2.
$

The _quasi negative binomial_ model has both scale and shape parameters, so that the mean/variance relationship is

$
"Var"[Y|X=x] = phi dot.c (E[Y|X=x] + alpha dot.c E[Y|X=x]^2).
$

In the negative binomial models, if $alpha=0$, the mean/variance relationship is the same as in the Poisson model or quasi-Poisson models. If $alpha > 0$, the variance increases faster as a convex non-linear function of the mean.

A _gamma GLM_ has a mean/variance relationship that is

$
"Var"[Y|X=x] = phi dot.c E[Y|X=x]^2.
$

Note that this implies that the coefficient of variation is

$
"CV"[Y|X=x] = "SD"[Y|X=x] / E[Y|X=x] = sqrt(phi),
$

which is constant. By contrast, in a Poisson GLM, the CV is $frac(1, E[Y|X=x]^(1/2), style: "horizontal")$, meaning that the conditional variance is smaller in relative terms when the conditional mean is larger.

An _inverse Gaussian GLM_ has

$
"Var"[Y|X=x] = phi dot.c E[Y|X=x]^3
$

as its mean/variance relationship. The conditional variance in this case grows
very fast with the mean.

The _Tweedie GLM_ interpolates the mean/variance relationship between a Poisson and a gamma GLM. It has the form

$
"Var"[Y|X=x] = phi dot.c E[Y|X=x]^p
$

where $1 <= p <= 2$ is a "power parameter".

= Mean/variance relationships and "overdispersion"

For data sets that involve counting the number of times that an event happens, it is natural to think in terms of a Poisson model. For example, if we have the number of car accidents that occur at a particular roadway intersection in a city during one year, we might imagine that this count is Poisson-distributed. Recall that a Poisson distribution describes the number of times that an event occurs over a long period of time when the event has a constant probability of happening in each small time interval, and if these occurrences are independent. If there is a fixed, small probability of a car accident occurring at a particular intersection on each day, then the number of accidents at that intersection per year might follow a Poisson distribution, as long as the accidents on different days are independent. If all the roadway intersections in a city have the same probability of having an accident per day, then the set of counts for all such intersections may behave like a sample from a (common) Poisson distribution.

If either the independence or homogeneity (equal probability) conditions for a Poisson process stated above do not hold, then the resulting counts may not be Poisson-distributed. For example, if there are differences in risk for different days of the week, or during different seasons of the year, then the total number of events in one year will generally not follow a Poisson distribution. For example, even if the count for each intersection is Poisson-distributed, the collection of counts for all intersections will not be a sample from a Poisson distribution unless the daily probabilities of an accident occurring at the different roadway intersections are the same.

To make use of Poisson regression, we do not need the data to be exactly Poisson-distributed. We only need the specified mean and variance structures to hold. In this case the model is being fit from a "quasi-likelihood" perspective, and this does not make the findings any less meaningful than if the underlying data were truly Poisson.

A key property of a Poisson distribution is that the mean is equal to the variance. "Overdispersion" occurs when the variance is greater than the mean (underdispersion is also possible but we will not consider that here). Overdispersion typically results from "heterogeneity" -- the total count being modeled is a sum of counts over temporal or spatial sub-units with different event rates. Overdispersion can be accommodated by using the quasi-Poisson, negative binomial, quasi-negative binomial, Tweedie, Gamma, or inverse Gaussian models as discussed above.

= Scale parameters

While the variance function $V(dot.c)$ must be specified when fitting a GLM, the conditional variance has the form $"Var"(Y|X=x) = phi dot.c V(mu)$, where the scale parameter $phi >= 0$ is a constant.  It turns out that estimation of $beta$ is orthogonal to estimation of $phi$.  That is, we do not need to worry about how to estimate $phi$ if our goal is to estimate $beta$.

A special case arises when working with binary data, where the variance is completely determined by the mean (i.e. there is only one possible conditional distribution for a binary random variable once its mean is known). In this case, the scale parameter must be set to $phi=1$.

== Scale parameter estimation

One challenge to understanding the variance structure is that, as noted above, a GLM defines a probability model for each value of the covariate $x$. That is, we observe one value $y$ drawn from the probability model defined by the GLM at $X=x$.

As noted above, it is the conditional distribution of $Y$ at a single value of $x$, not the marginal distribution of $Y$, that matches the "family" (Poisson, Gaussian, etc.). Unless there are ties (replicate observations at one $x$), we only observe one $y$ for each $x$, so we only observe one $y$ from each conditional distribution. With a sample of size 1 (at each $x$), it is not straightforward to assess the shape of these conditional distributions or even their
variances.

Fortunately, there is a way to estimate the scale parameter $phi$ even when there are no replicates. The most basic way of doing this is using the following estimator that pools information across the whole sample, adjusting for both the mean structure and the unscaled variance function $V(dot.c)$

$
hat(phi) = (n-p)^(-1) sum_i frac((y_i - hat(mu)_i)^2, V(hat(mu)_i), style: "horizontal")
$

where $hat(mu)_i = g^(-1)(hat(beta)' x_i)$ is the estimated mean
value for one observation, $n$ is the sample size, and $p$ is the number of
covariates.

== Robust estimation of the scale parameter

In some settings, we know or suspect that our data are contaminated with outliers, but we do not wish to come up with a rule for deciding which observations are problematic. In this, case we can use a so-called _robust_ estimate of the scale parameter that is less impacted by outliers than the standard estimate. Note that this does not address the impact of the outliers
on estimation of $beta$, it only addresses their impact on estimation of $phi$. In some situations, this is sufficient.

An effective way to robustly estimate the scale parameter is _Huber's proposal 2_. If $r_i = frac((y_i - hat(mu)_i), V(hat(mu)_i)^(1/2), style: "horizontal")$ is the Pearson standardized residual for observation $i$, the standard estimate of the scale parameter solves the equation

$
sum_i (frac(r_i, hat(phi), style: "horizontal"))^2 = n-p,
$

where $n$ is the sample size and $p$ is the number of covariates. Huber's proposal 2 modifies this estimating equation to

$
sum_i h_(c)(frac(r_i, hat(phi), style: "horizontal"))^2 = k dot.c (n-p),
$

where $h_c(x) = x dot.c bb(1)_(|x| <= c) + c dot.c "sign"(x) dot.c bb(1)_(|x| > c)$. The value of $c$ is a tuning parameter, often set to $c=1.345$. The value of $k$ is set to a value so that if the $r_i$ are Gaussian, the estimate of $hat(phi)$ coincides with the population variance.

= Variance diagnostics

It is possible to assess the goodness of fit of the variance model directly by considering the _Pearson residuals_, which are defined to be $r_i = frac((y_i - hat(y)_i), V(hat(y)_i)^(1/2), style: "horizontal")$. These quantities should be approximately independent with mean zero and variance $phi$.

A powerful way to detect miss-specification of the variance model is to consider the relationship between the Pearson residuals $r_i$ and the fitted values $hat(y)_i$. This is a generalization of the "residuals on fitted values plot" used in basic multiple regression analysis. We can plot $r_i$ on $hat(y)_i$ and visually assess non constant variance, which if present would imply that the heteroscedasticity in the data is not correctly captured by the
chosen variance function. A more focused way to detect heteroscedasticity is to plot $|r_i|$ against $hat(y)_i$ and use a scatterplot smoothing technique like Lowess to estimate the conditional mean relationship between these two quantities, which should be constant.

= GLM parameter estimation

While GLMs are often presented as a form of parametric regression model fit using maximum likelihood estimation, an alternative development based on _estimating equations_ does not refer to an explicit probability model.  This motivates us to seek a better understanding of what distributional properties are actually important for a GLM to be meaningful. In particular, since research aims often focus on the mean structure, we may ask:

- What orthogonality properties exist between the residuals $y_i - mu_i$ and the covariates?

- Can the estimates be viewed as minimizing the discrepancy between the observed data $y_i$ and the fitted means $mu_i$?

GLMs are fit by solving the system of _estimating equations_ or _score equations_

$
sum_(i=1)^n (frac(partial mu, partial beta, style: "horizontal"))(mu_i) dot.c frac((y_i - mu_i), V(mu_i), style: "horizontal") = 0_p,
$ <score1>

where $mu_i = mu(eta_i) = g^(-1)(x'_i beta)$ is the mean of observation $i$ implied by a given value of $beta$.  This can also be written as

$
sum_(i=1)^n X'_i dot.c (frac(d mu, d eta, style: "horizontal"))(eta_i) dot.c frac((y_i - mu_i), V(mu_i), style: "horizontal") = 0_p,
$ <score2>

where $eta_i = X'_i beta$ is the linear predictor.

Equivalently, we can write the score equations as

$
X'D'V^(-1)(y - mu),
$ <score4>

where $D$ is the $n times n$ diagonal matrix of the $(frac(d mu, d eta, style: "horizontal"))(eta_i)$, and $V$ is the $n times n$ diagonal matrix of the $V(mu_i)$.

Equation @score2 tells us that in a fitted GLM, the residuals are orthogonal to the fitted values when weighting by $frac((d mu) / (d eta)(eta_i), V(mu_i), style: "horizontal")$.  Since $mu(eta)$ is always chosen to be increasing in $eta$, these weights are non-negative.

Since $mu =& g^(-1)(eta)$, we have

$
(d mu) / (d eta) = (d g^(-1)(eta)) / (d eta)
= 1 / (g'(g^(-1)(eta))) = 1 / (g'(mu)).
$

Therefore, if $g'(mu) = frac(1, V(mu), style: "horizontal")$, the score equations (@score2) simplify greatly to

$
sum_(i=1)^n X'_i (y_i - mu_i) = 0_p,
$ <score3>

implying that the estimates are characterized by the property that the residuals are orthogonal to the covariates.  When $g'(mu) = frac(1, V(mu), style: "horizontal")$, the link is said to be _canonical_.  While there are some advantages to using the canonical link, it would be limiting to always require the link to be canonical, and arguably the weighted orthogonality in (@score2) may sometimes be more appropriate than the unweighted orthogonality in (@score3).

The parameter estimates obtained by solving (@score2) are equivalent to matching the data $y_i$ to the means $mu_i$ using a #link("https://en.wikipedia.org/wiki/Bregman_divergence")[Bregman divergence] loss function.  That is, $hat(beta)$ minimizes $sum_(i=1)^n D_(A^*)(y_i, mu_i)$, where $D_(A^*)(dot.c, dot.c)$ is a particular Bregman divergence.  The specific Bregman divergence is implied by the cumulant function $b$ such that $mu = b'$, $V = b''$, and so on, and $A^*$ is the #link("https://en.wikipedia.org/wiki/Convex_conjugate")[Legendre-Fenchel convex conjugate] of $b$.  At a high level this means that fitting a GLM involves matching the data to the fitted means, and the higher order properties of the assumed model only impact the estimates to the extent that they alter the loss function that is used for this mean matching.

Given that GLMs can be seen as a form of mean matching, not exclusively as maximum likelihood estimation of a fully and correctly-specified model, it is now unsuprising that the estimates obtained by solving (@score2) are consistent and asymptotically unbiased even if the variance model is incorrect.  Misspecification of the variance model can result in lower efficiency, and must be accounted for when doing statistical inference, but does not harm the fundamental validity of the estimates in terms of their being consistent and asymptotically unbiased.

== Quasi-likelihood

The _quasi-likelihood_ function

$
q(mu; y) = integral_0^mu (y - u)/V(u) "du"
$

satisfies

$
lr(d / (d mu) q(mu, y_i) |)_(mu=mu_i) = frac((y_i - mu_i), V(mu_i), style: "horizontal").
$

Now if we parameterize $mu_i = mu_(i)(beta)$ we get

$
lr(partial / (partial beta) q(mu(beta), y_i) |)_(mu=mu_(i)(beta)) = (frac(partial mu, partial beta, style: "horizontal"))(mu_i) dot.c frac((y_i - mu_i), V(mu_i), style: "horizontal"),
$

which recovers the GLM estimating equations (@score1).

The quasi-likelihood was first introduced by  #link("https://en.wikipedia.org/wiki/Robert_Wedderburn_(statistician)")[RWM Wedderburn].  It can be seen as providing an objective function for the estimating equations (@score1), that is available regardless of whether the link is canonical. This quasi-likelihood may not be an actual likelihood, but it can be used to produce quasi-likelihood counterparts to important quantities from likelihood analysis, including AIC, score testing, and log-likelihood ratio testing.

In summary, just as most of the favorable properties of linear regression do not depend on Gaussianity of the data, most of the favorable properties of GLMs do not require the data to follow the specific family used to define the model. For example, Poisson regression can be used even if the data do not (conditionally) follow a Poisson distribution. The primary requirement for attaining good performance with a GLM is that the specified mean and variance structures are correct. Other properties of the distribution such as higher moments do not need to match the specified parametric family. Furthermore, domain constraints can generally be ignored, so a Poisson GLM can be fit to data that includes non-integer values, as long as the mean and variance models hold.

= Inference for GLMs

The "naive" (model-based) variance covariance matrix for the GLM parameter estimates is

$
hat("cov")(hat(beta)) = hat(phi) dot.c (X' hat(W) X)^(-1),
$

where $hat(W) = "diag"(hat(w)_1, ..., hat(w)_n)$, with

$
hat(w)_i &= frac((lr((d mu) / (d eta) |)_(eta=eta_i))^2, V(hat(mu)_i), style: "horizontal")\
 &= frac(1, (V(hat(mu)_i) dot.c g'(hat(mu_i))^2), style: "horizontal").
$

In the case of a canonical link, this simplifies to $hat(w)_i = V(hat(mu)_i)$.  This approach works well as long as both the link function and the variance model are correctly specified.

== Robust inference

To allow the variance model to be misspecified, we can set $B = X' hat(W) X$, and $M = X' hat(W) hat(Sigma) hat(W) X$, where $hat(Sigma)$ is a diagonal matrix such that

$
hat(Sigma)_(i i) = [(y_i - mu_i) dot.c g'(mu_i)]^2
$

is the squared $i^"th"$ _working residual_.  The definition of the working residual is motivated by the linearization

$
g(y_i) - g(mu_i) approx g'(mu_i)(y_i - mu_i).
$

We then use

$
hat("cov")(hat(beta)) = B^(-1) M B^(-1).
$ <robust_cov>

The matrix $M$ can be constructed directly as $sum_(i=1)^n c_i^2 X_i X'_i$, where

$
c_i = (frac(d mu, d eta, style: "horizontal"))(eta_i) dot.c frac((y_i - mu_i), V(mu_i), style: "horizontal")
$ <score_coef>

is the contribution of observation $i$ to the estimating equations (@score1).

The covariance estimate (@robust_cov) is known as the "robust covariance estimate" or the "sandwich covariance estimate", and the matrices $B$ and $M$ are referred to as the "bread" and the "meat" respectively.

== Score testing

In the quasi-likelihood approach to GLMs, the primary interpretation of the estimates is through the score equations rather than through minimization of an objective function.  This leads to a useful technique for comparing models based on _score tests_.  Suppose we have a fitted GLM based on a design matrix $X$, and we have additional covariates in a second design matrix $Z$.  The null hypothesis is that $mu in "col"(X)$. We can test this null hypothesis based on score equations as follows.

Let $hat(mu)$ denote the fitted values under the null hypothesis, and let

$
Z_perp = Z - X(X' W X)^(-1)X' W Z
$

denote the design matrix $Z$ residualized with respect to $X$.  Then, define

$
U = Z'_perp D V^(-1)(y - hat(mu))
$ <score_test1>

as the score equations.  Recall that $X' D V^(-1)(y - hat(mu)) equiv 0$ always holds (equation @score4).  If there is no additional information in the $Z$ covariates, then (@score_test1) should be small.  A robust covariance matrix for $U$ is $Psi^(-1)$, where

$
Psi = sum_i c_i^2 z_(perp i) z'_(perp i),
$

with $c_i$ as in (@score_coef).

The score test statistic is $U' Psi^(-1) U$ and under the null hypothesis asymptotically has a $chi^2$ distribution with $q$ degrees of freedom, where $q$ is the number of columns of $Z$.

== Model selection

AIC and other information criteria (BIC, etc.) are very popular tools for model selection, but in their standard implementation can only be used for models fit with maximum likelihood. Since quasi-likelihood estimates are generally not maximum likelihood estimates, AIC cannot be used in this setting. As an alternative, there is a model selection statistic called _QAIC_, which is defined as

$
"QAIC"_j = -2 dot.c frac(Q_j, tilde(phi), style: "horizontal") + 2p_j,
$

where $j=1, ..., m$ indexes models that we wish to compared, $Q_j$ is the quasi-likelihood for model $j$, calculated using a scale parameter of $phi=1$, $tilde(phi)$ is the estimated scale parameter for the parent model of all $m$ models being compared, and $p_j$ is the number of parameters for model $j$.

The need to have a common parent model makes the QAIC a bit more difficult to
work with than the AIC, which can be computed for one model without reference
to what it will be compared to.

= Generalized Estimating Equations

Many studies involve data that are statistically dependent, and it is important to account for this in data analysis. There are many approaches to working with dependent data. Here we focus on a framework known as _Generalized Estimating Equations_ (GEE), which extends the GLM approach to accommodate dependent data.

First, we present two examples of data that would likely be dependent:

- Suppose we have a collection of subjects $i=1, ..., n$ and repeatedly measure a biomarker in each person's blood. Let $y_(i j)$ denote the $j^"th"$ repeated measure of the biomarker for subject $i$, with this measurement taking place at time $t_(i j)$, which does not need to be the same for all subjects. The repeated measures taken on one subject would likely be statistically dependent. This form of dependence is called _serial dependence_ because it results from having repeated measures taken over time. In most such cases we can assume that observations taken on different subjects are independent, so we can partition our data into blocks such that there may be dependence within blocks but not between blocks.

- Suppose that we consider the number of accidental deaths per day in each US county, considering days ranging over a span of several months. These counts could be statistically dependent within counties over time (which is serial dependence as discussed above), and moreover could also be dependent between different counties within the same state. This may be called _clustered dependence_.

Note that the data being dependent is mostly a property of the way in which the data were collected, rather than being intrinsic to one type of measurement. For example, biomarker data would generally be independent if collected in a cross-sectional study, but becomes dependent when collected in a longitudinal study.

Modeling multidimensional probability distributions is very challenging. The goal of GEE, like GLM, is to focus on the conditional mean, but while accommodating the presence of statistical dependence.

Above we used working variance models and quasi-likelihood analysis to justify the use of GLM fitting procedures in settings where the GLM probability model may not hold. We can now extend this idea to develop a GLM-like procedure to accommodate non-independent data.

The key idea here is that of a _working correlation structure_. This is a model for how observations in the dataset are thought to be related. Most GEE software arranges the data into _clusters_ (also called _groups_ or _blocks_). Two observations in different clusters are always independent, but two observations in the same cluster may be dependent. The working correlation structure is an attempt to specify how these dependencies are structured. Formally, $R_(i)(alpha)$ is the $n_i times n_i$ working correlation matrix for cluster $i$, where $n_i$ is the number of observations in cluster $i$. The vector $alpha$ contains parameters that determine the correlation structure (generally $alpha$ has low dimension).

We estimate the parameter vector $beta$ using score equations that are the natural multivariate generalization of (@score1)

$
sum_i (frac(partial mu_i, partial beta, style: "horizontal"))' dot.c Sigma_i^(-1) dot.c (y_i - mu_i) = 0,
$

and (@score2)

$
sum_i X'_i D'_i Sigma_i^(-1) (y_i - mu_i) = 0,
$

where $D_(i)(j, j) = lr((d mu) / (d eta) |)_(mu = mu_(i j))$ is a $n_i times n_i$ diagonal matrix of derivatives, and $y_i - mu_i$ is a $n_i$-dimensional vector of residuals.  The covariance matrix is expressed as $Sigma_i = V_i^(1/2)R_(i)(alpha)V_i^(1/2)$, where $V_i$ are diagonal matrices determined by the variance function, i.e. $V_i = "diag"(V(mu_(i 1)), ..., V(mu_(i n_i)))$. The matrix $R_(i)(alpha)$ is a "working correlation matrix" that does not need to be specified correctly.

The estimating equations can be solved using _Gauss-Seidel_ iterations. Let
$J_i = frac(partial mu_i, partial beta, style: "horizontal") = D_i X_i$ denote the Jacobian, an $n_i times p$ matrix of partial derivatives.  We use $hat(beta)$ from the previous iteration to define $Sigma_i$ and $J_i$. The resulting update is

$
hat(beta) -> hat(beta) + (sum_(i=1)^m J'_i Sigma_i^(-1)J_i)^(-1) sum_(i=1)^m J'_i Sigma_i^(-1)(y_i - hat(mu)_i).
$

The _robust_ covariance matrix of $hat(beta)$ has the form $B^(-1) M B^(-1)$ where,

$
B = sum_(i=1)^m J'_i Sigma_i^(-1)J_i
$

$
M = sum_(i=1)^m J'_i Sigma_i^(-1) r_i r'_i Sigma_i^(-1) J_i.
$

This covariance matrix is robust to miss-specification of the working correlation model $R_(i)(alpha)$, however the estimate of $beta$ will be more efficient if the working correlation structure is approximately correct.

There are many possible working correlation structures, but these are some of
the most common:

- _Independence_: In the independence working correlation structure, all observations are taken to be independent, both within and between clusters, so $R_(i)(alpha) = I_(n_i times n_i)$. The parameter estimates $beta$ will be the same as obtained in a GLM fit, but the standard errors will be different.

- _Exchangeable_: In an exchangeable working correlation structure, any two observations in the same cluster are correlated at level $rho$, which is a parameter to be estimated from the data. The correlation model is $R_(i)(alpha) = alpha dot.c bold(1)_(n_i times n_i) + (1-alpha) dot.c I_(n_i times n_i)$.

- _Autoregressive_: In an autoregressive working correlation structure, the observations within a cluster are ordered, and the correlation between observations $j$ and $k$ within cluster $i$ is $R_(i)(j,k) = alpha^(|j-k|)$, where $alpha$ is a parameter to be estimated from the data.

- _Stationary_: In a stationary working correlation structure, the observations within a cluster are ordered, and the correlation between cases $i$ and $j$ is $rho_(|i-j|)$ if $|i-j| <= m$, and the two cases are uncorrelated otherwise. The values of $alpha=(rho_1, ..., rho_m)$ are estimated from the data, and $m$ is a tuning parameter.

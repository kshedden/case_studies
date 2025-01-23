# Survival analysis

Survival analysis is a set of techniques for characterizing
[probability distributions](https://en.wikipedia.org/wiki/Probability_distribution)
and building
[statistical models](https://en.wikipedia.org/wiki/Statistical_model). Most
commonly, survival analysis methods are used with *time to event* or
*duration* data, where "time to event" refers to the duration of time from an
origin until some event of interest occurs. Many such examples consider
*failure times*, such as the duration of time from when a person is diagnosed
with a disease until they progress to a more advanced stage of the disease.
There are many other applications of survival analysis methods that have
nothing to do with "survival" or "failure", and these methods can even be
applied in settings where time does not play a central role (although this is
less common).

## Time origin

Consider a setting where we are observing a person over time, and we are
interested in the time $T$ at which some event of interest occurs. It is very
important to define the *origin* from which time is measured, i.e. what is the
meaning of time zero? In the case of human-centered studies, we may let $T$
denote the age of the person when the event occurs, in which case the time
origin is the date of birth. Note that in this case you will need to define
[units](https://en.wikipedia.org/wiki/Unit_of_measurement) for $T$, e.g.
months or years. Alternatively, there may be an event that must occur before
the event of interest, and the time of this event may make a more sensible
time origin. For example, if $T$ corresponds to graduating from university, we
may choose the date of university matriculation (enrollment) as the time
origin, so that e.g. if $T=4$ and the time units are years, then the person
graduated four years after beginning their studies.

## Event time distributions

The theoretical basis of conventional survival analysis is that we are
studying the probability distribution of a
[random variable](https://en.wikipedia.org/wiki/Random_variable) $T$,
corresponding to the time at which an event of interest occurs. In
conventional survival analysis, the event will always occur if we wait long
enough, so $P(T < \infty) = 1$. We note that in some situations this may not
be completely realistic and there is a subdomain of survival analysis called
"cure modeling" in which this assumption is not made. However this document
deals exclusively with the conventional setting where $P(T < \infty) = 1$ is
assumed.

## Censoring

Data in a survival analysis are often subject to "censoring", which means that
we only have partial information about the value of $T$ for many subjects.
Suppose we are analyzing data from a medical study where we are studying the
incidence of stroke following implantation of a cardiac stent. We use age in
years as our time scale and let $T$ denote the age when a subject first has a
stroke.

Typically, only a subset of the subjects will have a stroke during our study.
Other subjects will be followed for a period of time and will never be
observed to have a stroke. Let $R$ denote the last age at which the person is
observed. If $T<R$, we observe $T$ but if $T > R$ we do not know the value of $T$. More formally, we
observe the time $Y={\rm min}(T, R)$ and the _status indicator_
$\delta = {\cal I}(Y = T)$. This is called
[right censoring](<https://en.wikipedia.org/wiki/Censoring_(statistics)>) --
we know that the value of $T$ is greater than some known value, but we do not
know the exact value of $T$.  In this context, $R$ is known as the _right censoring time_.

Right censoring is the most commonly-encountered form of censoring, but in
some settings we may have *left censoring*, meaning that we only know that $T$
is less than some observed value. Also, there is *interval censoring* in which
we know, for example, that someone had a stroke between the age of 75 and 77
but we do not know the exact age at which the stroke occurred.

A common assumption in survival analysis is that of *independent censoring*.
We will define this here in the context of right censoring. As above, let $T$
be the event time and $R$ be the right censoring time. We never observe both
$T$ and $R$. Nevertheless, we
can imagine that both values exist, one being a "latent" value. Independent
censoring simply means that $T$ and $R$ are independent random variables.
Concretely, this means, for example, that people who are more prone to having
an early event (e.g. unhealthy people) do not have systematically different
censoring times than people are are prone to having late events.

Since we don't observe $T$ and $R$ together, independent censoring is usually
an untestable assumption, but in some cases, based on the study design or other
external information, there may be reason to accept it
and in other cases there may be good reason to doubt that independent
censoring holds. There are various methods for effectively handling with
dependently censored data, but that is an advanced topic that we will not
consider further here.

## Truncation

An important concept in survival analysis is the potential selection bias
induced by *truncation* or *delayed entry*. If units (e.g. people) are
selected into the sample conditionally on their event time, then this must be
taken into account. The most common form of truncation is *left truncation* in
which there is a value $E$ (which may be specific to each observation) such
that if the event occurs before time $E$ (i.e. if $T < E$) then the person
would not have been included in our sample.

For example, in our analysis of stroke following stent placement, people
"enter" the study once their stent is placed. Let $E$ denote the age at which
a stent is placed. If a person never has a stent or if they have a stroke
before the stent is placed, this person cannot belong to our sample.

More formally, if we have left truncation then we are working with the
conditional distribution $P(T | T\ge E)$, while if there is no left truncation
then we are working with the unconditional distribution $P(T)$. If we are
doing regression analysis with covariates $X$, then with left truncation we
are studying $P(T | T\ge E, X)$ and with no such truncation we are studying
$P(T | X)$.

## Competing risks

In a survival analysis there may be other events that "compete" with the event
of interest. For example, if we are studying the time $T$ at which a person
has a stroke, it is possible that the person dies of another cause before
having a stroke. Death unrelated to stroke is a *competing risk* for the event
of interest.

## The risk set

The *risk set* at a specific time $t$ is the subset of units (e.g. people) who
could possibly experience the event at time $t$. Anyone who has already had
the event before time $t$, has been right censored before time $t$, or has
experienced a competing risk before time $t$ is not in this risk set.

## Recurrent events

In conventional survival analysis, each subject experiences the event one time
(although we may not observe this occurrence in our study as there may be
censoring). However some events may be able to recur. For example, suppose
that we are studying arrests by police and the time $T$ is the duration until
a person is arrested. There can be subsequent arrests for the same subject and
it may be of interest to study the distribution of event times for all
"spells" between one arrest and the next arrest.

## Notation and terminology

Subjects who are right-censored have a maximal *follow-up time* which is the
greatest time at which they were observed and confirmed not to have yet
experienced the event of interest. In survival analysis with right censoring,
where $T_i$ is the (possibly unobserved) event time and $R_i$ is the time at
which we would no longer be able to observe the subject, we typically write
$Y_i = {\rm min}(T_i, R_i)$ as the "observed" time, which is either the
follow-up time for censored subjects or the event time for non-censored
subjects. Then, we define the *status indicator* $\delta_i$ such that
$\delta_i=1$ if the event is observed and $\delta_i=0$ if the event is not
observed. Note that when $\delta_i=1$, then $Y_i=T_i$, and if $\delta_i=0$
then $T_i>Y_i$.

A more general notation that is often encountered is that for each subject we
have an interval $[L_i, R_i)$ such that the event is known to occur within
this interval. For right censored subjects, $R_i=\infty$. For non-censored
subjects, $L_i = R_i$. An *interval censored* subject has
$0 < L_i < R_i < \infty$, and a left-censored subject has
$0 = L_i < R_i< \infty$.

## Parametric and non-parametric methods

As in other areas of statistics, survival analysis can be conducted using
parametric or non-parametric methods. Moreover, "semi-parametric" methods play
an important role in survival analysis. While parametric methods can be
useful, survival analysis tends to emphasize non-parametric and
semi-parametric methods over parametric methods.

The most elementary parametric distribution used for survival times is the
[exponential distribution](https://en.wikipedia.org/wiki/Exponential_distribution),
although this is usually too simplistic of a model to use in practice. The
most commonly-encountered parameterized distribution in survival analysis is
arguably the
[Weibull distribution](https://en.wikipedia.org/wiki/Weibull_distribution),
and Gamma and log-normal distributions are also encountered.

## Estimation of the survival function

The *survival function* of a random variable $T$ is defined to be
$S(t) \equiv P(T>t)$. It is closely related to the
[cumulative distribution function](https://en.wikipedia.org/wiki/Cumulative_distribution_function)
(CDF), defined to be $F(t) = P(T\le t)$ since $S(t) = 1 - F(t)$. In words, the
survival function at time $t$ is the probability that the event has not
occurred by time $t$. Another name for the survival function is the _complementary CDF_. We may
also refer to this as the *marginal survival function* to emphasize that it is
not conditioned on any covariates.

The empirical CDF (eCDF) is one of the fundamental objects in statistics.
Based on an independent and identically distributed (IID) sample from some
distribution, the eCDF is defined as $\hat{F}(t) = \\#\\{T_i \le t\\} / n$. If
there is no truncation and if the values of $T_i$ are all observed (i.e. there
is no censoring) then we can estimate the survival function as
$\hat{S}(t) = 1 - \hat{F}(t)$.

As noted above, in survival analysis we usually have censoring and/or
truncation. We will consider here only the important subcase where there is
right censoring and no truncation. In this setting there is a simple estimator
of the survival function $S(t)$ known as the *product limit* estimator or the
[Kaplan-Meier](https://en.wikipedia.org/wiki/Kaplan%E2%80%93Meier_estimator)
estimator.

Let $t_1 < t_2 < \cdots < t_m$ denote the distinct times at which events are
observed to occur, let $d_i$ denote the number of events that occur at time
$t_i$, and let $n_i$ denote the size of the risk set just before time $t_i$.
Roughly speaking, the probability of passing through time $t_i$ without
experiencing the event is $1 - d_i/n_i$. Further, the probability of making it
from time 0 to time $t$ without experiencing the event is

$$ \hat{S}(t) \equiv \prod_{i:t_i\le t}(1 - d_i/n_i).  $$

This is the product limit estimator of the survival function.

There are many methods for statistical inference relating to survival
functions. The [log rank test](https://en.wikipedia.org/wiki/Logrank_test) is
a formal
[hypothesis test](https://en.wikipedia.org/wiki/Statistical_hypothesis_test)
of the null hypothesis that two survival functions are equal, i.e. the null
hypothesis $S_0(t) \equiv S_1(t)$. It is also possible to put confidence
intervals around the estimated survival function $\hat{S}(t)$ to convey the
precision with which it is estimated.

## Hazard functions

The [hazard function](https://en.wikipedia.org/wiki/Failure_rate) is a way of
mathematically representing a probability distribution that is commonly used
in survival analysis. The hazard function is defined to be

$$ h(t) \equiv \lim_{\delta\downarrow 0}P(T\le t+\delta|T\ge t)/\delta
= \lim_{\delta\downarrow 0}(S(t) - S(t+\delta)) / (\delta S(t)).
$$

The hazard function can be interpreted as the "instantaneous event rate". It
has units of 1/time so is not dimensionless (the time units matter).

It takes some practice to understand how to interpret this limit. If the time
unit is "days" and the hazard is 0.001 at day 100, then this means that
approximately 0.1% of the subjects at risk on day 100 will experience the
event on that day. Note that this is an approximate statement since we are not
actually taking a limit here. This approximate statement is closer to being
true over time intervals where the hazard function is approximately constant.

Note that the hazard function is not a probability and can be greater than 1
(but it must be non-negative). It is also not a density, although if the
density exists it can be determined from the hazard function.

The hazard function is related to the survival function in a fairly simple
way. First, define the *cumulative hazard function* as

$$ H(t) = \int_0^t h(s)ds.  $$

Then the identity $S(t) = \exp(-H(t))$ holds. If $T$ has a density $f$, then
$f = F^\prime = -S^\prime$. Thus

$$
f(t) = h(t)\exp(-H(t)).
$$

Note that this also implies that when densities exist, $h(t) = f(t)/S(t)$,
giving another natural view of the hazard funtion.

In many applications, the hazard function has a natural interpretation and may
be easier to interpret than the survival function. A common consideration is
whether the hazard function is increasing, decreasing, approximately constant,
or has some other shape like a "U" ("bathtub") shape. In terms of parametric
distributions, the exponential distribution has a constant hazard function,
whereas the Weibull distribution can have either an increasing or decreasing
hazard function depending on its parameters.

One common application of survival analysis is in the setting of failures of
manufactured products, e.g. how likely is it that your car will break down at
a particular point in time, given that it is currently operational? If the
hazard function is constant, then the car is equally likely to break down on
every day that you own it. If the hazard function is increasing, then as the
car gets older it becomes more likely to break down. This could be due to the
parts of the car wearing out and failing with use (e.g. due to material
fatigue).

If the hazard function is decreasing then as the car gets older it becomes
less likely to break down. This could occur if most failures are due to
manufacturing flaws. A car consists of thousands of parts, and suppose that in
any given car, a few parts may be flawed at the point of manufacture, or are
installed incorrectly. These flaws may not be sufficient to cause the car to
break down immediately, hence the car is sold in apparent working order.
However these flawed parts may tend to fail very early. Many cars need to be
serviced early in their lifespans for this reason. However, once you have
driven your car say 5000 miles, all of these flaws (if any) have been resolved
and the risk of a break-down becomes lower with time.

With manufactured products, both of the above mechanisms are likely to exist,
so the hazard function may exhibit a "bathtub" shape. This means that the
hazard function is higher for small $t$ due to manufacturing errors, then the
hazard function is lower once the car has been driven long enough to identify
and resolve any such issues. But as the car gets older and parts fatigue, a
different type of failure becomes more likely and the hazard function begins
to increase again.

A similar phenomenon exists with human lifespans, whereby the hazard of dying
is greater for infants and very young children (up to around age 3-5) and then
becomes very low for several decades before beginning to increase again around
age 50.

## Hazard ratios and hazard proportionality

A *hazard ratio* is the ratio between two hazard values. For example, we may
have two groups of subjects (e.g. people exposed or not exposed to a risk
factor), with each group having a hazard function $h_k(t)$ for $k=0, 1$
corresponding to not exposed and exposed people, respectively. The hazard
ratio at time $t$ is $h_1(t)/h_0(t)$. This is a very useful measure of the
"risk" associated with an exposure. For example, if the hazard ratio is 2 then
(roughly speaking) exposed people have twice the risk of experiencing the
event as non-exposed people.

Under an assumption of *proportional hazards* all hazard ratios are constant,
meaning in the present example that $h_1 \propto h_0$. As we will see below,
many popular methods for survival analysis assume proportional hazards, but it
is important to note that this assumed proportionality may not hold in any
particular setting.

## Estimating marginal hazard functions

The cumulative hazard function is easier to estimate than the hazard function,
but is more difficult to interpret. The most basic non-parametric estimate of
the cumulative hazard function under right censoring with no truncation is the
[Nelson-Aalen](https://en.wikipedia.org/wiki/Nelson-Aalen_estimator) estimator

$$ \hat{H}(t) = \sum_{i: t_i \le t} d_i/n_i, $$

using the same notation as used above when presenting the product-limit
estimator of the survival function. Since the hazard function $h$ is the
derivative of the cumulative hazard function $H$, it is possible to estimate
$h$ by numerically differentiating an estimate of $H$.

## Hazard regression

*Survival regression* is any method that aims to model conditional
distributions $P(T|X)$, where $T$ is an event time variable possibly subject
to censoring and/or truncation, and $X$ is a vector of explanatory variables.
If $T$ is fully observed, specialized techniques for survival regression are not needed. For
example, we may regress $T$ or $\log(T)$ on $X$ using least squares or a
generalized linear model (GLM). This type of direct approach has been extended
to accommodate censoring and truncation, leading to the so-called *accelerated
failure time* (AFT) models that we will not discuss further here.

The most widely-used approach to survival regression is arguably the
semi-parametric *proportional hazards* (PH) model, often known as the "Cox
model". This model is explicitly expressed in terms of the hazard function:

$$ h(t|X=x) = \exp(\beta^\prime x) h_0(t).  $$

The *baseline hazard function* $h_0$ is unknown and arbitrary, i.e. it is not
assumed to follow any parametric family. This is therefore a *semi-parametric*
model since is has a finite-dimensional parameter of interest $\beta$ and an
infinite-dimensional nuisance parameter $h_0$. It turns out that it is
possible to estimate $\beta$ using a type of maximum-likelihood technique
without simultaneously estimating $h_0$. This makes the PH model feel in
practice more like a conventional parametric model estimated using maximum
likelihood. The cumulative baseline hazard function can be estimated in a
separate step if desired, using a modified version of the Nelson-Aalen
estimator discussed above.

When interpreting the results of a PH model, the fact that it is based on
proportionality of the hazard function is key. Thus, a given regression slope
$\beta_j$ is the *log hazard ratio* when comparing the hazard functions for
two individuals who differ by one unit on variable $X_j$, and have identical
values for all other variables. The estimated hazard ratio for the
$j^{\rm th}$ covariate is simply $\exp(\hat{\beta}_j)$. Since the PH model
assumes proportionality of the hazard functions, this hazard ratio does not
depend on $t$ (although the true hazard ratio may depend on $t$ if the PH
model is incorrectly-specified).

The PH model is essentially a single-index model fit with maximum likelihood
techniques. Thus, once the concept of the hazard function and proportionality
of hazard functions is understood, familiar strategies for regression modeling
can be employed. For example, we can include interactions, splines, covariate
transformations, and conduct step-wise model searches. There are also versions
of information-based criteria such as AIC and BIC for model selection with PH
models. The PH model can be extended to accommodate left truncation and
competing risks.

As indicated above, proportionality of the conditional hazard functions is a
critical assumption in the PH model. This is not always an easy assumption to
check, but there are some methods based on residuals that can be employed.

## Time-varying covariates

In the survival regression model discussed above, all variables are defined at
"baseline". That is, every covariate $X_j$ is known at time zero and its value
cannot change. There are various approaches to survival regression that can
accommodate *time-varying covariates*, e.g. if a subject's status changes in a
way that changes their risk for the event of interest. The Cox PH regression
model discussed above can be extended to accommodate time-varying covariates
but we do not discuss that further here.

## Pseudo-observations

A "pseudo-observation" is a synthetic datapoint that combines the observed
time $Y$ and censoring status $\delta$ into a single real number. The
resulting value can then be used in many forms of statistical analysis, for
example as an independent or dependent variable in a regression, or in a
multivariate analysis such as PCA. Importantly, as discussed in more detail
below, the pseudo-observations are approximately independent and their
standard deviation reflects the underlying uncertainty in estimating the
parameter of interest. Using pseudo-observations allows survival analysis
(both estimation and inference) to be conducted using general-purpose
statistical methods instead of requiring specialized methods.

Here we will discuss pseudo-observations for the survival probability (i.e.
the function $S(t)$ evaluated at a specific time $t$), but note that it is
also possible to construct pseudo-observations for other quantities such as
the mean restricted life or the cumulative hazard.

Pseudo-observations are closely related to the
[jackknife](https://en.wikipedia.org/wiki/Jackknife_resampling). To motivate
the technique, let $\bar{X} = (X_1 + \cdots + X_n)/n$ be the sample mean of
$n$ observations from a common distribution. Let $\bar{X} _ {-i}$ denote the
"deleted" version of this statistic (the sample mean with the $i^{\rm th}$
observation deleted). These statistics satisfy the identity
$X_i = n\bar{X} - (n-1)\bar{X} _ {-i}$. Now consider the more general setting
where we have a statistic $\hat{\theta} _ {n}$ based on the full sample of
size $n$, and then we compute this statistic while deleting observation $i$,
to yield $\hat{\theta}_{-i}$. The pseudo-observation is defined to be

$$u_i \equiv n\hat{\theta} _ n - (n-1)\hat{\theta} _ {-i},$$

and is interpreted as the contribution of $X_i$ to the statistic of interest,
$\hat{\theta}_n$.

It can be shown that the $u_i$ are approximately independent, that
$\rm Avg\\{u_i\\} \approx \hat{\theta}_n$, and that

$${\rm SD}[u_i]/\sqrt{n} \approx {\rm SE}[\hat{\theta}_n].$$

Thus, the pseudo-observations approximately convert estimation and inference
for $\theta$ into a linear inference problem, analogous to estimating the
population mean with the sample mean. This idea could be useful in many
settings, but is particularly useful in the setting of survival analysis since
we can now treat the pseudo-observations $u_i$ like any other collection of
independent quantitative measurements, and analyze them using a wide variety
of statistical methods that are not otherwise adapted to survival analysis.

The most common construction of pseudo-observations for survival analysis is
based on the Kaplan-Meier (product limit) estimate of the marginal survival
function. We can compute $\hat{S}(t)$ using all data, and then we can compute
$\hat{S} _ {-i}(t)$ by deleting observation $i$. There are fast approximations
for doing this on large samples without repeating the full calculation for
each $i$. The pseudo-observation is
$u_i(t) = n\hat{S}(t) - (n-1)\hat{S}_{-i}(t)$. These can be used, for example,
in a regression analysis, regressing $u_i$ on covariates $x_i$, since it can
be shown that $E[u|x]$ can be interpreted as the probability of surviving to
time $t$ when the covariates are equal to $x$. Pseudo-observation conditional
variances, ${\rm var}[u|x]$, are generally not constant in $x$ (i.e. there is
heteroscedasticity), so typically regressions involving pseudo-observations
are fit with robust regression techniques such as using the Huber-White type
of inference.

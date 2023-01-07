# Survival analysis

Survival analysis is a set of techniques for characterizing distributions
and bulding models.  Traditionally, survival analysis methods are used
with *time to event* data, such as assessing the distribution of *failure
times* (e.g. the duration of time from when a person is diagnosed with
a disease until they progress to a more advanced stage of the disease).
There are many other applications of survival analysis methods that
have nothing to do with "survival" or "failure", and these methods can
even be applied in settings where "time" does not play a central role
(although this is less common).

## Time origin

Consider a setting where we are observing a person over time, and
we are interested in the time $T$ at which some event of interest
occurs.  It is very important to define the *origin* from which time
is measured, i.e. what is the meaning of time zero?  In the case of
human-centered studies, we may let $T$ denote the age of the person when
the event occurs, in which case the time origin is the date of birth.
Alternatively, there may be an event that must occur before the event of
interest, and the time of this event may make a more sensible time origin.
For example, if $T$ corresponds to graduating from university, we may
choose the date of matriculation as the time origin, so that e.g. $T=4$
means that the person graduated four years after beginning their studies.

## Event time distributions

The theoretical basis of conventional survival analysis is that we are
studying the distribution of a random variable $T$, corresponding to
the time at which an event of interest occurs.  In conventional survival
analysis, the event will eventually happen to everyone (i.e. $P(T<\infty)
= 1$).  We note that in some situations this may not be completely
realistic and there is a subdomain of survival analysis ("cure modeling")
in which this assumption is not made.  However this document deals
exclusively with the conventional setting where $P(T<\infty) = 1$
is assumed.

## Truncation

An important concept in survival analysis is the potential selection
bias induced by *truncation*.  If units (e.g. people) are selected into
the sample conditionally on their event time, then this must be taken
into account.  The most common form of truncation is *left truncation*
in which there is a value $Q$ (which may be specific to each person)
such that if the event occurs before time $Q$ (i.e. if $T < Q$) then
the person would not have been included in our sample.

More formally, if we have left truncation then we are working with
the conditional distribution $P(T | T\ge Q)$, while if there is no left
truncation then we are working with the unconditional distribution $P(T)$.
If we are doing regression analysis with covariates $X$, then with left
truncation we are studying $P(T | T\ge Q, X)$ and with no such truncation
we are studying $P(T | X)$.

An example of left truncation would be a medical study where we are
studying the incidence of stroke.  Suppose we use age as our time scale
(so the time origin is a person's date of birth), and we recruit older
people into our study based on some eligibility criteria.  Usually we
would not recruit people into our study if they have already had a stroke.
Suppose we recruit someone into the study who is currently 65 years
of age.  We are selecting this subject conditionally on the event $T\ge 65$,
and $Q=65$ is the subject's left truncation time.

## Censoring

Data in a survival analysis are often subject to "censoring", which
means that we only have partial information about the value of $T$ for
many subjects.  Continuing with the example where stroke is the event
of interest, only a subset of the subjects will have a stroke during
our study.  Suppose a person is recruited into the study at age 65 (their
left truncation time), and presently they are 75 years old and have not
yet had a stroke. In the conventional survival analysis this imaginary subject will have
a stroke at some age $T<\infty$, but since that event occurs in the
future we only know that $T > 75$.  This is called *right censoring* --
we know that the value of $T$ is greater than some known value, but we
do not know the exact value of $T$.

Right censoring is the most commonly-encountered form of censoring,
but in some settings we may have *left censoring*, meaning that we
only know that $T$ is less than some observed value.  Also, there is
*interval censoring* in which we know, for example, that someone had a
stroke between the age of 75 and 77 but we do not know the exact age
at which the stroke occurred.

A common assumption in survival analysis is that of *independent
censoring*.  We will define this here in the context of right censoring.
Let $T$ be the event time and $R$ be the right censoring time.
Note that only one of $T, R$ is observed, but we can
imagine that the other value exists as a "latent" value.
Independent censoring simply means that $T$ and $R$ are independent
random variables. Concretely, this means, for example, that people who
are more prone to having an early event (e.g. unhealthy people) do not
have systematically different censoring times than people are are prone
to having late events.

Since we don't observe $T$ and $R$ together, independent censoring is
usually an untestable assumption, but in some cases there may be reason
to accept it and in other cases there may be good reason to doubt that
independent censoring holds.  There are various methods for dealing with
dependent censored data, but that is an advanced topic that we
will not consider further here.

## Competing risks

In a survival analysis there may be other events that "compete" with
the event of interest. For example, if we are studying the time $T$
at which a person has a stroke, it is possible that the person dies of
another cause before having a stroke.  Death unrelated to stroke
is a *competing risk* for the event of interest.

## The risk set

The *risk set* at a specific time $t$ is the subset of units
(e.g. people) who could possibly experience the event at time $t$.
Anyone who has already had the event before time $t$, has been right
censored before time $t$, or has experienced a competing risk before
time $t$ is not in this risk set.

## Recurrent events

In conventional survival analysis, each subject experiences the event
one time (although we may not observe this occurrence in our study as
there may be censoring).  However some events may be able to recurr.
For example, suppose that we are studying arrests by police and the time
$T$ is the duration until a person is arrested.  There can be subsequent
arrests for the same subject and it may be of interest to study the
distribution of event times for all "spells" between one arrest and the
next arrest.

## Notation and terminology

Subjects who are right-censored have a maximal *follow-up time* which
is the greatest time at which they were observed and confirmed not to
have yet experienced the event of interest. In survival analysis with
right censoring, where $T_i$ is the (possibly unobserved) event time, we
typically write $Y_i$ as the observed time which is either the follow-up
time for censored subjects or the event time for non-censored subjects.
Then, we define the *status indicator* $\delta_i$ such that $\delta_i=1$
if the event is observed and $\delta_i=0$ if the event is not observed.
Note that when $\delta_i=1$, then $Y_i=T_i$, and if $\delta_i=0$ then
$T_i>Y_i$.

A more general notation that is often encountered is that for each
subject we have an interval $[L_i, R_i)$ such that the event is known to
occur within this interval.  For right censored subjects, $R_i=\infty$.
For non-censored subjects, $L_i = R_i$.  An *interval censored* subject
has $0 < L_i < R_i < \infty$, and a left-censored subject has
$0 = L_i < R_i< \infty$.

## Parametric and non-parametric methods

As in other areas of statistics, survival analysis can be conducted using
parametric or non-parametric methods.  Moreover, "semi-parametric" methods
play an important role in survival analysis.  While parametric methods
can be useful, survival analysis tends to emphasize non-parametric and
semi-parametric methods over parametric methods.

The most elementary parametric distribution in survival analysis is the
exponential distribution, although this is usually much too simplistic of
a model to use in practice. The most commonly-encountered parameterized
distribution in survival analysis is arguably the Weibull distribution,
and Gamma and log-normal distributions are also encountered.

## Estimation of the survival function

The *survival function* of a random variable $T$ is defined to be $S(t)
\equiv P(T>t)$.  It is closely related to the cumulative distribution
function (CDF) $F(t) = P(T\le t)$ since $S(t) = 1 - F(t)$.  In words,
the survival function at time $t$ is the probability that the event has
not occured by time $t$.

The empirical CDF (eCDF) is one of the fundamental objects in statistics.
Based on an independent and identically distributed (IID) sample from
some distribution, the eCDF is defined as $F(t) = \\#\\{T_i \le t\\} / n$.
If there is no truncation and if the values of $T_i$ are all observed
(i.e. there is no censoring) then we can estimate the survival function
as $\hat{S}(t) = 1 - \hat{F}(t)$.

As noted above, in survival analysis we usually have censoring and/or
truncation.  We will consider here only the important subcase where there
is right censoring and no truncation.  In this case there is a simple
estimator of the survival function $S(t)$ known as the *product limit*
estimator or the *Kaplan-Meier* estimator.  It is worth understanding
how this estimator is constructed.

Let $t_1 < t_2 < \cdots < t_m$ denote the distinct times at which events
are observed to occur, let $d_i$ denote the number of events that occur
at time $t_i$, and let $n_i$ denote the size of the risk set just before
time $t_i$.  Roughly speaking, the probability of passing through time
$t_i$ without experiencing the event is $1 - d_i/n_i$.  Further, the
probability of making it from time 0 to time $t$ without experiencing
the event is

$$ \hat{S}(t) \equiv \prod_{i:t_i\le t}(1 - d_i/n_i).  $$

This is the product limit estimator of the survival function.

There are many methods for statistical inference relating to survival
functions.  The *log rank test* is a formal hypothesis test of the
null hypothesis that two survival functions are equal, e.g. the null
hypothesis $S_0(t) \equiv S_1(t)$.  It is also possible to put confidence
intervals around the estimated survival function $\hat{S}(t)$ to convey
the precision with which it is estimated.

## Hazards

The *hazard* function is a way of mathematically representing a
probability distribution that is commonly used in survival analysis.
The hazard function is defined to be

$$ h(t) \equiv \lim_{\delta\downarrow 0}P(T\le t+\delta|T\ge t)/\delta.
$$

The hazard function can be interpreted as the "instantaneous event
rate".  It takes some practice to understand how to interpret this
limit. If the time unit is "days" and the hazard is 0.001 at day 100,
then this means that approximately 0.1%% of the subjects at risk on
day 100 will experience the event on that day.  Note that this is an
approximate statement since we are not actually taking a limit here.
This approximate statement is closer to being true over time intervals
where the hazard function is approximately constant.

Note that the hazard function is not a probability and can be greater
than 1 (but it must be non-negative).

The hazard function is related to the survival function in a fairly
simple way.  First, define the *cumulative hazard function* as

$$ H(t) = \int_0^t h(s)ds.  $$

Then the identity $S(t) = \exp(-H(t))$ holds.

The cumulative survival function is easier to estimate than the
hazard function, but is more difficult to interpret.  The most basic
non-parametric estimate of the cumulative hazard function under right
censoring with no trunction is the *Nelson-Aalen* estimator:

$$ \hat{H}(t) = \sum_{i: t_i \le t} d_i/n_i, $$

using the same notation as used above when presenting the product-limit
estimator of the survival function.  Since the hazard function $h$ is
the derivative of the cumulative hazard function $H$, it is possible to
estimate $h$ by numerically differentiating an estimate of $H$.

A *hazard ratio* is the ratio between two hazard values.  For example,
we may have two groups of subjects (e.g. people exposed or not exposed
to a risk factor), with each group having a hazard function $h_k(t)$ for
$k=0, 1$ corresponding to not exposed and exposed people, respectively.
The hazard ratio at time $t$ is $h_1(t)/h_0(t)$.  This is a very useful
and widely-utilized measure of the "risk" associated with an exposure.
For example, if the hazard ratio is 2 then roughly speaking exposed people
have twice the risk of experiencing the event as non-exposed people.

Under an assumption of *proportional hazards* all hazard ratios are
constant, meaning in the present example that $h_1 \propto h_0$. As
we will see below, many popular methods for survival analysis
assume proportional hazards, but it is important to note that this
proportionality may not hold in a particular setting.

## Hazard regression

*Survival regresion* is any method that aims to model conditional
distributions $P(T|X)$, where $T$ is an event time variable possibly
subject to censoring and/or truncation, and $X$ is a vector of explanatory
variables.

If $T$ is fully observed, specialized techniques may not be needed.
For example, we may regress $T$ or $\log(T)$ on $X$ using least squares
or a generalized linear model (GLM).  This type of direct approach
has been extended to accommodate censoring and truncation, leading to
the so-called *accelerated failure time* (AFT) models that we will not
discuss further here.

The most widely-used approach to survival regression is arguably the
semi-parametric *proportional hazards* (PH) model, often known as the
"Cox model".  This model is explicitly expressed in terms of the hazard
function:

$$ h(t|X=x) = \exp(\beta^\prime x) h_0(t).  $$

The *baseline hazard function* $h_0$ is unknown and arbitrary, i.e. it
is not assumed to follow any parametric family.  This is therefore a
*semi-parametric* model since is has a finite-dimensional parameter of
interest $\beta$ and an infinite-dimensional nuiscance parameter $h_0$.
It turns out that it is possible to estimate $\beta$ using a type of
maximum-likelihood technique without simultaneously estimating $h_0$.
This makes the PH model feel in practice like a conventional parametric
model estimated using maximum likelihood.  The cumulative baseline hazard
function can be estimated in a separate step if desired, using a modified
version of the Nelson-Aalen estimator discussed above.

When interpreting the results of a PH model, the fact that it is based on
proportionalithy of the hazard function is key.  Thus, a given regression
slope $\beta_j$ is the *log hazard ratio* when comparing the hazard
functions for two individuals who differ by one unit on variable $X_j$,
and have identical values for all other variables.  The estimated hazard
ratio for the $j^{\rm th}$ covariate is simply $\exp(\hat{\beta}_j)$.
Since the PH model assumes proportionality of the hazard functions, this
hazard ratio does not depend on $t$ (although the true hazard ratio may
depend on $t$ if the PH model is incorrectly-specified).

The PH model is essentially a single-index model fit with maximum
likelihood techniques.  Thus, once the concept of the hazard function and
proportionality of hazard functions is understood, familiar strategies
for regression modeling can be employed.  For example, we can include
interactions, splines, covariate transformations, conduct step-wise
model searches.  There are also versions of information-based criteria
such as AIC and BIC for model selection with PH models. The PH model
can be extended to accommodate left truncation and competing risks.

As indicated above, proportionality of the conditional hazard functions
is a critical assumption in the PH model.  This is not always an easy
assumption to check, but there are some residual-based methods that can
be employed.

## Time-varying covariates

In the survival regression model discussed above, all variables are
defined at "baseline".  That is, every covariate $X_j$ is known at
time zero and its value cannot change.  There are various approaches
to survival regression thqat can accommodate *time-varying covariates*,
e.g. if a subject experiences some other event that changes their risk
for the event of interest.  The Cox PH regression model discussed above
can be extended to accommodate time-varying covariates but we do not
discuss that further here.

#import "@preview/ilm:2.1.1": *

#show link: set text(fill: blue)

#set text(lang: "en")

#show: ilm.with(
  title: [Survival analysis and state transition models],
  authors: "Kerby Shedden",
  figure-index: (enabled: true),
  table-index: (enabled: true),
  listing-index: (enabled: true),
  chapter-pagebreak: false
)

A broad definition of _survival analysis_ would include any research setting in which we are monitoring units who can transition between discrete states over time.  The research aim of the survival analysis could be to gain an understanding any aspect of these transitions, such as the rate at which transitions occur, the dwell time in each state, statistical dependence between transitions, and how any of these phenomena are related to covariates.

Most commonly, survival analysis methods are used with _time to event_ or _duration_ data, where "time to event" refers to the duration of time from an origin until some event of interest occurs, such as as state transition. In the most basic example, there are usually only two states, typically described as "alive" and "dead".  All subjects begin in the "alive" state and eventually transition to the "dead" state, which is "absorbing".

A typical example with more than two states would be one in which subjects can be "healthy", "ill", or "dead", with all transitions allowed except that the "dead" state is absorbing -- once reached, there is no possibility to transition to another state.

An important aspect of survival analysis is that there is often only partial information about each subject's status, usually due to incomplete monitoring time.  For example, in a health study we would rarely be able to wait until all subjects have died, so some subjects will end the study in the "alive" state.  This is an example of "censoring", which will be defined more formally below.

= Key concepts

== Time origin

Consider a setting where we are monitoring a person over time, and are interested in the time $T$ at which some event of interest occurs. To begin, it is important to explicitly define the time
#link("https://en.wikipedia.org/wiki/Unit_of_measurement")[units], e.g. days, months, or years, and the _origin_ from which time is measured, i.e. what is the meaning of time $T=0$? For example, if $T$ denotes the age of a person when an event of interest occurs, the time origin is the date of birth. Alternatively, there may be an event that must occur before the event of interest, and the time of this event may make a more sensible time origin. For example, if $T$ corresponds to graduating from university, we may choose the date of university matriculation (first enrollment) as the time origin, so that e.g. if $T=4$ and the time units are years, then the person graduated four years after beginning their studies.

== Event time distributions

The theoretical basis of conventional survival analysis is that we are studying the probability distribution of a
#link("https://en.wikipedia.org/wiki/Random_variable")[random variable] $T$, corresponding to the time at which an event of interest occurs. In conventional survival analysis, the event will always occur if we wait long enough, so $P(T < infinity) = 1$. In some situations this may not be realistic and there is a subdomain of survival analysis called _cure modeling_ in which $P(T < infinity) < 1$ is allowed. However here we deal exclusively with the conventional setting where $P(T < infinity) = 1$ is assumed.

== Censoring

Data in a survival analysis are often subject to "censoring", which means that we only have partial information about the value of $T$ for many subjects. Suppose we are analyzing data from a medical study where we are studying the incidence of a complication following a medical procedure.  For example, if the population consists of people who have had a cardiac stent implanted, then a common complication would be stent thrombosis (a clot forming near the stent). In this case, we may use age in years as our time scale and let $T$ denote the age when a subject first has stent thrombosis.

Typically, only a subset of the subjects will have stent thrombosis during our study. Other subjects will be followed for a period of time and will never be observed to have stent thrombosis. Let $R$ denote the last age at which the person is observed. If $T < R$, we observe $T$ but if $T > R$ we do not know the value of $T$. More formally, we observe the time $Y = "min"(T, R)$ and the _status indicator_ $delta = "I"(Y = T)$. This is called #link("https://en.wikipedia.org/wiki/Censoring_(statistics)")[right censoring] -- we know that the value of $T$ is greater than some known value, but we do not know the exact value of $T$ for many subjects. In this context, $R$ is known as the _right censoring time_.

Right censoring is the most commonly encountered form of censoring, but in some settings we may have _left censoring_, meaning that we only know that $T$ is less than some observed value. Also, there is _interval censoring_ in which we know, for example, that someone had stent thrombosis between the age of 75 and 77 but we do not know the exact age at which the stent thrombosis occurred.

A common assumption in survival analysis is that of _independent censoring_. We will define this here in the context of right censoring. As above, let $T$ be the event time and $R$ be the right censoring time. We never observe both $T$ and $R$. Nevertheless, we can imagine that both values exist, one being a "latent" value. Independent censoring simply means that $T$ and $R$ are independent random variables (in a regression analysis, we may require the weaker condition that $T$ and $R$ are independent given the covariates). This independence implies, for example, that people who are prone to having an early event (e.g. unhealthy people) do not have systematically different censoring times than people who typically have late events.

Since we don't observe $T$ and $R$ together, independent censoring is almost always an untestable assumption. In some cases, based on the study design or other external information, there may be reason to accept it and in other cases there may be good reason to doubt that independent censoring holds. There are various methods for effectively handling dependently censored data, but that is an advanced topic that we will not consider further here.

== Truncation

An important concept in survival analysis is the potential selection bias induced by _truncation_ or _delayed entry_. If units (e.g. people) are selected into the sample conditionally on their event time, then this must be taken into account. The most common form of truncation is _left truncation_ in which there is a value $L$ (which may be specific to each observation) such that if the event occurs before time $L$ (i.e. if $T < L$) then the person would not have been included in our sample.

For example, suppose in our analysis of stent thrombosis following stent placement that we only have access to data from one health insurance provider. People may join this health insurance program at any age, and we can let $L$ denote the age at enrollment. If a person had a stent placed before joining this insurance company's plan, we would not know the exact age when the stent was placed. Therefore, we may choose to eliminate from the analysis all people who already had a stent when joining the health insurance plan. Thus, not having a stent at the beginning of the insurance record is a requirement for being in our study, so the age at which a person began their health insurance coverage would be a left truncation time.

More formally, if we have left truncation then we are working with the conditional distribution $P(T | T >= L)$, while if there is no left truncation then we are working with the unconditional distribution $P(T)$. If we are doing regression analysis with covariates $X$, then with left truncation we are studying $P(T | T >= L, X)$ and with no such truncation we are studying $P(T | X)$.

== Competing risks

In a survival analysis there may be other events that "compete" with the event of interest. For example, if we are studying the time $T$ at which a person has a stroke, it is possible that the person dies of another cause before having a stroke. Death unrelated to stroke is a _competing risk_ for the event of interest.  This can be viewed as another example of a multi-state transition model, with states "healthy", "had stroke", and "dead".  Each person has a stroke age $T$, a death age $D$, and a censoring age $R$.  Conceptually, there is a joint distribution for $(R, D, T)$, but if $D < T$ we will not know the value of $T$.  People can transition from "healthy" to "had stroke", from "had stroke" to "dead", and from "healthy" to "dead".

Superficially, $D$ appears to be another form of censoring, so we may consider redefining $R$ as $R' = "min"(R, D)$.  But doing this may create dependent censoring. Even if $R$ is independent of $T$, $R'$ may not be, for example if less healthy subjects are at greater risk for both stroke and death.

== The risk set

The _risk set_ at a specific time $t$ is the subset of units (e.g. people) who could possibly experience the event at time $t$. Anyone who has already had the event before time $t$, has been right censored before time $t$, has experienced a competing risk before time $t$, or is truncated after time $t$ is not in the risk set at time $t$.

== Recurrent events

In conventional survival analysis, each subject experiences the event exactly once (although we may not observe this occurrence in our study as there may be censoring). However some events may be able to recur. For example, suppose that we are studying arrests by police and the time $T$ is the duration until a person is arrested. There can be subsequent arrests for the same individual and it may be of interest to study the distribution of event times for all "spells" between one arrest and the next arrest.

This setting can be viewed as another example of a multi-state model.  In the example, we would have a state for "never arrested", "arrested 1 time", "arrested 2 times", etc.  Subjects can transition from each state to the next one.

== Counting process notation

Subjects who are right-censored have a maximal _follow-up time_ which is the greatest time at which they were observed and confirmed not to have yet experienced the event of interest. In survival analysis with right censoring, where $T_i$ is the (possibly unobserved) event time and $R_i$ is the time at which we would no longer be able to observe the subject, we typically write $Y_i = "min"(T_i, R_i)$ as the "observed" time, which is either the follow-up time for censored subjects or the event time for non-censored subjects. Then, we define the _status indicator_ $delta_i$ such that $delta_i=1$ if the event is observed and $delta_i=0$ if the event is not observed. Note that when $delta_i=1$, then $Y_i=T_i$, and if $delta_i=0$ then $T_i > Y_i$.

A more general notation that is often encountered is that for each subject we have an interval $[L_i, R_i)$ such that the event is known to occur within this interval. For right censored subjects, $R_i = infinity$. For non-censored subjects, $L_i = R_i$. An _interval censored_ subject has $0 < L_i < R_i < infinity$, and a left-censored subject has $0 = L_i < R_i< infinity$.

For more general recurrent event data or when there are more than two states, we can decompose each individual's history into a series of disjoint records $L_j$, $Y_j$, $delta_j$, where $L_j$ denotes the time at the beginning of the record, $Y_j$ denotes the time at the end of the record, and $delta_j$ denotes the state at the end of the record.  Here, instead of $delta_j$ being a binary variable, it is an indicator of the state, with one possible state being "censored".

#link("https://arxiv.org/pdf/2210.07114")[Here] is a thorough treatment of the use of counting process notation in survival analysis.

= Parametric and non-parametric methods

As in other areas of statistics, survival analysis can be conducted using parametric or non-parametric methods. Moreover, "semi-parametric" methods play an important role in survival analysis. While parametric methods can be useful, survival analysis tends to emphasize non-parametric and semi-parametric methods over parametric methods.

The most elementary parametric distribution used for survival times is the #link("https://en.wikipedia.org/wiki/Exponential_distribution")[exponential distribution], although this is usually too simplistic of a model to use in practice. The most commonly-encountered parameterized distribution in survival analysis is arguably the
#link("https://en.wikipedia.org/wiki/Weibull_distribution")[Weibull distribution], and Gamma and log-normal distributions are also encountered.

= Estimation of the survival function

The _survival function_ of a random variable $T$ is defined to be
$S(t) equiv P(T>t)$. It is closely related to the
#link("https://en.wikipedia.org/wiki/Cumulative_distribution_function")[cumulative distribution function] (CDF), defined to be $F(t) = P(T<= t)$ since $S(t) = 1 - F(t)$. In words, the survival function at time $t$ is the probability that the event has not occurred by time $t$. Another name for the survival function is the _complementary CDF_. We may also refer to this as the _marginal survival function_ to emphasize that it is not conditioned on any covariates.

The empirical CDF (eCDF) is one of the fundamental objects in statistics. Based on an independent and identically distributed (IID) sample from some distribution, the eCDF is defined as $hat(F)(t) = frac(\#{T_i <= t}, n, style: "horizontal")$. If there is no truncation and if the values of $T_i$ are all observed (i.e. there is no censoring) then we can estimate the survival function as $hat(S)(t) = 1 - hat(F)(t)$.

As noted above, in survival analysis we usually have censoring and/or truncation. We will consider here only the important subcase where there is right censoring and no truncation. In this setting there is a simple estimator of the survival function $S(t)$ known as the _product limit_ estimator or the
#link("https://en.wikipedia.org/wiki/Kaplan%E2%80%93Meier_estimator")[Kaplan-Meier] estimator.

The Kaplan-Meier estimator focuses exclusively on the observed event times. Let $t_1 < t_2 < dots.c < t_m$ denote the distinct times at which events are observed to occur, let $d_i$ denote the number of events that occur at time $t_i$, and let $n_i$ denote the size of the risk set just before time $t_i$. The estimated probability of passing through time $t_i$ without experiencing the event is $1 - d_i/n_i$. Thus, the estimated probability of making it from time 0 to time $t$ without experiencing the event is

$
hat(S)(t) equiv product_(i:t_i<= t)(1 - d_i/n_i).
$

This is the product limit estimator of the survival function. Note that the actual survival function can be any non-increasing right continous function, and thus in general $S(t)$ will change at infinitely many values of $t$. However the Kaplan-Meier estimate of the survival function is a step function that only changes at the observed values of $t$ where an event occurs (just as the eCDF only changes at the observed data values).  A basic exercise is to prove that the product limit estimate $hat(S)$ is the complement of the eCDF, i.e. $hat(S) = 1 - hat(F)$.

There are many methods for statistical inference relating to survival functions. The #link("https://en.wikipedia.org/wiki/Logrank_test")[log rank test] is a formal
#link("https://en.wikipedia.org/wiki/Statistical_hypothesis_test")[hypothesis test] of the null hypothesis that two survival functions are equal, i.e. the null hypothesis $S_0(t) equiv S_1(t)$. It is also possible to put both pointwise and simultaneous confidence intervals (bands) around the estimated survival function $hat(S)(t)$ to convey the precision with which it is estimated.

= Hazard functions

The #link("https://en.wikipedia.org/wiki/Failure_rate")[hazard function] is a way of mathematically representing a probability distribution that is commonly used in survival analysis. The hazard function is defined to be

$
h(t) equiv lim_(delta arrow.b 0)P(T <= t + delta | T >= t) / delta = lim_(delta arrow.b 0)(S(t) - S(t+delta)) / (delta S(t)).
$

The hazard function can be interpreted as the "instantaneous event rate". It has units of 1/time so is not dimensionless (the time units matter).

It takes some practice to understand how to interpret this limit. If the time unit is "days" and the hazard is 0.001 at day 100, then this means that approximately 0.1% of the subjects at risk on day 100 will experience the event on that day. Note that this is an approximate statement since we are not actually taking a limit here. This approximate statement is closer to being true over time intervals where the hazard function is approximately constant.

Note that the hazard function is not a probability and can be greater than 1 (but it must be non-negative). It is also not a density, although if the density exists it can be determined from the hazard function.

When the survival function is smooth, the hazard function is also the logarithmic derivative of the survival function:

$
h(t) = -"d"/"dt" log S(t).
$

Another important quantity to understand is the _cumulative hazard function_

$
H(t) = integral_0^t h(s)"ds".
$

The identity $S(t) = exp(-H(t))$ holds. If $T$ has a density $f$, then $f = F' = -S'$. Thus

$
f(t) = h(t) exp(-H(t)).
$

Note that this also implies that when densities exist, $h(t) = frac(f(t), S(t), style: "horizontal")$, giving another natural view of the hazard function.

In many applications, the hazard function may be easier to interpret than the survival function. A common consideration is whether the hazard function is increasing, decreasing, approximately constant, or has some other shape like a "U" ("bathtub") shape. In terms of parametric distributions, the exponential distribution has a constant hazard function, whereas the Weibull distribution can have either an increasing or decreasing hazard function depending on its parameters.

One common application of survival analysis is in the setting of failures of manufactured products, e.g. how likely is it that your car will break down at a particular point in time, given that it is currently operational? If the hazard function is constant, then the car is equally likely to break down on every day that you own it. If the hazard function is increasing, then as the car gets older it becomes more likely to break down. This could be due to the parts of the car wearing out and failing with use (e.g. due to material fatigue).

If the hazard function is decreasing then as the car gets older it becomes less likely to break down. This could occur if most failures are due to manufacturing flaws. A car consists of thousands of parts, and suppose that in any given car, a few parts may be flawed at the point of manufacture, or are installed incorrectly during manufacture. These flaws may not be sufficient to cause the car to break down immediately, hence the car is sold in apparent working order. However these flawed parts may tend to fail very early. Many cars need to be serviced early in their lifespans for this reason. However, once you have driven your car say 5000 miles, these flaws have been discovered and resolved and the risk of a break-down becomes lower with time (while the other type of failure due to wear out and fatigue has not yet become common).

With manufactured products, both of the above mechanisms are likely to exist, so the hazard function may exhibit a "bathtub" shape. This means that the hazard function is higher for small $t$ due to manufacturing errors, then the hazard function is lower once the car has been driven long enough to identify and resolve any such issues. But as the car gets older and parts fatigue, a different type of failure becomes more likely and the hazard function begins to increase again.

A similar phenomenon exists with human lifespans, whereby the hazard of dying is greater for infants and very young children (up to around age 3-5) and then becomes very low for several decades before beginning to increase again around age 40.  Prior to the advent of modern medicine, most deaths were either before age 5 (infant mortality) or after around age 40 (when diseases of aging begin to appear).

== Hazard ratios and hazard proportionality

A _hazard ratio_ is the ratio between two hazard values. For example, we may have two groups of subjects (e.g. people exposed or not exposed to a risk factor), with each group having a hazard function $h_k (t)$ where $k=0, 1$ corresponds to not exposed and exposed people, respectively. The hazard ratio at time $t$ is $frac(h_1(t), h_0(t), style: "horizontal")$. This is a very useful measure of the "risk" associated with an exposure. For example, if the hazard ratio is 2 then (roughly speaking) exposed people have twice the risk of experiencing the event as non-exposed people.

Under an assumption of _proportional hazards_ all hazard ratios are constant, meaning in the present example that $h_1 prop h_0$. As we will see below, many popular methods for survival analysis assume proportional hazards, but it is important to note that this assumed proportionality may not always hold in practice.

== Estimating marginal hazard functions

The cumulative hazard function is easier to estimate than the hazard function, but is more difficult to interpret. The most basic non-parametric estimate of the cumulative hazard function under right censoring with no truncation is the #link("https://en.wikipedia.org/wiki/Nelson-Aalen_estimator")[Nelson-Aalen] estimator

$
hat(H)(t) = sum_(i: t_i <= t) frac(d_i, n_i, style: "horizontal")
$

using the notation introduced above. Since the hazard function $h$ is the derivative of the cumulative hazard function $H$, it is possible to estimate $h$ by numerically differentiating a smooth estimate of $H$.

== Propotional hazards regression

_Survival regression_ is any method that aims to model conditional distributions $P(T|X)$, where $T$ is an event time variable possibly subject to censoring and/or truncation, and $X$ is a vector of explanatory variables. If $T$ is fully observed, specialized techniques for survival regression are not needed. For example, we may regress $T$ or $log(T)$ on $X$ using least squares or a generalized linear model (GLM). This type of direct approach has been extended to accommodate censoring and truncation, leading to the so-called _transformation models_ and _accelerated failure time_ (AFT) models that we will not discuss further here.

The most widely-used approach to survival regression is arguably the semi-parametric _proportional hazards_ (PH) model, often known as the "Cox model". This model is explicitly expressed in terms of the hazard function:

$
h(t|X=x) = exp(beta^prime x) h_0(t).
$

The _baseline hazard function_ $h_0$ is unknown and arbitrary, i.e. it is not assumed to follow any parametric family. This is therefore a _semi-parametric_ model since is has a finite-dimensional parameter of interest $beta$ and an infinite-dimensional nuisance parameter $h_0$. It turns out that it is possible to estimate $beta$ using a type of maximum-likelihood technique without simultaneously estimating $h_0$. This makes the PH model feel in practice more like a conventional parametric model estimated using maximum likelihood. The cumulative baseline hazard function can be estimated in a separate step if desired, using a modified version of the Nelson-Aalen estimator discussed above.

When interpreting the results of a PH model, remembering that it is based on proportionality of the hazard function is key. Thus, a given regression slope $beta_j$ is the _log hazard ratio_ that compares the hazard functions for two individuals who differ by one unit on variable $X_j$, and have identical values for all other variables. The estimated hazard ratio for the $j^"th"$ covariate is simply $exp(hat(beta)_j)$. Since the PH model assumes proportionality of the hazard functions, this hazard ratio does not depend on $t$ (although the true hazard ratio may depend on $t$ if the PH model is incorrect).

The PH model is essentially a single-index model fit with maximum likelihood techniques. Thus, once the concept of the hazard function and proportionality of hazard functions is understood, familiar strategies for regression modeling can be employed. For example, we can include interactions, splines, covariate transformations, and conduct step-wise model searches. There are also versions of information-based criteria such as AIC and BIC for model selection with PH models. The PH model can be extended to accommodate left truncation and competing risks.

As indicated above, proportionality of the conditional hazard functions is a critical assumption in the PH model. This is not always an easy assumption to check, but there are some methods based on residuals that can be employed.  Independent censoring is also an important assumption for the PH model to be meaningful.

== Additive hazards regression

Although much less common than the proportional hazards approach, another approach to survival regression based on hazard functions exists which is useful in some particular situations.  In the _Aalen additive hazard regression model_, the population hazard function has the form

$
h(t | x) = h_0(t) + sum_(j=1)^p alpha_j (t) x_j (t),
$

where the $x_j (t)$ are time-varying covariates, and the $alpha_j (t)$ are time-varying coefficients.  The presence of both forms of time-varying structure makes this model very flexible.  The cumulative coefficients $beta_j (t) = integral_(s=0)^t alpha_j (s)"ds"$ play an important role in the estimation and interpretation for this model.  If a covariate function $x_j (dot.c)$ increases by one unit from time $0$ to time $t$, then the cumulative hazard of the event being modeled increases by $beta_j(t)$.

Notably, this model can be fit using ordinary least squares.  For each time $t$ that is an observed event time ($t in {t_1, t_2, ..., t_m}$), let $X^t$ denote the design matrix based on the covariates at time $t$, including only subjects who are at-risk at time $t$.  Similarly, let $y^t$ denote a vector indicating all subjects who have the event at time $t$, also restricted to subjects at risk.  Let $hat(alpha)_t$ denote the ordinary least squares slopes for the linear regression of $y^t$ on $X^t$.  Then define $hat(beta)_t = sum_(s <= t) hat(alpha)_s$.  These $hat(beta)_t$ are estimates of the cumulative coefficient functions $beta_j (t)$.

There is very little information available to estimate the $alpha_j (t)$ parameters -- if there are no ties, then only one failure occurs at each time $t_j$.  Therefore, the $hat(alpha)_j$ do not concentrate well around the $alpha_j (t)$.  However the cumulative coefficients behave much better -- the $hat(beta)_j (t)$ do concentrate around the $beta_j (t)$.  The $hat(beta)_j$ are asymptotically normal and confidence bands can be calculated for them, although the procedure for doing so is somewhat complex.  Further, it is possible to gain a small amount of additional power by using weighted least squares, taking account of the heteroscedasticity inherent in the fact that $y^t$ has binary elements.  A semiparametric version of the method models some of the coefficient functions as constant in time.

In practice, the additive hazard model is most often used when there are one or more time-varying covariates with frequently changing values, and/or when the hazards are unlikely to be proportional, or when they change over time.  Also, the additive hazards model is sometimes used without inference to provide a quick visual check of the shape of the hazard functions and to assess covariate effects.

== Cause-specific hazard regression

If there is a competing risk, it may be meaningful to fit a _cause-specific hazard regression_, treating subjects as being censored when they experience the competing risk.  As noted above when discussing the Kaplan-Meier estimator, it is important to be careful when treating a competing risk as a form of censoring.  When estimating the marginal survival function, it is rarely advisable to do this, due to the risk of introducing dependent censoring.  However in the case of survival regression (proportional or additive), it may be plausible that we have independent censoring after conditioning on covariates.

= Time-varying covariates

In the survival regression model discussed above, all variables are defined at "baseline". That is, every covariate $X_j$ is known at time zero and its value cannot change. There are various approaches to survival regression that can accommodate _time-varying covariates_, e.g. if a subject's status changes in a way that changes their risk for the event of interest. The Cox PH regression model discussed above can be extended to accommodate time-varying covariates, but we do not discuss that further here.

= Cumulative incidence functions

Suppose we are in the competing risks setting, where we have a primary endpoint of interest that occurs at time $T$, and another "competing" event that occurs at time $D$.  A person is censored at time $R$.  The _cumulative incidence_ of the primary endpoint is a function $"CI"(t)$, telling us the probability of having the primary event before time $t$.  It is defined as

$
"CI"(t) = integral_0^t h_0(t)S_a (t)"dt",
$

where $S_a(t)$ is the _all cause survival function_ and $h_0$ is the _cause-specific hazard function_.  The all cause survival function is the probability of not having experienced either the primary or the competing event by time $t$.  It can be estimated by defining $Y = "min"(T, D)$ and using the product limit estimator of the survival function.  The cause-specific hazard function is the hazard for having the primary event at time $t$, given that you are at risk at time $t$.  It can be estimated using the Nelson-Aalen estimator of the hazard function, treating death as a form of censoring.

Cumulative incidence is a convenient way to handle competing risks and avoids the challenges of working with so-called _sub-distribution hazards_, which many people find confusing.  It is also possible to introduce covariates into the analysis, to assess how subjects with different characteristics have different cumulative incidences.

Cumulative incidence is most useful when you are interested in the burden (e.g. social or economic) of a condition on individuals or on society.  For example, if the primary event is dementia and the competing risk is death, we may want to know the cumulative incidence (probability) of getting dementia by, say, age 85.  This probability includes people who get dementia and subsequently die before age 85 (i.e. it is not the proportion of living 85 year old people wiith dementia).  On the other hand, if someone dies young, say at age 60, then they are "saved" from getting dementia by their death.  As a result, high risk subpopulations like smokers may appear to have a lower cumulative incidence of dementia, even if smoking increases both the risk of dementia and the risk of death.  In this setting, the cumulative incidence function would be telling us that smokers are less likely to experience dementia in their lifetimes than non-smokers.

= Pseudo-observations

A _pseudo-observation_ (or _pseudo-value_) is a synthetic datapoint that combines the observed time $Y$ and censoring status $delta$ into a single real number. The resulting value can then be used in many forms of statistical analysis, for example as an independent or dependent variable in a regression, or in a multivariate analysis such as PCA. Importantly, as discussed in more detail below, the pseudo-observations are approximately independent and their standard deviation reflects the underlying uncertainty in estimating the parameter of interest. Using pseudo-observations allows survival analysis (both estimation and inference) to be conducted using general-purpose statistical methods instead of requiring specialized methods.

Here we will discuss pseudo-observations for the survival probability (i.e. the function $S(t)$ evaluated at a specific time $t$), but note that it is also possible to construct pseudo-observations for other quantities such as the mean restricted life or the cumulative hazard.

Pseudo-observations are closely related to the #link("https://en.wikipedia.org/wiki/Jackknife_resampling")[jackknife]. To motivate the technique, let $overline(X) = frac((X_1 + ... + X_n), n, style: "horizontal")$ be the sample mean of $n$ observations from a common distribution. Let $overline(X)_(-i)$ denote the "deleted" version of this statistic (the sample mean with the $i^"th"$ observation deleted). These statistics satisfy the identity $X_i = n overline(X) - (n-1)overline(X)_(-i)$. Now consider the more general setting where we have a statistic $hat(theta)_n$ based on the full sample of size $n$, and then we compute this statistic while deleting observation $i$, to yield $hat(theta)_(-i)$. The pseudo-observation is defined to be

$
u_i equiv n hat(theta) _ n - (n-1)hat(theta)_(-i),
$

and is interpreted as the contribution of $X_i$ to the statistic of interest, $hat(theta)_n$.

It can be shown that the $u_i$ are approximately independent, that

$
"Avg"{u_i} approx hat(theta)_n,
$

and that

$
"SD"(u_i) / sqrt(n) approx "SE"[hat(theta)_n].
$

Thus, the pseudo-observations approximately convert estimation and inference for $theta$ into a linear inference problem, analogous to estimating the population mean with the sample mean. This idea could be useful in many settings, but is particularly useful in the setting of survival analysis since we can now treat the pseudo-observations $u_i$ like any other collection of independent quantitative measurements, and analyze them using a wide variety of statistical methods that are not otherwise adapted to survival analysis.

The most common construction of pseudo-observations for survival analysis is based on the Kaplan-Meier (product limit) estimate of the marginal survival function. We can compute $hat(S)(t)$ using all data, and then we can compute $hat(S)_(-i)(t)$ by deleting observation $i$. There are fast approximations for doing this on large samples without repeating the full calculation for each $i$. The pseudo-observation is $u_i (t) = n hat(S)(t) - (n-1)hat(S)_(-i)(t)$. These can be used, for example, in a regression analysis, regressing $u_i$ on covariates $x_i$, since it can be shown that $E[u|x]$ can be interpreted as the probability of surviving to time $t$ when the covariates are equal to $x$.

There are a few limitations of the pseudo-observation approach, but for each there are work-arounds.  Pseudo-observation conditional variances, $"var"[u|x]$, are generally not constant in $x$, i.e. there is heteroscedasticity. Therefore, typically regressions involving pseudo-observations are fit with robust regression techniques such as using the Huber-White type of inference.  A second issue is that pseudo-observation regressions require unconditional independent censoring ($T_i$ is independent of $R_i$ without conditioning on the covariate $X_i$).  This limitation can be resolved by using inverse probability weighting to account for censoring that depends on measured covariates.

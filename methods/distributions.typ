#import "@preview/ilm:2.1.1": *

#show link: set text(fill: blue)

#set text(lang: "en")

#show: ilm.with(
  title: [Characterizing Distributions],
  authors: "Kerby Shedden",
  figure-index: (enabled: true),
  table-index: (enabled: true),
  listing-index: (enabled: true),
  chapter-pagebreak: false
)

= Introduction

Probability distributions are the central object of interest in probability theory and statistics. The most general treatment of probability distributions involves #link("https://en.wikipedia.org/wiki/measure_(mathematics)")[measure theory]. In applied work, we usually are able to use more elementary representations of probability distributions, including the #link("https://en.wikipedia.org/wiki/Probability_density_function")[probability density function] (pdf), the #link("https://en.wikipedia.org/wiki/Cumulative_distribution_function")[cumulative distribution function] (cdf), the _quantile function_, and the #link("https://en.wikipedia.org/wiki/Moment-generating_function")[moment generating function] (mgf). The space of probability distributions is infinite-dimensional, so in practice we cannot estimate every aspect of a distribution with a finite amount of data, and instead work with finite dimensional numerical summaries such as the mean, variance, and quantiles.

The field of statistics (as opposed to the field of probability) focuses on using data to estimate either the full probability distribution or a summary quantity describing a certain aspect of a probability distribution. The full distribution for a numerical random variable or vector can be estimated using the #link("https://en.wikipedia.org/wiki/Empirical_distribution_function")[empirical cdf] (an estimator of the cdf), or the #link("https://en.wikipedia.org/wiki/Histogram")[histogram] (an estimator of the pdf). Summary measures often have a natural estimator, such as the sample mean (an estimator of the population mean), and the sample median (an estimator of the population median).

The discussion here focuses on univariate distributions of a quantitative random variable. The settings of categorical (non-quantitative) and multivariate data involve different methods. The two most common characteristics used to summarize probability distributions on the real line are #link("https://en.wikipedia.org/wiki/Moment_(mathematics)")[moments] and #link("https://en.wikipedia.org/wiki/Quantile")[quantiles]. Both of these approaches can provide us with measures of _location_ (also known as _centrality_ or _central tendency_), and measures of #link("https://en.wikipedia.org/wiki/Statistical_dispersion")[dispersion] (also known as _scale_). For example, a moment-based measure of location is the mean, while a quantile-based measure of location is the median. A moment-based measure of dispersion is the
#link("https://en.wikipedia.org/wiki/Standard_deviation")[standard deviation], while quantile-based measures of dispersion include the
#link("https://en.wikipedia.org/wiki/Interquartile_range")[inter-quartile range] (IQR) and the #link("https://en.wikipedia.org/wiki/Median_absolute_deviation")[median absolute deviation]. "Higher order" summary characteristics of a distribution such as
#link("https://en.wikipedia.org/wiki/Skewness")[skewness] and #link("https://en.wikipedia.org/wiki/Kurtosis")[kurtosis] are less commonly encountered but are of great interest in certain settings.

Below we summarize some less familiar characteristics of probability
distributions, and discuss how to estimate these characteristics from data.

= Extremes, heavy-tailed distributions, and tail parameter estimation

Many important research questions involve the frequency of "extreme" events, for example major earthquakes, large movements in financial markets, extremely long human lifespans, or extreme rainfall events. The study of extremes naturally leads us to focus on the right tail of a probability distribution. In some cases the extremes of interest may lie in the left tail, but in that case we can flip the distribution (multiply the values by -1 so the left and right tails are swapped). Therefore by convention methods for studying extremes focus on the right tail.

In the statistical study of extremes, we do not attempt to classify individual data points as being "extreme" or "non-extreme". In some cases there may be non-statistical reasons to define a threshold beyond which an observation is extreme (for example, we may describe a hurricane as extreme if the wind speed exceeds 160 miles per hour). However, there is no objective statistical basis for defining a single observation as being extreme. Instead, we usually study extremes by characterizing the tail of the probability distribution according to its asymptotic rate of decay.  We will also consider an alternative approach to thinking about extremes based on a measure of kurtosis.

Recall that the cumulative distribution function (cdf) of a random variable $X$ is the function $F(t) = P(X <= t)$, viewed as a function of $t in RR$. The _complementary cdf_ (ccdf), also known as the
_survival function_, is the right tail probability $S(t) = P(X > t) = 1 - F(t)$. To understand the frequency of extreme (large) values, we can consider how rapidly the tail probability converges to zero as
$t$ increases. In many familiar distributions, the tails are _exponential_ meaning that

$
P(X > t) = L(t) dot exp(-frac(t, mu, style: "horizontal")),
$

where $L(t)$ is a #link("https://en.wikipedia.org/wiki/Slowly_varying_function")[slowly varying function] and $mu$ is a scale parameter. If $L(t)$ is constant, we have the exponential distribution, but the property of having exponentially decaying tails is much more general, and includes, for example, all of the gamma distributions. The normal distribution has tails that are thinner than exponential ("light tailed") since $P(X > t) = L(t) dot exp(-frac(t^2, mu, style: "horizontal"))$ for appropriate choices of $L(t)$ and $mu$.

In a #link("https://en.wikipedia.org/wiki/Heavy-tailed_distribution")[heavy tailed distribution], the tail probabilities shrink more slowly than an exponential rate, which means that for all $k > 0$,

$
lim_(t arrow.r.long infinity) exp(k dot t) dot P(X > t) = infinity.
$

The prototypical heavy-tailed distribution has a #link("https://en.wikipedia.org/wiki/Power_law")[power law] tail, meaning that

$
P(X > t) = frac(L(t), t^alpha, style: "horizontal"),
$

where $alpha$ is the _tail index_ (you may also encounter the _shape parameter_ defined as $xi = frac(1, alpha, style: "horizontal")$). In such distributions, the $k^"th"$ moment only exists and is finite if $k < alpha$ (the greater the value of $alpha$, the more moments exist).  The non-existence of moments turns out to drive many of the consequences of a distribution having a heavy tail, so it is natural to view $alpha$ as a measure of how heavy the tail is.

The simplest family of distributions with power law tails is the
#link("https://en.wikipedia.org/wiki/Pareto_distribution")[Pareto distribution]. The sample space of the Pareto distribution is $[1, infinity)$, and the CCDF is $P(X > t) = frac(1, t^alpha, style: "horizontal")$. This distribution has a tail index of $alpha$ as defined above. If $U$ follows a uniform distribution, then $U^(-frac(1, alpha, style: "horizontal"))$ is Pareto. Alternatively, if $Y$ follows a standard exponential distribution, then $Z = exp(frac(Y, alpha, style: "horizontal"))$ follows a Pareto distribution.

= Exceedances

The Pareto and exponential distributions are one-parameter families and will not fit many datasets well. Furthermore, when focusing on extreme values we usually don't want to become distracted by the structure of the center of the distribution. Therefore, we need a more flexible way to model the tail of a probability distribution.

One way to focus on the tail is to convert the data to _exceedances_. This means that we select a parameter $T$ and replace the dataset ${X_i}$ with the dataset ${X_i-T | X_i >= T}$.

If $T$ is appropriately selected then the exceedances may follow a Pareto or exponential distribution, even though these models are a poor fit to the full dataset. However we will want to use a more flexible two-parameter family of models in most cases.

= Tail plots

Before considering formal estimation and inference for the tail of a
distribution, we will discuss some graphical approaches that capture the structure of the tail of a distribution. These approaches consider the upper _order statistics_ of a sample of data and plot them in log space to best reflect the shape of the tail. Recall that the $j^"th"$ order statistic is the $j^"th"$ sorted value in our data, sorted in increasing order.

Let $X_((j))$ denote the $j^"th"$ order statistic either of our data, or of the exceedances derived from our data. This order statistic corresponds to a _probability point_ or "plotting position". A plotting position that takes the probability represented by the $j^"th"$ order statistic to be

$
frac(j-a, n+1-2a, style: "vertical")
$

for a parameter $0 < a < 1$. For plotting positions with $a=0$,

$
P(X <= X_((j))) approx frac(j, (n+1), style: "horizontal").
$

Suppose that the tail of $X$ is a power law with tail index $alpha$. Then we have

$
P(X > X_((j))) &=
frac(c, X_((j))^alpha, style: "horizontal")\
& approx 1 - frac(j, (n+1), style: "horizontal").
$

Therefore in log-space we have

$
log(1 - frac(j, (n+1), style: "horizontal")) approx log(c) - alpha log(X_((j))).
$

This implies that if we plot the probability points $1 - frac(j, (n+1), style: "horizontal")$ against the order statistics $X_((j))$ in log-space, we obtain an approximate linear relationship with slope $-alpha$. This type of probability plot can be used both as a means for estimating $alpha$, and as a diagnostic for whether the tails actually follow a power law.

Alternatively, if our tails are exponential we have the relationship

$
P(X > X_((j))) & = c dot exp(-frac(X_((j)), mu, style: "horizontal")) \
& approx 1 - frac(j, (n+1), style: "horizontal"),
$

therefore

$
log(1 - frac(j, (n+1), style: "horizontal")) approx log(c) - frac(X_((j)), mu, style: "horizontal").
$

In a _semi-log_ plot (log transforming the probability points but not the order statistics), when the distribution has exponential tails we obtain a linear relationship with slope $-frac(1, mu, style: "horizontal")$.

Using such probability tail plots, we can estimate distributional parameters ($alpha$ or $mu$) by fitting a least squares regression line to the points in the plot. The number of points used in the least squares fit is a tuning parameter that must be selected, typically in the range $20$ to $200$. These estimators are convenient, intuitive, and _distributionally robust_ (since they depend on the assumed form of the tail but do not require a complete specification of a probability model). However these estimators may not be very efficient (i.e. they may have high estimation variance). Alternative estimators that may be more efficient are discussed below.

= The Hill estimate of the tail parameter

If a distribution has a power-law tail $P(X > t) = frac(c, t^alpha, style: "horizontal")$, we can solve for the upper quantiles, yielding the quantile function

$
Q(p) = (c/(1-p))^(1/alpha).
$

Since the order statistics estimate quantiles, we have

$
X_((j)) approx (c / (1 - frac(j, (n+1), style: "horizontal")))^(1/alpha).
$

An estimator known as the #link("https://en.wikipedia.org/wiki/Heavy-tailed_distribution#Hill's_tail-index_estimator")[Hill estimator]
begins by considering the ratios of upper order statistics

$
X_((n-j)) / X_((n-k)) approx ((j+1)/(k+1))^(1/alpha).
$

Thus

$
log X_((n-j)) / X_((n-k)) approx -alpha^(-1) log((j+1)/(k+1)).
$

If we hold $k$ fixed and average these log ratios for $1 <= j < k$, we get

$
hat(A) equiv (k-1)^(-1) sum_(j=1)^(k-1) log X_((n-j)) / X_((n-k)) approx -alpha^(-1) sum_(j=1)^(k-1) log((j+1)/(k+1)).
$

This establishes a relationship between the statistic $hat(A)$ and the quantity of interest $alpha$. The constant of proportionality turns out to be nearly equal to $1$:

$
sum_(j=1)^(k-1) log((j+1)/(k+1)) approx integral_0^k log((x+1)/(k+1))"dx" arrow.r.long -1.
$

Therefore, the _Hill estimate_ of the tail parameter is

$
hat(alpha)_"Hill" = frac(1, hat(A), style: "horizontal").
$

In the Hill estimate, the value of $k$ is a tuning parameter. To select $k$, we usually calculate $hat(alpha)_"Hill"$ for various values of $k$ (typically $k approx 20-200$) and choose a value that corresponds to a stable range of values of the estimate.

If the data are exactly Pareto, the maximum likelihood estimate (MLE) of $alpha$ is the Hill estimate using $k=n$.

= The generalized Pareto distribution

As noted above, the one-parameter Pareto distribution may not fit many data sets well, and further has an awkward sample space of $[1, infinity)$. To address these issues, the #link("https://en.wikipedia.org/wiki/Generalized_Pareto_distribution")[generalized Pareto distribution] was developed, which has sample space $[0, infinity)$ and complementary CDF

$
P(X > t) = sigma^(-1)(1 + t/(sigma alpha))^(-alpha).
$

Note that the ratio of the CCDF's for the Pareto and generalized Pareto satisfies

$
frac(sigma^(-1)(1 + frac(t, (sigma alpha), style: "horizontal"))^(-alpha), t^(-alpha)) arrow.r.long sigma^(-1)(sigma alpha)^alpha != 0,
$

so the generalized Pareto and Pareto distributions both have power law tails with the same tail parameter $alpha$.

As $alpha arrow.r.long infinity$ (or equivalently $xi = frac(1, alpha, style: "horizontal") arrow.r.long 0$), the generalized Pareto distribution becomes the exponential distribution.

The famous #link("https://en.wikipedia.org/wiki/Pickands%E2%80%93Balkema%E2%80%93De_Haan_theorem")[Pickands-Balkema-De Haan theorem] demonstrates that with appropriate choice of threshold $T$, the exceedances
for many distributions approximately follow a generalized Pareto distribution. This theorem plays the role of the central limit theorem in the study of extremes, since it allows us to use a specific parametric model to study data that may follow a large range of distributions.

= Block maxima and the generalized extreme value distribution

Another way to approach extremes is to partition the data into blocks, and calculate the maximum observed value within each block. According to the #link("https://en.wikipedia.org/wiki/Fisher%E2%80%93Tippett%E2%80%93Gnedenko_theorem")[Fisher-Tippett-Gnedenko theorem],
the distribution of these values should be well approximated by a #link("https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution")[generalized extreme value distribution] (GEV), which is a three-parameter distribution. This is an example of a central limit theorem-like result for extremes, since a wide variety of populations have block maxima that are well-approximated by the GEV distribution.

The block maxima approach is often used with serially observed data (time series), and the block is a coarse resolution of time. For example, if our time series consists of daily values, we might choose a block size of one year. To estimate a GEV, we need to have a sufficient number of blocks. If there are too few blocks, say 50 or fewer, then the GEV parameter estimates will be very uncertain.

One advantage of working with block-wise maxima is that they are less sensitive to _positive serial dependence_ that causes clusters of extreme values to occur in close proximity to each other. These clusters will generally occur within one block and only the largest among them will influence the results. On the other hand, fitting a generalized Pareto distribution to the exceedances may produce biased estimates due to such dependence.

= Likelihood-based estimation

Likelihood-based estimation is generally more statistically efficient than moment or quantile-matching estimates such as considered above. The most well-known likelihood-based estimator is the maximum likelihood estimator (MLE), which is asymptotically fully efficient under standard conditions. However, the MLE can be non-unique and difficult to compute. Additionally, in the case of the GEV the support of the distribution depends on the parameters, violating the standard conditions on which theoretical guarantees about the MLE are based.

For the standard one-parameter Pareto distribution, the MLE of the tail index is simply

$
hat(alpha)_"MLE" = 1/("Avg"(log(X_i))).
$

and the MLE of the $mu$ in the exponential distribution is simply
$hat(mu) = overline(X_i)$. These MLE's are not very useful in practice since few populations can be well-approximated by these one-parameter distributions.

For families with more than a single parameter, the MLE is usually computed numerically, and results may depend strongly on starting values. The generalized Pareto distribution has two parameters and the generalized extreme value distribution has three parameters. Both of these distributions are challenging to work with numerically. Alternative likelihood-based estimators have been developed including the empirical Bayes estimator discussed #link("https://www.jstor.org/stable/40586625")[here].

Maximum likelihood estimates for the generalized extreme value (GEV)
distribution can be calculated numerically, but good starting values are needed to obtain robust convergence. One way to get good starting values is using the _probability weighted moments_ approach discussed
#link("https://www.stat.cmu.edu/technometrics/80-89/VOL-27-03/v2703251.pdf")[here].

= Return levels

The _m-observation return level_ (or simply $m$-return) is a value $T$ that is expected to be exceeded once out of every $m$ observations. In other words, $T$ is the m-observation return level if $I_i = "I"(X_i > T)$ and $E[I_1 + ... + I_m] = 1$, where the $X_i$ are identically distributed random variables. Since $E[I_1 + ... + I_m] = m E[I]$, where $I$ has the same distribution as the $I_i$, the $m$-return can be inferred from the equation $m E[I] = 1$ or $E[I] = frac(1, m, style: "horizontal")$. Since $E[I] = P(X > t) = 1 - P(X <= t)$, if $F(x)$ is the cumulative distribution function (CDF) of a random variable $X$, then the m-observation return level (for $m$ independent copies of $X$) is the solution to $F(x) = 1 - frac(1, m, style: "horizontal")$. Thus, the m-observation return level is the $1 - frac(1, m, style: "horizontal")$ quantile of $X$.

If we are working with exceedances, we need to take account of the observations that were excluded due to falling below the threshold $T$. Let $q$ denote the proportion of the full dataset that exceeds $T$. Then the m-observation return level is the $1 - frac(1, (q dot m), style: "horizontal")$ quantile calculated using the distribution of exceedances.

= Interpretation of moments

A classical statistical moment is defined to be the expected value of a random variable raised to a power. For example, the raw first and second moments are $E[X]$ (the expected value) and $E[X^2]$. In practice we usually work with the _centered moments_ (or _central moments_), for example $E[(X - E X)^2]$ is the centered second moment, which is better known as the variance.

In most situations, if you know all the moments of a distribution then you know everything about the distribution (there are some technical conditions for this claim to be true, as it is based on the invertibility of the moment generating function). But this fact is not very useful in practice because it is nearly impossible to estimate high order moments $E[(X - E X)^k]$ for large values of $k$. The sample estimator of this moment is $n^(-1) sum (X_i - overline(X))^k$, and this estimator is consistent and asymptotically unbiased, but if $k > 2$ it has huge mean squared error for any practically realistic sample size.

It is important to have fluent ways to discuss what we learn by studying moments in a data analysis. The _location_ captures the "central value" of a distribution, but exactly what this means depends on what measure of location is being used. Although the mean is so common that we take it for granted, its interpretation follows from a somewhat opaque physical analogy, the balancing point for mass distributed along a line.  The median can be interpreted as the "deepest" point in a distribution in that it is "surrounded by other values on all sides".  The mode is the most common value in the distribution.  The mean can fail to exist or be infinite, and the median and mode may not be uniquely defined.

_Dispersion_ captures how far observations tend to fall from the central value, and _skewness_ captures whether the most extreme values tend to fall more often on one side of the central value than the other. The interpretation of _kurtosis_ has been debated.  #link("https://www.tandfonline.com/doi/full/10.1080/00031305.2014.917055")[This] reference argues that kurtosis should be viewed solely as a measure of the heaviness of the tails.

A natural question is whether two distributions can always be compared in terms of measures of location, scale, etc.  Bickel and Lehman #link("https://projecteuclid.org/journals/annals-of-statistics/volume-4/issue-6/Descriptive-Statistics-for-Nonparametric-Models-III-Dispersion/10.1214/aos/1176343648.full")[answered] this question in the negative.  For example, they introduced the notion of _dispersive order_ to assess when probability distributions can be compared based on dispersion.  If $Q_1$ and $Q_2$ are the quantile functions of two distributions, then if

$
Q_1(p) - Q_1(q) < Q_2(p) - Q_2(q)
$

for all $q < p$, the distribution represented by $Q_2$ is more dispersed than the distribution represented by $Q_1$.

Later we will consider how these ideas can be generalized to multivariate data.  For example, in searching for a meaningful multivariate generalization of quantiles, the research area of #link("https://link.springer.com/chapter/10.1007/978-1-4613-0045-8_4")[data depth] emerged.

= Quantile analogues to moments

Moments and quantiles are fundamentally different -- at the sample level, moments involve averaging whereas quantiles involve sorting.  But there are several bridges between moments and quantiles, such as the identity $E[X] = integral_0^1 Q(p)"dp"$, where $Q(dot.c)$ is the quantile function.

Moments are often used to characterize properties of a distribution such as location, dispersion, and skewness.  There are familiar quantile-based approaches for achieving this same goal.  The median is a quantile-based measure of location and the inter-quartile range (IQR) is a quantile-based measure of dispersion.  But is there a general way to define quantile analogues to all possible moments?

Taking skewness as an example, measures of the form

$
(Q(p) + Q(1-p) - 2Q(frac(1, 2, style: "horizontal"))) / (Q(p) - Q(1-p)) =

((Q(p) - Q(frac(1, 2, style: "horizontal"))) - (Q(frac(1, 2, style: "horizontal")) - Q(1-p))) / (Q(p) - Q(1-p))
$

have been proposed, using, e.g. $p=frac(3, 4, style: "horizontal")$ or $p=frac(9, 10, style: "horizontal")$.

Kurtosis can be defined as a ratio of a measure of "tail dispersion" to a measure of "inner dispersion", such as

$
(Q(p) - Q(1 - p)) / (Q(q) - Q(1-q))
$

where $p > q > 1/2$.

These approaches are practically useful, but don't provide a fully satisfying solution since they include tuning parameters (e.g. $p, q$), and do not immediately provide constructions for moments with order greater than 4.

= L-moments

As noted earlier, many descriptive statistics are either moments or quantiles. If high order moments are hard to estimate, perhaps there is a quantile-based analogue to these moments that is easier to estimate? This idea led to the development of #link("https://en.wikipedia.org/wiki/L-moment")[L-moments] which are linear combinations of order statistics (order statistics in turn are
estimates of quantiles).

The definition of an L-moment of arbitrary order will be shown below, but first we present the first four L-moments as special cases.

The first L-moment $lambda_1$ is the same as the usual mean.

The second L-moment of a distribution represented through the random variable $X$ is defined as $lambda_2 = frac((E X_(2:2) - E X_(1:2)), 2, style: "horizontal")$. Here, $X_(j:k)$ is defined to be the random variable obtained by sampling $k$ independent values from the distribution of $X$, and then taking the $j^"th"$ largest among them.

The third L-moment is

$
lambda_3 &= frac((E X_(3:3) - 2 E X_(2:3) + E X_(1:3)), 3, style: "horizontal") \
&= frac(((E X_(3:3) - E X_(2:3)) - (E X_(2:3) - E X_(1:3))), 3, style: "horizontal").
$

From the second expression we can see that the third L-moment measures the asymmetry between the upper and lower order statistic differences.

Finally, the fourth L-moment is

$
lambda_4 &= frac((E X_(4:4) - 3E X_(3:4) + 3E X_(2:4) - E X_(1:4)),  4, style: "horizontal")\
&= frac(((E X_(4:4) - E X_(3:4)) + (E X_(2:4) - E X_(1:4)) - 2(E X_(3:4) - E X_(2:4))),  4, style: "horizontal").
$

From the second expression we can see that the fourth L-moment measures the asymmetry between the outer and inner order statistic differences.

Often we work with the standardized third and fourth L-moments, $lambda_3^s = frac(lambda_3, lambda_2, style: "horizontal")$ and $lambda_4^s = frac(lambda_4, lambda_2, style: "horizontal")$.

Note that these standardized L-moments are _scale invariant_ meaning that their value is not changed by scaling the data. All L-moments except for the first L-moment are _translation invariant_, meaning that their values are not changed by adding a constant to all data values. Scale and translation
invariance (also known as _affine invariance_) are important because they imply that the result does not depend on the units or origin of the measurement scale.

L-moments are useful descriptive statistics that capture the shape of distributions. They are more robust (less sensitive to contamination) than the classical moments, and one can estimate higher order L-moments than is practical with classical moments.

A connection has emerged between L-moments and the study of heavy tails.  The standardized fourth L-moment $lambda_4^s$ is a measure of tail heaviness. It is not directly equivalent to the tail index $alpha$, but arguably captures a similar characteristic of a distribution. It has been argued (see #link("https://www.tandfonline.com/doi/full/10.1080/00031305.2024.2402898#abstract")[here]) that if the standardized fourth L-moment exceeds 0.35, then the tails are sufficiently heavy that conventional statistical inference is "disrupted".

== Estimation of L-moments

Population L-moments are defined in terms of the values $E X_(j:k)$ as defined above.  An unbiased estimate of this quantity based on iid data $X_1, ..., X_n$ can be obtained by taking all subsets of size $k$ from the data, selecting the $j^"th"$ largest value from each subset, and averaging these values.  It turns out that this is a linear combination of the order statistics $X_((i))$ of the observed data.  Specifically, $X_((i))$ will be the $j^"th"$ largest of $k$ values exactly

$
binom(i-1, j-1) dot.c binom(n-i, i-j)
$

times. If we let

$
c_(i j k) = frac(binom(i-1, j-1) dot.c binom(n-i, i-j), binom(n, k), style: "horizontal"),
$

then

$
sum_i c_(i j k) X_((i))
$

is an unbiased estimator of $X_(j:k)$.

== A general definition for population values of L-moments

Above we gave definitions for the first four L-moments as special cases.  In general, an L-moment is a linear functional of the quantile function as given by

$
lambda_m = integral_0^1 Q(p) tilde(P)_(n-1)(p) "dp",
$

where $tilde(P)_n$ is the $n^"th"$ #link("https://en.wikipedia.org/wiki/Legendre_polynomials#Shifted_Legendre_polynomials")[shifted Legendre polynomial]

$
tilde(P)_n = (-1)^n sum_(k=0)^n binom(n, k) binom(n+k, k) (-x)^k.
$

These polynomials form an orthogonal basis on $[0, 1]$.  By studying the graphs of these polynomials, it becomes clear why they are capturing features of a probability distribution that can be interpreted as location, dispersion, skewness, and kurtosis.

== L-moment relationships

Suppose we have a population that can be meaningfully stratified into many subpopulations, such as county of residence for US adults.  We can then estimate summary statistics such as L-moments within each subpopulation.  In many cases two summary statistics, e.g. measuring location and scale, will be related in informative ways.

When working with classical moments, it is often noted that the mean and variance are related.  This is known as a _mean/variance relationship_.  One possible way this might arise if the distributions follow a family such as the Poisson family, where the variance is equal to the mean.  In other settings, we may find that the variance has a different fixed relationship to the mean, such as the variance being proportional to the mean (_quasi-Poisson_), proportional to the square of the mean (_quasi-gamma_), or a linear combination of the mean and its square (_quasi negative binomial_).

== L-comoments

Just as covariance is a bivariate analogue to the univariate variance, the univariate L-moments can be extended to define bivariate measures of association known as #link("https://www.sciencedirect.com/science/article/pii/S0047259X07000103")[L-comoments]. Given a bivariate random vector $(X, Y) in RR^2$, the goal of an L-comoment is to assess the extent to which the value of $X$ predicts where $Y$ falls relative to its marginal distribution $F_Y$.  Unlike conventional covariance and correlation, the L-comoments are not symmetric in $X$, $Y$.

Population L-comoments are defined as

$
"Cov"(X, tilde(P)_(n-1)(Y)),
$

where as above the $tilde(P)_n$ are shifted Legendre polynomials.

The L-covariance, much like the conventional covariance, tells us whether small values of $X$ tend to co-occur with values of $Y$ falling in one tail of $F_Y$, while large values of $X$ tend to co-occur with values of $Y$ falling in the other tail of $F_Y$.  The L-coskewness tells us whether large values of $X$ tend to co-occur with values of $Y$ falling in either tail of $F_Y$.  The L-cokurtosis tells us whether large values of $X$ tend to co-occur with values of $Y$ that can occur in the far portion of either tail of $F_Y$.

An alternative representation of population L-comoments is based on _concomitants_.  Let $X^((Y))_(j:k)$ denote the value obtained by sampling $k$ iid copies of $(X, Y)$, and taking the value $X_i$ such that $Y_i = Y_(j:k)$.  That is, we sort the pairs $(X_i, Y_i)$ with respect to $Y$ to obtain $Y_(j:k)$, and then take the value of $X$ that accompanies $Y_(j:k)$.  The $k^"th"$ L-comoment of $X$ with respect to $Y$ is a linear combination of expected concomitants

$
lambda^(X Y)_k = k^(-1) sum_(j=0)^(k-1) (-1)^j binom(k-1, j) E X^((Y))_(k-j:k).
$

Unbiased estimates of the L-comoments can be obtained by estimating the concomitants of $X$ with respect to $Y$ in the full sample of size $n$.  Define the coefficients

$
w_(k r n) = sum_(j=0)^("min"(r-1, k-1)) (-1)^(k-1-j) frac(binom(k-1, j) binom(k-1+j, j) binom(r-1, j), binom(n-1, j), style: "horizontal").
$

then

$
hat(lambda)^(X Y)_k = n^(-1) sum_(r=1)^n w_(k r n) X^((Y))_(r:n).
$

Comoments based on classical moments also exist, for example, the coskewnesses can be defined as $"Cov"(X, Y^2)$ and $"Cov"(X^2, Y)$. These suffer from some of the same practical challenges as other statistics based on higher order moments.


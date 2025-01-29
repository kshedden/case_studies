# Characterizing distributions

Probability distributions are the central object of interest in probability
theory and statistics. The most general treatment of probability distributions
involves
[measure theory](<https://en.wikipedia.org/wiki/measure_(mathematics)>). In
applied work, we usually are able to use more elementary representations of
probability distributions, including the
[probability density function](https://en.wikipedia.org/wiki/Probability_density_function)
(pdf), the
[cumulative distribution function](https://en.wikipedia.org/wiki/Cumulative_distribution_function)
(cdf), and the
[moment generating function](https://en.wikipedia.org/wiki/Moment-generating_function)
(mgf). The space of probability distributions is infinite-dimensional, so in
practice we often work with finite dimensional numerical summaries such as the
mean and variance.

The field of statistics (as opposed to the field of probability) focuses on
using data to estimate either the full probability distribution or a summary
quantity describing certain aspects of a probability distribution. The full
distribution can be estimated using the
[empirical cdf](https://en.wikipedia.org/wiki/Empirical_distribution_function)
(an estimator of the cdf), or the
[histogram](https://en.wikipedia.org/wiki/Histogram) (an estimator of the
pdf), among other tools not discussed here. Summary measures often have a
natural estimator, such as the sample mean (an estimator of the population
mean), and the sample median (an estimator of the population median).

The discussion here focuses on univariate distributions of a quantitative
random variable. The settings of categorical (non-quantitative) and
multivariate data involve different methods.

The two most common characteristics used to summarize univariate distributions
of a quantitative value (i.e. probability distributions on the real line) are
[moments](<https://en.wikipedia.org/wiki/Moment_(mathematics)>) and
[quantiles](https://en.wikipedia.org/wiki/Quantile). Both of these approaches
can provide us with measures of *location* (also known as *centrality* or
*central tendency*), and measures of
[dispersion](https://en.wikipedia.org/wiki/Statistical_dispersion) (also known
as *scale*). For example, a moment-based measure of location is the mean,
while a quantile-based measure of location is the median. A moment-based
measure of dispersion is the
[standard deviation](https://en.wikipedia.org/wiki/Standard_deviation), while
quantile-based measures of dispersion include the
[inter-quartile range](https://en.wikipedia.org/wiki/Interquartile_range)
(IQR) and the [MAD](https://en.wikipedia.org/wiki/Median_absolute_deviation)
(median absolute deviation).

"Higher order" summary characteristics of a distribution such as
[skew](https://en.wikipedia.org/wiki/Skewness) and
[kurtosis](https://en.wikipedia.org/wiki/Kurtosis) are less commonly
encountered but are of great interest in certain settings.

Below we summarize some less familiar characteristics of probability
distributions, and discuss how to estimate these characteristics from data.

## Extremes, heavy-tailed distributions, and tail parameter estimation

Many important research questions involve the frequency of "extreme" events,
for example major earthquakes, large movements in financial markets, extremely
long human lifespans, or extreme rainfall events. The study of extremes
naturally leads us to focus on the right tail of a probability distribution.
In some cases the extremes of interest may lie in the left tail, but in that
case we can flip the distribution (multiply the values by -1 so the left and
right tails are swapped). Therefore by convention methods for studying
extremes focus on the right tail.

In the statistical study of extremes, we do not attempt to classify individual
data points as being "extreme" or "non-extreme". In some cases there may be
non-statistical reasons to define a threshold beyond which an observation is
extreme (for example, we may describe a hurricane as extreme if the wind speed
exceeds 160 miles per hour). However, there is no objective statistical basis
for defining a single observation as being extreme. Instead, we study extremes
by characterizing the tail of the probability distribution according to its
asymptotic rate of decay.

Recall that the cumulative distribution funtion (cdf) of a random variable $X$
is the function $F(t) = P(X \le t)$, viewed as a function of
$t \in {\mathbb R}^+$. The *complementary cdf* (ccdf), also known as the
*survival function*, is the right tail probability
$S(t) = P(X > t) = 1 - F(t)$. To understand the frequency of extreme (large)
values, we should consider how rapidly the tail probability converges to zero
as $t$ increases. In many familiar distributions, the tails are *exponential*
meaning that $P(X > t) = L(t)\cdot \exp(-t/\mu)$, where $L(t)$ is a
[slowly varying function](https://en.wikipedia.org/wiki/Slowly_varying_function)
and $\mu$ is a scale parameter. For our purposes, we can treat $L(t)$ as a
constant, so having exponential tails implies that
$P(L > t) \propto \exp(-t/\mu)$.

In a
[heavy tailed distribution](https://en.wikipedia.org/wiki/Heavy-tailed_distribution),
the tail probabilities do not shrink exponentially fast, which means that for
all $k > 0$,

$$
\lim_{t\rightarrow \infty} \exp(k\cdot t) \cdot P(X > t) = \infty.
$$

Many heavy-tailed distributions have
[power law](https://en.wikipedia.org/wiki/Power_law) tail, meaning that
$P(X > t) = L(t) / t^\alpha$, where $\alpha$ is the *tail index* or *shape
parameter* (note that in some settings, the shape parameter is defined to be
$1/\alpha$). In such distributions, the $k^{\rm th}$ moment only exists and is
finite if $k < \alpha$ (the greater the value of $\alpha$, the more moments
exist).

The simplest family of distributions with power law tails is the
[Pareto distribution](https://en.wikipedia.org/wiki/Pareto_distribution). The
sample space of the Pareto distribution is $[1, \infty)$, and the CCDF is
$P(X > t) = 1/t^\alpha$. This distribution has a tail index of $\alpha$ as
defined above. If $U$ follows a uniform distribution, then $U^{-1/\alpha}$ is
Pareto. Alternativey, if $Y$ follows a standard exponential distribution, then
$Z = \exp(Y / \alpha)$ follows a Pareto distribution.

### Exceedances

The Pareto and exponential distributions are one-parameter families and will
not fit many datasets well. Furthermore, when focusing on extreme values we
usually don't want to become distracted by the structure of the center of the
distribution. Therefore, we need a more flexible way to model the tail of a
probability distribution.

One way to focus on the tail is to convert the data to *exceedances*. This
means that we find a parameter $T$ and replace the dataset $\\{X_i\\}$ with
the dataset $\\{X_i-T | X_i \ge T\\}$.

If $T$ is appropriately selected then the exceedances may follow a Pareto or
exponential distribution, even though these models are a poor fit to the full
dataset. However we will want to use a more flexible two-parameter family of
models in most cases.

### Tail plots

Before considering formal estimation and inference for the tail of a
distribution, we will discuss some graphical approaches that capture the
structure of the tail of a distribution. These approaches consider the upper
*order statistics* of a sample of data and plot them in log space to best
reflect the shape of the tail. Recall that the $j^{\rm th}$ order statistic is
the $j^{\rm th}$ sorted value in our data, sorted in increasing order.

Let $X_{(j)}$ denote the $j^{\rm th}$ order statistic either of our data, or
of the exceedances derived from our data. This order statistic corresponds to
*probability point* or "plotting position" $j/(n+1)$ (there is a more general
definition of plotting position that takes the probability to be
$(j-a)/(n+1-2a)$ for a parameter $0 < a < 1$). For the plotting position with
$a=0$, $P(X \le X_{(j)}) \approx j/(n+1)$ and
$P(X > X_{(j)}) \approx 1 - j/(n+1)$.

Suppose that the tail of $X$ is a power law with tail index $\alpha$. Then we
have

$$
P(X > X_{(j)}) = c/X_{(j)}^\alpha \approx 1 - j/(n+1).
$$

Therefore in log-space we have

$$
\log(1 - j/(n+1)) \approx \log(c) - \alpha\log(X_{(j)}).
$$

This implies that if we plot the probability points $1 - j/(n+1)$ against the
order statistics $X_{(j)}$ in log-space, we obtain an approximate linear
relationship with slope $-\alpha$. This type of probability plot can be used
both as a means for estimating $\alpha$, and as a diagnostic for whether the
tails actually follow a power law.

Alternatively, if our tails are exponential we have the relationship

$$
P(X > X_{(j)}) = c\cdot\exp(-X_{(j)}/\mu) \approx 1 - j/(n+1),
$$

therefore

$$
\log(1 - j/(n+1)) \approx \log(c) - X_{(j)}/\mu.
$$

In a *semi-log* plot (log transforming the probability points but not the
order statistics), when the distribution has exponential tails we obtain a
linear relationship with slope $-1/\mu$.

Using such probability tail plots, we can estimate distributional parameters
($\alpha$ or $\mu$) by fitting a least squares regression line to the points
in the plot. The number of points used in the least squares fit is a tuning
parameter that must be selected, typically in the range $20-200$. These
estimators are convenient, intuitive, and *distributionally robust* (since
they depend on the assumed form of the tail but do not require a complete
specification of a probability model). However these estimators may not be
very efficient (i.e. they may have high estimation variance). Several more
efficient estimators will be discussed below.

### The Hill estimate of the tail parameter

If a distribution has a power-law tail, we can solve for the upper quantiles,
yielding the quantile function

$$
Q(p) = (c/(1-p))^{1/\alpha}.
$$

Since the order statistics estimate quantiles, we have

$$
X_{(j)} \approx (c / (1 - j/(n+1)))^{1/\alpha}.
$$

An estimator known as the
[Hill estimator](https://en.wikipedia.org/wiki/Heavy-tailed_distribution#Hill's_tail-index_estimator)
begins by considering the ratios of upper order statistics

$$
X_{(n-j)} / X_{(n-k)} \approx ((j+1)/(k+1))^{1/\alpha}.
$$

Thus

$$
\log X_{(n-j)} / X_{(n-k)} \approx -\alpha^{-1} \log((j+1)/(k+1)).
$$

If we hold $k$ fixed and average these ratios for $1 \le j < k$, we get

$$
\hat{A} \equiv (k-1)^{-1}\sum_{j=1}^{k-1} \log X_{(n-j)} / X_{(n-k)} \approx -\alpha^{-1} \sum_{j=1}^{k-1}\log((j+1)/(k+1)).
$$

This establishes a relationship between the statistic $\hat{A}$ and the
quantity of intetrest $\alpha$. The constant of proportionality turns out to
be nearly equal to $1$:

$$
\sum_{j=1}^{k-1}\log((j+1)/(k+1)) \approx \int_0^k \log((x+1)/(k+1))dx \rightarrow -1.
$$

Therefore, the *Hill estimate* of the tail parameter is

$$
\hat{\alpha}_{\rm Hill} = 1/\hat{A}.
$$

In the Hill estimate, the value of $k$ is a tuning parameter. To select $k$,
we usually calculate $\hat{\alpha}_{\rm Hill}$ for various values of $k$
(typically $k \approx 20-200$) and choose a value that corresponds to a stable
range of values of the estimate.

If the data are exactly Pareto, the maximum likelihood estimate (MLE) of
$\alpha$ is the Hill estimate using $k=n$.

### The generalized Pareto distribution

As noted above, the one-parameter Pareto distribution may not fit many data
sets well, and further has an awkward sample space of $[1, \infty)$. To
address these issues, the
[generalized Pareto distribution](https://en.wikipedia.org/wiki/Generalized_Pareto_distribution)
was developed, which has sample space $[0, \infty)$ and complementary CDF

$$
P(X > t) = \sigma^{-1}(1 + t/(\sigma\alpha))^{-\alpha}.
$$

Note that the ratio of the CCDF's for the Pareto and generalized Pareto
satisfies

$$
\frac{\sigma^{-1}(1 + t/(\sigma\alpha))^{-\alpha}}{t^{-\alpha}} \rightarrow \sigma^{-1}(\sigma\alpha)^{\alpha} \ne 0,
$$

so the generalized Pareto and Pareto distributions both have power law tails
with the same tail parameter $\alpha$.

As $\alpha\rightarrow\infty$ (or equivalently $\xi=1/\alpha\rightarrow 0$),
the generalized Pareto distribution becomes the exponential distribution.

The famous
[Pickands-Balkema-De Haan theorem](https://en.wikipedia.org/wiki/Pickands%E2%80%93Balkema%E2%80%93De_Haan_theorem)
demonstrates that with appropriate choice of threshold $T$, the exceedances
for many distributions approximately follow a generalized Pareto distribution.
This theorem plays the role of the central limit theorem in the study of
extremes, since it allows us to use a specific parametric model to study data
that may follow a large range of distributions.

## Block maxima and the generalized extreme value distribution

Another way to approach extremes is to partition the data into blocks, and
calculate the maximum observed value within each block. According to the
[Fisher-Tippett-Gnedenko theorem](https://en.wikipedia.org/wiki/Fisher%E2%80%93Tippett%E2%80%93Gnedenko_theorem),
the distribution of these values should be well approximated by a
[generalized extreme value distribution](https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution),
(GEV), which is a three-parameter distribution. This is another example of a
central limit theorem-like result for extremes, since a wide variety of
populations have block maxima that are well-approximated by the GEV
distribution.

The block maxima approach is often used with serially observed data (time
series), and the block is a coarse resolution of time. For example, if our
time series consists of daily values, we might choose a block size of one
year. To estimate a GEV, we need to have a sufficient number of blocks. If
there are too few blocks, say 50 or fewer, then the GEV parameter estimates
will be very uncertain.

One advantage of working with block-wise maxima is that they are less
sensitive to *positive serial dependence* that causes clusters of extreme
values to occur in close proximity to each other. These clusters will
generally occur within one block and only the largest among them will
influence the analysis results. Fitting a generalized Pareto distribution to
the exceedances may produce biased estimates due to such dependence.

### Likelihood-based estimation

Likelihood-based estimation is generally more efficient than moment or
quantile-matching estimates such as considered above. The most well-known
likelihood-based estimator is the maximum likelihood estimator (MLE), which is
asymptotically fully efficient under standard conditions. However, the MLE can
be non-unique and difficult to compute. Additionally, in the case of the GEV
the support of the distribution depends on the parameters which violates the
standard conditions on which theoretical guarantees about the MLE are based.

For the standard one-parameter Pareto distribution, the MLE of the tail index
is simply

$$
\hat{\alpha}_{\rm MLE} = 1/({\rm Avg}(\log(X_i))).
$$

and the MLE of the $\mu$ in the exponential distribution is simply
$\hat{\mu} = \bar{X_i}$. These MLE's are not very useful in practice since few
populations can be well-approximated by these one-parameter distributions.

For families with more than a single parameter, the MLE is usually computed
numerically. Note that the generalized Pareto distribution has two parameters
and the generalized extreme value distribution has three parameters. Both of
these distributions are challenging to work with numerically. Alternative
likelihood-based estimators have been developed including the empirical Bayes
estimator discussed [here](https://www.jstor.org/stable/40586625).

Maximum likelihood estimates for the generalized extreme value (GEV)
distribution can be calculated numerically, but good starting values are
needed to obtain robust convergence. One way to get good starting values is
using the *probability weighted moments* approach discussed
[here](https://www.stat.cmu.edu/technometrics/80-89/VOL-27-03/v2703251.pdf).

### Return levels

The *m-observation return level* (or simply $m$-return) is a value $T$ that is
expected to be exceeded once out of every $m$ observations. In other words,
$T$ is the m-observation return level if $I_i = {\cal I}(X_i > T)$ and
$E[I_1 + \cdots + I_m] = 1$, where the $X_i$ are identically distributed
random variables. Since $E[I_1 + \cdots + I_m] = mE[I]$, where $I$ has the
same distribution as the $I_i$, the $m$-return can be inferred from the
equation $mE[I] = 1$ or $E[I] = 1/m$. Since
$E[I] = P(X > t) = 1 - P(X \le t)$, if $F(x)$ is the cumulative distribution
function (CDF) of a random variable $X$, then the m-observation return level
(for $m$ independent copies of $X$) is the solution to $F(x) = 1 - 1/m$. Thus,
the m-observation return level is the $1 - 1/m$ quantile of $X$.

If we are working with exceedances, we need to take account of the
observations that were excluded since they fell below the threshold $T$. Let
$q$ denote the proportion of the full dataset that exceeds $T$. Then the
m-observation return level is the $1 - 1/(q\cdot m)$ quantile calculated using
the distribution of exceedances.

## L-moments

A classical statistical moment is defined to be the expected value of a random
variable raised to a power. For example, $E[X]$ is the expected value and
$E[X^2]$ is the raw second moment. In practice we usually work with the
*centered moments*, for example $E[(X-EX)^2]$ is the centered second moment,
which is better known as the variance.

In principle, if you know all the moments of a distribution, then you know
everything about the distribution (there are some technical conditions for
this to be literally true, as it is based on the invertibility of the moment
generating function). But this fact is not very useful in practice because it
is nearly impossible to estimate high order moments $E[(X-EX)^k]$ for large
values of $k$. The sample estimator of this moment is
$n^{-1}\sum (X_i - \bar{X})^k$, and this estimator is consistent and
asymptotically unbiased, but if $k > 2$ it has huge mean squared error for
practically realistic sample sizes.

As noted earlier, many descriptive statistics are either moments or quantiles.
If high order moments are hard to estimate, perhaps there is a quantile-based
analogue to these moments that is easier to estimate? This idea led to the
development of [L-moments](https://en.wikipedia.org/wiki/L-moment) which are
linear combinations of order statistics (order statistics in turn are
estimates of quantiles).

The definition of an L-moment of arbitrary order is complex, so we focus here
on the first four L-moments.

The first L-moment $\lambda_1$ is the same as the usual mean.

The second L-moment of a distribution represented through the random variable
$X$ is defined as $\lambda_2 = (EX_{2:2} - EX_{1:2}) / 2$. Here, $X_{j:k}$ is
defined to be the random variable obtained by sampling $k$ independent values
from the distribution of $X$, and then taking the $j^{\rm th}$ largest among
them.

The third L-moment is

$$
\lambda_3 = (EX_{3:3} - 2EX_{2:3} + EX_{1:3}) / 3.
$$

Finally, the fourth L-moment is

$$
\lambda_4 = (EX_{4:4} - 3EX_{3:4} + 3EX_{2:4} + EX_{1:4}) / 4.
$$

Often we work with the standardized third and fourth L-moments,
$\lambda_3^s = \lambda_3/\lambda_2$ and $\lambda_4^s = \lambda_4/\lambda_2$.
Note that these standardized L-moments are *scale invariant* meaning that
their value is not changed by scaling the data. All L-moments except for the
first L-moment are *translation invariant*, meaning that their values are not
changed by adding a constant to all data values. Scale and translation
invariance (also known as *affine invariance*) are important becuase they
imply that the result does not depend on the units or origin of the
measurement scale.

L-moments are useful descriptive statistics that capture the shape of
distributions. They are more robust (less sensitive to contamination) than the
classical moments, and one can estimate much higher order L-moments than is
practical with classical moments.

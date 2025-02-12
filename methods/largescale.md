# Large scale inference

Modern data analysis often involves estimation of many parameters $\theta_j$,
and in some cases these parameters are weakly coupled so that it is possible
to estimate each parameter without referring to the others. For example, this
occurs when the $\theta_j$ are _local parameters_ referring to different
subpopulations of the overall population.

In many cases it is desirable to integrate the parameters $\theta_j$ to
understand the population as a whole. This contrasts with the strategy of
building a single all-encompassing model using the entire dataset. Situations
where the former strategy is useful commonly arise due to modern innovations
that allow for collection of data that are both high-dimensional (many
attributes measured on each unit) and heterogeneous (the units can be
clustered into subpopulations with greater heterogeneity between than within
subpopulations).

[Large scale inference](https://efron.ckirby.su.domains/other/2010LSIexcerpt.pdf)
is a body of statistical methods developed to accommodate data analysis in
this type of setting. In one common form, large scale inference facilitates
analysis conducted in two stages, where we first obtain a large collection of
local findings, e.g. parameter estimates $\hat{\theta}_j$ and their
accompanying standard errors $\hat{s}_j$. These results are then integrated to
obtain a more holistic understanding of the entire population. A central
statistic of inference is the Z-score $\hat{\theta}_j/\hat{s}_j$ that measures
evidence against the null hypothesis (when the null hypothesis is
$\theta_j=0$).

__Classical multiple inference__

Large scale inference has connections to the more classical area of *multiple
inference* (also known as *multiple testing* or *simultaneous inference*). The
usual setting for multiple inference is that we have $k$ hypotheses and
conduct a formal hypothesis test for each of them. This means that for the
j'th hypothesis we have a p-value $p_j$, and we make a decision to reject this
hypothesis based on the value of $p_j$, with smaller values of $p_j$ providing
more evidence against the null hypothesis. It is important to note that the
p-values are random variables, and if the j'th null hypothesis is true, then
$p_j$ follows a uniform distribution on the unit interval $(0, 1)$.

If the p-values are correctly calibrated and we only have one test ($k=1$),
then we typically proceed by choosing $0 < \alpha < 1$, and then reject the
test if $p_1 < \alpha$. This yields a *type-1 error rate* or *false positive
rate* of $\alpha$, which is often set to $\alpha=0.05$.

If $k > 1$ and we reject the j'th hypothesis whenever $p_j < \alpha$, then
each individual test has false positive rate equal to $\alpha$, but the
probability of rejecting at least one of the hypotheses incorrectly, known as
the
[family-wise error rate](https://en.wikipedia.org/wiki/Family-wise_error_rate),
is greater than $\alpha$. Thus, if people conduct multiple tests in the same
research project and only report the tests that have small p-values, the
evidence underlying the reported claims is exaggerated. This issue has been
known for a very long time, and recently has sometimes been termed
*p-hacking*.

An early method for addressing multiple hypotheses is the
[Bonferroni method](https://en.wikipedia.org/wiki/Bonferroni_correction),
based on the
[union bound](https://en.wikipedia.org/wiki/Boole%27s_inequality). Under the
Bonferroni method, to achieve a family-wise error rate of $\alpha$, the
individual hypotheses are carried out using a modified value of $\alpha$,
usually set to $\alpha/k$ (where $k$ is the number of tests). More generally,
the j'th test can be rejected when $p_j < \alpha_j$, with the $\alpha_j$
satisfying $\sum_j \alpha_j = \alpha$. Equivalently, we can multiply each
p-value by the number of tests, yielding adjusted p-values $k\cdot p_j$, and
assess them relative to an unadjusted $\alpha$.

A key consideration in multiple testing is the extent to which the p-values
are dependent random variables. The p-values are usually derived from the test
statistics $\hat{\theta}_j$, so if the test statistics are dependent, the
p-values will be as well. If the p-values are positively dependent, the
Bonferroni approach is conservative, and a number of approaches have been
devised to address this, which we do not consider here. If the p-values are
independent then the Bonferroni approach is only very slightly conservative.

__False Discovery Rates__

Multiple testing analysis based on the family-wise error rate takes the point
of view that incorrectly rejecting even a single null hypothesis is a bad
outcome whose probability should be controlled. This would be appropriate in
some settings. For example, suppose we are testing a medical treatment on 10
different human subpopulations (e.g. subpopulations based on demographic
characteristics). If the treatment has a serious adverse outcome in even one
of these subpopulations, the treatment as a whole may be considered
unsuitable.

However in many settings where we are aiming to make discoveries, the
family-wise error rate is arguably inappropriate. For example, suppose we are
trying to identify drug candidates using "screening experiments" that do not
involve human subjects. Any findings from this screening stage of research
will subsequently be tested on animals and then on humans. In this setting, it
may be appropriate to control the expected proportion of drug candidates
identified in the screening experiments that will go on to fail in human or
animal trials, say at 10%. Note that it is very different to identify a list
of drug candidates of which 10% are expected to be false positives, than to
identify a list of drug candidates such that there is a 0.1 probability of
even one of them being a false positive. The former perspective, discussed
further below, is known as
[false discovery rates](https://en.wikipedia.org/wiki/False_discovery_rate)
(FDR), while the latter is the family-wise error rate (FWER).

To formally define the FDR, suppose that we are interested in many hypotheses.
Let $H_j=1$ if the null hypothesis for the j'th hypothesis is true, and let
$R_j=1$ if we decide to reject the j'th hypothesis. The false positive rate,
or type-1 error rate for the j'th hypothesis, is $P(R_j=1 | H_j=1)$. The
family-wise error rate for all $k$ hypotheses is
$P(\sum H_jR_j > 0 | H_1, \ldots, H_k)$. The *false discovery rate* is
$E [\sum H_jR_j / \sum R_j | \sum R_j > 0]$. This is the expected proportion
among all rejected hypotheses that are falsely rejected.

There are several ways to estimate false discovery rates and we only briefly
introduce two of them here. A very basic approach is to start with a Z-scored
test statistic for the j'th hypothesis. This will often be the Wald test
statistics $z_j = \hat{\theta}_j / \hat{s}_j$. Under the null hypothesis,
$Z_j$ follows a standard normal distribution. If the null distribution is not
true, the expected value of $Z_j$ is non-zero, and the variance of $Z_j$
typically will either remain equal to 1 or become greater than 1. Suppose that
we decide to reject all null hypotheses for which $|Z_j| > T$, for some chosen
threshold $T$, that is, $R_j = I(|Z_j| > T)$. Further suppose that most of the
null hypotheses are false. In this case, the expected number of false
discoveries is upper bounded by $2k(1 - F(T))$, where $F$ is the standard
normal CDF (assuming that the reference distribution for the Z-scores is
standard normal). The estimated FDR at a given threshold $T$ is then
$2k(1 - F(T)) / \sum R_j$. We can adjust $T$ to achieve a desired FDR.

Another approach known as the *local FDR* works with densities instead of tail
probabilities. The local FDR at $t$ is the ratio $p_0f(t)/g(t)$, where $f$ is
the reference density for the Z-scores (often standard normal), $g$ is the
empirical density of all test statistics, and $p_0$ is the proportion of null
hypotheses that are true. A conservative setting for $p_0$ is $1$, but there
are ways to estimate $p_0$. A common approach for estimating $g$ is *Lindsey's
method* which involves binning the test statistics and estimating the log
probability in each bin (proportional to the density) using Poisson regression
with polynomial basis functions. There are also approaches employing an
"empirical null hypothesis" to estimate $f$ from the data, rather than relying
on asymptotic or other theoretical justifications.

__Other strategies__

At a high level, large scale inference can describe any approach to study a
collection of Z-scores or p-values that result from many analyses conducted on
different components of a single dataset. It can be very effective to consider
Z-scores in relation to auxiliary quantities such as the sample size of each
component. This is related to the concept of a "funnel plot" that considers
whether evidence against the null hypothesis (large Z-scores and small
p-values) mainly occurs when the sample size is large. In this case the
analysis may be limited by power, in which some of the partitions have too
little data to yield any findings. It can also be useful to plot Z-scores (or
p-values) against other auxiliary quantities that characterize the chosen data
components.

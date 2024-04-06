# Large scale inference

A useful strategy for statistical analysis can be to partition a dataset
so that relatively simple statistical analysis can be conducted on
each component of the data.  This
contrasts with the strategy of building a single model using the entire dataset.
This situation arises quite common due to modern innovations that allow for collection
of data that is both high-dimensional (many attributes measured on each unit),
and heterogeneous (the units can be clustered into subpopulations with much
greater heterogeneity between than within subpopulations).

As a concrete example, suppose we have a dataset with many observations and many
variables, e.g. with the observations being people in the US who live in the
50 American states.  For each person we observe their income, and we also know their
expenditure for 100 different spending categories (electricity bill, water bill, restaurant food,
etc.).  One possible analysis would be to build a single high-dimensional regression model, perhaps
using principal components regression or another approach for handling the high-dimensionality,
while simultaneously accounting for the hetereogeneity among the 50 states (e.g. using multilevel
regression).

Alternatively, we could pursue one of the following strategies:

* Regress income on each expenditure separately, accounting for state effects (e.g. by clustering the standard errors).
If we focus on the coefficient for each expenditure category, we have 100 estimates and their standard errors.

* Regress income on each expenditure separately, also stratifying by state.  Now we have $50\times 100 = 5000$
regression coefficients and standard errors.

In both cases, we get a large volume of results of the form $\hat{\theta}_j, \hat{s}_j$, where $\hat{\theta}_j$
is a point estimate and $\hat{s}_j$ is its estimated standard error.  It may also be possible to obtain
a p-value $p_j$ based on the reference distribution of $\hat{\theta}_j/\hat{s}_j$.  These quantities are usually not
fully independent between the hypotheses, but in many cases the dependencies may be weak.

*Large scale* inference is a body of statistical methods developed to accommodate analysis in this setting.
Specifically, large scale inference facilitiates analysis in two stages where we first obtain a large
collection of "local" estimates $\hat{\theta}_j$, $\hat{s}_j$, $\hat{p}_j$, then integrate these results
to obtain a more holistic understanding of the entire population.

_Classical multiple testing_

Large scale inference has connections to the more classical area of *mutliple inference* (also known as *simultaneous inference*).  The
canonical setting for multiple inference is that we have $k$ hypotheses and conduct a formal hypothesis
test for each of them.  This means that for the j'th hypothesis we have a p-value $p_j$, and we make a decision to
reject this hypothesis based on the value of $p_j$, with smaller values of $p_j$ providing more evidence
against the nyull hypothesis.  It is important to note that the p-values are random variables, and
if the j'th null hypothesis is true, then $p_j$ follows a uniform distribution on the unit interval $(0, 1)$.

If the p-values are correctly calibrated and we only have one test ($k=1$),
then we typically proceed by choosing $0 < \alpha < 1$, and then reject the test if $p_1 < \alpha$.  This yields a *type-1 error rate* or
*false positive rate* of $\alpha$, which is often set to $\alpha=0.05$.

If $k > 1$ and we reject the j'th hypothesis
whenever $p_j < \alpha$, then each individual test has false positive rate equal to $\alpha$, but the probability
of rejecting at least one of the hypotheses incorrectly, known as the *family-wise error rate*, is greater than $\alpha$.
Thus, if people conduct multiple tests in the same research project and only report the tests that have small
p-values, the evidence underyling the reported claims is exaggerated.  This issue has been known for a very long
time, and recently has sometimes been termed *p-hacking*.

An early method for addressing nultiple hypotheses is the *Bonferroni method*, based on the *union bound*.
Under the Bonferroni method, to achieve a family-wise error rate of $\alpha$, the individual hypotheses
are carried out using a modified value of $\alpha$, usually set to $\alpha/k$ (where $k$ is the number
of tests).  More generally, the j'th test can be rejected when $p_j < \alpha_j$, with the $\alpha_j$
satisfying $\sum_{j=1}^k\alpha_j = \alpha$.

A key consideration in multiple testing is the extent to which the p-values
are dependent random variables.  If the p-values are positively dependent, the Bonferroni approach is conservative,
and a number of approaches have been devised to address this, which we do not consider here. If the p-values
are independent then the Bonferroni approach is only very slighty conservative.

_False Discovery Rates_

Multiple testing analysis based on the family-wise error rate takes the point of view that incorrectly
rejecting even a single null hypothesis is a very bad outcome whose probabiity should be controlled.  This
would be appropriate in some settings.  For example, suppose we are testing a medical treatment on
10 different human subpopulations (e.g. subpopulations based on demographics).  If the treatment has a serious
adverse outcome in even one of these subpopulations, the treatment
as a whole may be considered unsuitable.

However in many settings where we are aiming to make discoveries,
the family-wise error rate is arguably inappropriate.  For example, suppose we are trying to identify
drug candidates using "screening experiments" that do not involve human subjects.  Any findings from this stage
of research will subsequently be tested on animals and then on humans.  In this setting, it may be
appropriate to control the expected proportion of drug candidates identified in the screening experiments that will
go on to fail in human or animal trials, say at 10%.  It is very different to identify a list of drug candidates
of which 10% are expected to be false positives, than to identify a list of drug candidates such that there is a
0.1 probability of even one of them being a false positive.  The former perspective, discussed
fiurther below, is known as _false discovery rates_ (FDR).

Suppose that we are interested in many hypotheses.  Let $H_j=1$ if the null hypothesis for the
j'th hypothesis is true, and let $R_j=1$ if we decide to reject the j'th hypothesis.  The
false positive rate, or type-1 error rate for the j'th hypothesis, is $P(R_j=1 | H_j=1)$.
The family-wise error rate for all $k$ hypotheses is $P(\sum H_jR_j > 0 | H_1, \ldots, H_k)$.  The
*false discovery rate* is $E [\sum H_jR_j / \sum R_j | \sum R_j > 0]$.  This is the expected proportion
among all rejected hypotheses that are falsely rejected.

There are several ways to estimate false discovery rates and we only briefly introduce them here.  A very basic
approach is to start with a Z-scored test statistic for the j'th
hypothesis.  This could be $Z_j = \hat{\theta}_j / \hat{s}_j$, or a transformed version of it intended to improve the normality
of $Z_j$.  Under the null hypothesis, $Z_j$ follows a standard normal distribution, and otherwise the expected value of $Z_j$
is non-zero, and the variance of $Z_j$ typically will either remain equal to 1 or become greater than 1.  Suppose that we decide
to reject all null hypotheses for which $|Z_j| > T$, for some chosen threhsold $T$, that is, $R_j = I(|Z_j| > T)$.  Further
suppose that most of the null hypotheses are false, then the expected number of false discoveries is upper bounded by
$2k(1 - F(T))$, where $F$ is the standard normal CDF (assuming that the referfence distribution for the Z-scores is
standard normal).  The estimated FDR is then $2k(1 - F(T)) / \sum R_j$.

Another approach known as the *local FDR* works with densities instead of tail probabilities.  The local FDR at
$t$ is the ratio $f(t)/g(t)$, where $f$ is the reference density for the Z-scores (often standard normal) and
$g$ is the empirical density of all test statistics.


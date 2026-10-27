#import "@preview/ilm:2.1.1": *

#show link: set text(fill: blue)

#set text(lang: "en")

#show: ilm.with(
  title: [Large scale Inference],
  authors: "Kerby Shedden",
  figure-index: (enabled: true),
  table-index: (enabled: true),
  listing-index: (enabled: true),
  chapter-pagebreak: false
)

= Introduction

Modern data analysis often involves estimation of many parameters which we will denote here $theta_j$, for $j=1,..., p$. In some cases these parameters are completely uncoupled or only weakly coupled, so that it is possible to estimate each parameter without referring to the others. For example, if we are interested in the median credit card debt for residents of each US state, we can estimate the value $theta_j$ for state $j$ without referring to the data from other states.

On the other hand, some parameters are deeply connected so that it is impossible to estimate one parameter without simultaneously estimating all the others. For example, suppose we wish to estimate the average blood pressure for a 50 year old, using data for people of age 18 to 80. Simplistically, we might use a linear model of the form $E[Y|A] = alpha + beta(A - 50)$, where $A$ is age and $Y$ is blood pressure. The parameter of interest here is $alpha$, but we cannot estimate $alpha$ without simultaneously estimating $beta$ -- these two parameters are strongly coupled (unless the average age happens to be close to $50$).

Even when it is possible to estimate the parameters $theta_j$ in isolation, in many cases it is desirable to integrate our understanding of these parameters to obtain a better understanding of the population as a whole. This contrasts with the strategy of building a single all-encompassing model using the entire dataset. Situations where the former strategy is useful commonly arise due to modern innovations that allow for collection of data that are both high-dimensional (many attributes measured on each unit) and heterogeneous (the units can be clustered into subpopulations with greater variation between than within subpopulations).

#link("https://efron.ckirby.su.domains/other/2010LSIexcerpt.pdf")[Large scale inference] is a body of statistical methods developed to accommodate data analysis in this type of setting. In one common form, large scale inference facilitates analysis conducted in two stages, where we first obtain a large collection of local findings, e.g. (approximately) unbiased parameter estimates $hat(theta)_j$ of $theta_j$ and their accompanying standard errors $hat(s)_j$. These results are then integrated to obtain a more holistic understanding of the entire population, often focusing on a distribution function $P_theta$ for the $theta_j$.

= Measuring evidence

A central statistic of interest is the Z-score $Z_j = frac(hat(theta)_j, hat(s)_j, style: "horizontal")$, that measures evidence against the null hypothesis for the $j^"th"$ parameter when the null hypothesis is $theta_j=0$.  The more general form $Z_j = frac((hat(theta)_j - theta_j^0), hat(s)_j, style: "horizontal")$ can be used when the null hypothesis is $theta_j = theta^0_j$ for a fixed $theta^0_j$.  A Z-score is a measure of evidence against the null hypothesis that is calibrated to a standard normal distribution (as long as the parameter estimate is asymptotically and approximately normal).

When finite sample behavior persists or there are large estimation errors in the $hat(s)_j$, we might anticipate that the null distribution of $Z_j$, rather than being Gaussian, is better approximated by a different distribution, such as a Student $t$ distribution.  In that case we can map $Z_j$ to a standard normal distribution using the transformation $Z_j -> F_Z^(-1)(F_(T)(Z_j))$.

Under standard conditions, $Z_j$ is approximately equal in distribution to $sqrt(n)frac((theta_j - theta^0_j), s_j^*, style: "horizontal") + E_j$, where $E_j$ follows a standard normal distribution, and $s_j^*$ is the asymptotic standard deviation of $hat(theta)_j$, i.e. $sqrt(n)hat(s_j) -> s_j^*$.  Thus, if the null hypothesis is true ($theta_j = theta^0_j$), then $Z_j$ is standard normal, and if the null hypothesis is false, then $Z_j$ diverges at rate $sqrt(n)$.

The preceding analysis also shows that a Z-score combines information about effect size and precision/significance.  The value of $frac((theta_j - theta^0_j), s_j^*, style: "horizontal")$ is a standardized effect size.  It measures how far the true parameter $theta_j$ falls from the hypothesized null parameter $theta^0_j$ in units of the asymptotic standard error $s_j^*$.  The factor $sqrt(n)$ accounts for statistical power, so that the same standardized effect size produces a larger Z-score when the sample size (power) is greater.

An alternative very popular measure of evidence against the null hypothesis is a #link("https://en.wikipedia.org/wiki/P-value")[p-value], which is a statistic that follows a uniform distribution on $[0, 1]$ when the null hypothesis is true.  There is often a direct correspondence between p-values and Z-scores. A Z-score can be converted to a one-sided p-value using the transformations $F(-Z)$ or $F(Z)$, where $F$ is the standard normal CDF.  A Z-score can be converted to a two-sided p-value using the transformation $2 dot.c F(-|Z|)$.  While Z-scores can measure either one-sided or two-sided evidence, a p-value is always one-sided.

Recently, the notion of an #link("https://en.wikipedia.org/wiki/E-values")[E-value] has been popularized, where $T$ is an E-value if $E[E] <= 1$ when the null hypothesis is true.  Since E-values only constrain the expectation under the null, they are more robust and have some other desirable properties for purposes such as combining dependent sources of evidence.

In summary, a Z-score is a measure of evidence against the null hypothesis that follows a standard normal distribution when the null hypothesis is true, and diverges when the null hypothesis is false.  A p-value is a measure of evidence against the null hypothesis that follows a uniform distribution on $[0, 1]$ when the null hypothesis is true, and becomes concentrated near zero when the null hypothesis is false.

= Classical multiple inference

Large scale inference has connections to the more classical area of _multiple inference_ (also known as _multiple testing_ or _simultaneous inference_). The usual setting for multiple inference is that we have $k$ hypotheses and conduct a formal hypothesis test for each of them. This means that for the $j^"th"$ hypothesis we have a p-value $p_j$, and we make a decision to reject this hypothesis based on the value of $p_j$, with smaller values of $p_j$ providing more evidence against the null hypothesis. It is important to note that the p-values are random variables, and if the $j^"th"$ null hypothesis is true, then $p_j$ follows a uniform distribution on the unit interval $(0, 1)$.

If the p-values are correctly calibrated and we only have one test ($k=1$), then we typically proceed by choosing $0 < alpha < 1$, and then reject the test if $p_1 < alpha$. This yields a _type-1 error rate_ or _false positive rate_ of $P_0(p_1 < alpha) = alpha$, which is often set to $alpha=0.05$.

If $k > 1$ and we reject the $j^"th"$ hypothesis whenever $p_j < alpha$, then each individual test has false positive rate equal to $alpha$, but the probability of rejecting at least one of the hypotheses incorrectly, known as the #link("https://en.wikipedia.org/wiki/Family-wise_error_rate")[family-wise error rate], is greater than $alpha$. Thus, if people conduct multiple tests in the same research project and only report the tests that have small p-values, the evidence underlying the reported claims is exaggerated. This issue has been known for a very long time, and recently has sometimes been termed _p-hacking_.

An early method for addressing multiple hypotheses is the
#link("https://en.wikipedia.org/wiki/Bonferroni_correction")[Bonferroni method], based on the #link("https://en.wikipedia.org/wiki/Boole%27s_inequality")[union bound]. Under the Bonferroni method, to achieve a family-wise error rate of $alpha$, the individual hypotheses are carried out using a modified value of $alpha$, usually set to $frac(alpha, k, style: "horizontal")$ (where $k$ is the number of tests). More generally, the $j^"th"$ test can be rejected when $p_j < alpha_j$, with the $alpha_j$ satisfying $sum_j alpha_j = alpha$. Equivalently, we can multiply each p-value by the number of tests, yielding adjusted p-values $k dot.c p_j$, and assess them relative to an unadjusted probability threshold $alpha$.

A key consideration in multiple testing is the extent to which the p-values are dependent random variables. The p-values are usually derived from the test statistics $hat(theta)_j$, so if the test statistics are dependent, the p-values will be as well. If the p-values are positively dependent, the Bonferroni approach is conservative, and a number of approaches such as the step-down approaches have been devised to address this, which we do not consider here. If the p-values are independent then the Bonferroni approach is only very slightly conservative.

= False Discovery Rates

Multiple testing analysis based on the family-wise error rate takes the point of view that incorrectly rejecting even a single null hypothesis is a bad outcome whose probability should be controlled. This would be appropriate in some settings. For example, suppose we are testing a medical treatment on 10 different human subpopulations (e.g. subpopulations based on demographic characteristics). If the treatment has a serious adverse outcome in even one of these subpopulations, the treatment as a whole may be considered unsuitable.

However in many settings where we are aiming to make discoveries, the family-wise error rate is arguably inappropriate. For example, suppose we are trying to identify drug candidates using "screening experiments" that do not involve human subjects. Any findings from this screening stage of research will subsequently be tested on animals and then on humans. In this setting, it may be appropriate to control the expected proportion of drug candidates identified in the screening experiments that will go on to fail in human or animal trials, say at 10%. Note that it is very different to identify a list of drug candidates of which 10% are expected to be false positives, than to identify a list of drug candidates such that there is a 0.1 probability of even one of them being a false positive. The former perspective, discussed further below, is known as #link("https://en.wikipedia.org/wiki/False_discovery_rate")[false discovery rates] (FDR), while the latter is the family-wise error rate (FWER).

To formally define the FDR, suppose that we are interested in many hypotheses. Let $H_j=1$ if the $j^"th"$ null hypothesis is true, and let $R_j=1$ if we decide to reject the $j^"th"$ null hypothesis. The false positive rate, or type-1 error rate for the $j^"th"$ null hypothesis, is $P(R_j=1 | H_j=1)$. The family-wise error rate for all $k$ hypotheses is $P(sum H_j R_j > 0 | H_1, ..., H_k)$. The _false discovery rate_ is

$
E [sum frac(H_j R_j, sum R_j, style: "horizontal") | sum R_j > 0].
$

This is the expected proportion among all rejected hypotheses that are falsely rejected.

A basic approach for estimating the FDR is to start with a collection of Z-scored test statistics for the $k$ hypotheses. These will often be the Wald test statistics of the form $z_j = frac(hat(theta)_j, hat(s)_j, style: "horizontal")$. Under the null hypothesis, each $Z_j$ follows a standard normal distribution. If the null distribution is not true, the expected value of $Z_j$ is non-zero, and the variance of $Z_j$ typically will either remain equal to 1 or become greater than 1. Suppose we decide to reject all null hypotheses for which $|Z_j| > T$, for some chosen threshold $T$, that is, $R_j = I(|Z_j| > T)$. Further suppose that most of the null hypotheses are false. In this case, the expected number of false discoveries is upper bounded by $2k(1 - F(T))$, where $F$ is the standard normal CDF (assuming that the reference distribution for the Z-scores is standard normal). The estimated FDR at a given threshold $T$ is then $frac(2k(1 - F(T)), sum_j R_j, style: "horizontal")$. We can adjust $T$ to achieve a desired FDR.

A related but distinct notion is that of the _local FDR_, which works with densities instead of tail probabilities. The local FDR at $t$ is the ratio $frac(p_0f(t), g(t), style: "horizontal")$, where $f$ is the reference null density for the Z-scores (often standard normal), $g$ is the density of all test statistics, and $p_0$ is the proportion of null hypotheses that are true.  The density $g$ is a two-component mixture of $f$ and the density of $Z$-scores for non-null hypotheses. A conservative setting for $p_0$ is $1$, but there are ways to estimate $p_0$ if the non-null proportion is thought to be large.

A common approach for estimating $g$ is _Lindsey's method_ which involves binning the test statistics and estimating the log probability in each bin (proportional to the density) using Poisson regression with polynomial basis functions. There are also approaches employing an #link("https://www.tandfonline.com/doi/abs/10.1198/016214504000000089")[empirical null hypothesis] to estimate $f$ from the data, rather than relying on asymptotic or other theoretical justifications.

The conventional FDR targets a characteristic of a set of claims ${j: Z_j>T}$ based on test statistics exceeding a given threshold.  At a population level it is the proportion of this list that is false.  The conventional FDR is straightforward to estimate and there are non-asymptotic guarantees about its performance.  However, due to integrating over the tail, some hypotheses where $Z_j$ is close to the threshold $T$ may have a small FDR due mainly to the influence of other members of the set that have large Z-scores, even if they individually are much less likely to be true non-nulls.

In contrast, the local FDR estimates quantities that specifically refer to each individual hypothesis $Z_j$.  If the local FDR for the $j^"th"$ hypothesis is small, then that specific hypothesis is very unlikely to be null. Since the local FDR involves density estimation and these estimates are only correct asymptotically  there are no finite sample guarantees about its performance.

= Deconvolution

Suppose we estimate a collection of parameters $(theta_j)_(j=1)^m$ using estimates $(hat(theta)_j)_(j=1)^m$ with (estimated) standard errors $(hat(s)_j)_(j=1)^m$.  In many settings we can adopt the model $hat(theta)_j = theta_j + s_j E_j$, where $E_j$ follows a standard normal distribution given $theta_j$.  Here we are treating the $theta_j$ as a random sample from a distribution $P_(theta)(dot.c)$, which becomes the primary object of interest.

The distribution of the $hat(theta)_j$, denoted $P_(hat(theta))(dot.c)$ is a convolution of $P_theta$ with (heteroscedastic) Gaussian estimation errors.  The act of recovering $P_theta$ from the empirical distribution $hat(P)_(hat(theta))$ and knowledge of the estimation standard errors $hat(s)_j$ is a form of _deconvolution_, and is a canonical example of a hard inverse problem.  Nevertheless, there are methods for recovering $P_theta$ all or in part.

One approach that can work well is to approximate $P_theta$ as a finite mixture over a collection of grid points $g_1, ..., g_m$.  Thus, we have the approximation

$
P_theta approx sum_(j=1)^m pi_j bb(1)_(g_j).
$

This now becomes a finite mixture model where the _emission distributions_ $P(hat(theta)_j | theta=g_k) = N(g_k, s^2)$ are known, so we only need to estimate the probabilities $(pi_j)_(j=1)^m$.  This can be accomplished in many ways, but one straightforward approach is to use the EM algorithm, optionally with entropic regularization to partially counter the ill-posedness of the inverse problem.

= Empirical likelihood methods

While the deconvolution approaches give useful point estimates, they are difficult to use when uncertainty assessment is desired.  A powerful set of approaches for uncertainty assessment in these settings are the empirical likelihood methods, especially the variant known as _jackknife empirical likelihood_.

Let $hat(theta)_j$ denote a parameter estimate based on an iid sample of data, and let $hat(theta)_(j)[-i]$ denote the estimate with observation $i$ omitted.  There are generalizations of this approach that can be used when the data are not iid, but we do not discuss them here.  The _pseudo-observations_ (or _pseudo-values_) are

$
eta_j = n hat(theta) - (n-1)hat(theta)[-i].
$

Now suppose that we place an empirical distribution on the $eta_j$.  That is, we have probabilities $pi_j$ such that $sum pi_j bb(1)_(eta_j)$ is an empirical distribution. For any candidate value $theta^*$ of $theta$, we can estimate the $pi_j$ using a constrained maximum entropy approach. Impose the constraints

- $pi_j >= 0$ for $j=1, ..., m$
- $sum_j pi_j = 1$
- $sum_j pi_j eta_j = theta^*$.

Then, maximize $sum_j log(pi_j)$ subject to these constraints, yielding estimates $hat(pi)_(j)(theta^*)$.  Let $L(theta^*) = sum_j log(hat(pi)_(j)(theta^*))$, and let $hat(theta) = "argmax"_theta L(theta)$.  Finally, we can conduct a non-parametric likelihood ratio test employing a Wilk's type theorem that $2(L(theta) - L(theta^*))$ follows a $chi^2_1$ distribution.  This allows us to obtain a p-value for the null hypothesis that $theta = theta^*$, and to obtain confidence sets by inverting these tests.

= Exploratory analysis of Z-scores

Large scale inference can describe any approach for studying a collection of Z-scores or p-values that result from many analyses conducted on different components of a single dataset. Here are some examples of such analyses, focusing on visualization-based approaches.

- Plot the Z-scores against sample sizes, variance inflation factors, design effects, or any other parameter that captures how powerful the data are for estimating each parameter of interest. This is related to the concept of a "funnel plot" that considers whether evidence against the null hypothesis (large Z-scores and small p-values) mainly occurs when the sample size is large. In this case the analysis may be limited by power, in which some of the partitions have too little data to yield any findings.
- Plot Z-scores derived from two related statistics against each other. Since the two statistics will have good power against different alternatives, this can be used to reveal whether the setting in-hand is enriched for alternatives of a specific type.
- Plot Z-scores for paired subpopulations against each other.  For example, if we have $n$ traits and Z-scores for $Z^f_j$, $Z^m_j$ for female and male subjects, respectively, it can be informative to plot $Z^f_j$ against $Z^m_j$.

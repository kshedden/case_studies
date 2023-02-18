# Power analysis and study design

Empirical research is the process of using data to address
questions about the natural world.  Empirical research always
starts with a well-formulated question.  After sepcifying
a question, one can consider
what data and analytic methods may be used to address the question.

Designing an empirical research study involves many factors
including ethical constraints, costs, and logistical matters.
A key consideration in any research study
is *statistical power*.  Statistical power is often
considered when research is conducted in the setting
of a falsifiable null hypothesis.
Nearly always, the null hypothesis is a "straw man" that we want
to refute (reject).  We do not know whether the null
hypothesis is false or true, therefore, we wish to conduct our study
according to the following two principles: (i) if the
null hypothesis is true we wish to control the probability
of incorrectly asserting it to be false; (i) if the null
hypothesis is false, we want to have a high probability
of rejecting it.

Before proceeding, note that there are settings where we
wish to consider statistical power, but the research does
follow the "null hypothesis significance testing"
paradigm.  For example, we may wish to design a study to
estimate a parameter of interest, and we want the precision
of the estimate to be such that a confidence interval has
a certain width.

## Paired t-test

Most of the basic ideas can be considered in the setting of the
paired t-test.  Suppose we observe paired data $(x_i, y_i)$
for each of $n$ subjects, and the observations are made
independently and follow the same distribution (i.e. we have IID
data).  Let $\mu_x, \mu_y$ denote the population means of
$x$ and $y$, and let $\sigma_x$, $\sigma_y$ denote the corresponding
standard deviations.  Further, let $\rho$ denote the correlation
between $x$ and $y$.

A paired t-test considers the differences $d_i = x_i - y_i$
and conducts a one-sample test of the null hypothesis
$E[d] = 0$.  Under the null hypothesis, the $d_i$ have
expected value 0 and variance equal to
$\tau^2  = \sigma_x^2 + \sigma_y^2 - 2\rho\sigma_x\sigma_y$.
We can form a test statistic using the Z-score $T \equiv \sqrt{n}\bar{d}/\hat{\tau}$,
where $\hat{\tau}$ is the sample standard deviation of the
$d_i$.

Under the null hypothesis, $T$ approximately follows a standard
normal distribution, although if the sample size is small, the Student-t
distribution with $n-1$ degrees of freedom might provide a
better approximation.  For now we use the normal reference
distribution.  If $T$ is standard normal then we can reject
the null hypothesis when $T > 2$, focusing without loss of
generality on alternatives in which $\mu_x > \mu_y$. This
test approximately controls the type-I error rate at 5\%.

To consider the power, suppose that $E[d] = \theta$,
i.e. the effect size is $\theta / \tau$. The power (treating
$\tau$ as known) is:

$$
P(T > 2) = P((\bar{d} - \theta)/\tau > 2 - \theta/\tau).
$$

Under the alternative hypothesis, $(\bar{d} - \theta)/\tau$ follows a
standard normal distribution.
Conventionally (but somewhat arbitrarily) we seek 80% power.
Since the 0.2 quantile of a standard normal distribution is -0.84, we
achieve 80\% power when $\theta/\tau = 2.84$. That is, the effect
size must be at least 2.84 to achieve 80\% power.  If sample sizes are small,
then the standard normal distribution is underdispersed relative
to the true sampling distribution, so it would be better to replace 2.84 with 3
to obtain a good rule of thumb for 80\% power.  In words, the
effect must be three times the standard error to have good power.

This analysis reveals which factors influence the power.  For example,
the power is greater when $\theta/\tau$ is greater.  This in turn
happens when $\rho$ is greater.  Thus, in this setting we benefit
from strong correlations between $x$ and $y$.

To conduct a thorough power analyis in this setting, we would
consider plausible values for $n$, $\mu_x - \mu_y$, $\sigma_x$,
$\sigma_y$, and $\rho$.  Then we can assess which combinations
of these values yield power that is deemed adequate.


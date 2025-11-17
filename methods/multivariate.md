# Additional methods for multivariate analysis

This document discusses several modern methods for analyzing multivariate
data. These methods were developed since around 1970, whereas the classical
multivariate methods such as PCA, CCA, etc. were developed in the 1930's or
earlier. The methods discussed here make less use of factor analysis/dimension
reduction using the singular value decomposition, and instead make more use of
point set geometry. Methods such as PCA work best for
[elliptically distributed](https://en.wikipedia.org/wiki/Elliptical_distribution)
data, whereas the methods discussed here can in some cases overcome this
limitation.

## Functional data

A particular type of multivariate data is known as *functional data*, in which
we observe vectors $v \in {\mathbb R}^d$ that arise from evaluating a function
on a grid of points, i.e. $v = [f(t_i)]_{i=1}^d$ for a grid
$t_1 < t_2 < \cdots < t_d$. If the function $f$ is smooth then the elements of
each $v_i$ will reflect this smoothness.
[Functional Data Analysis (FDA)](https://en.wikipedia.org/wiki/Functional_data_analysis)
encompasses many methods for analyzing functions as data. In practice we never
actually observe a function in its entirety, thus the data we work with in FDA
are finite dimensional vectors, and (superficially) have the same form as
other types of quantitative multivariate data. But since the data are
considered to arise by evaluating smooth functions, we can view the population
of interest as consisting of functions, even while the data are vectors. A
distinct set of methods have been developed to take advantage of this
property.

## Data Depth

There are various ways to measure the [depth](https://arxiv.org/abs/1207.4988)
of a point $z \in {\mathbb R}^d$ relative to a distribution or collection of
points $\\{x_i\in {\mathbb R}^d; i=1,\ldots,n\\}$. Formally, a *depth measure*
on ${\mathbb R}^d$ is a function from ${\mathbb R}^d\rightarrow{\mathbb R}^+$
that quantifies the depth of each point as a non-negative real number.

"Deep" points are surrounded in all directions by many other points, while
"shallow" points lie near the surface or exterior of the point set. Deep
points are often referred to as having high "centrality" while shallow points
have low centrality or high "outlyingness". Data depth can be viewed as a
multivariate generalization of the notion of a quantile, with the deepest
point of a distribution being a type of multivariate median.

Below are several examples of depths.

### Halfspace depth

The original definition of depth was the
[halfspace depth](https://en.wikipedia.org/wiki/Tukey_depth) introduced by
John Tukey in 1975. The definition of the halfspace depth is simple to
describe graphically and a bit more difficult to define formally. To calculate
the halfspace depth of a single point $z\in {\mathbb R}^d$ with respect to a
collection of points $\{x_i; i=1, \ldots, n\}$, with each
$x_i \in {\mathbb R}^d$, let $U$ denote the set of all unit vectors in
${\mathbb R}^d$ and define the halfspace depth as

$$
D_{HS}(z; \\{x_i\\}) = {\rm min}_{u\in U} \sum_i {\cal I}(u^T (x_i - z) > 0)
$$

What we are doing here is searching for a hyperplane passing through $z$ that
places the greatest fraction of the $x_i$ on one side of the hyperplane. If
$z$ falls at the geometric center of a collection of symmetrically distributed
points, then $z$ is as deep as possible and will have halfspace depth
approximately equal to 1/2. At the other extreme there is a line passing
through $z$ such that all of the $x_i$ are on the same side of this line. In
this case the point $z$ is as shallow as possible and its halfspace depth is
approximately equal to zero.

The halfspace depth is geometrically natural but expensive to compute exactly
except in two dimensions.

### Spatial depth

The spatial depth has a simple definition that is relatively easy to compute
in high dimensions:

$$
D_S(z; \\{x_i\\}) = 1 - \lVert {\rm Avg}_i \\{(x_i-z)/ \lVert x_i-z \lVert \\} \lVert
$$

Note that $(x_i-z) / \lVert x_i-z \lVert$ is a unit vector pointing in the
direction from $z$ to $x_i$. If a point $z$ is "shallow" then most of the unit
vectors $(x_i-z) / \lVert x_i-z \lVert$ point in roughly the same direction,
and therefore their average value will have large magnitude. If a point $z$ is
"deep" then these unit vectors will point in many different directions and
their average value will have small magnitude.

### L2 depth

The $L_2$ depth also has a simple definition and is easy to compute:

$$
D_{L_2}(z; \\{x_i\\}) = 1 / (1 + {\rm Avg}_i\\{\lVert x_i-z \lVert \\}).
$$

### Band depth

Another way to measure depth that is especially useful for vectors that are
sampled functions is the *band depth*. In one common version of band depth, to
calculate the depth of $z$ relative to the reference set $\\{x_i\\}$, we
consider all triples of three distinct values in the reference set
$\\{x_i\\}$, say $x_{j_1}$, $x_{j_2}$, $x_{j_3}$ where $j_1 < j_2 < j_3$. For
each such triple we calculate the proportion of indices $k$ such that $z(k)$
lies between the maximum and minimum of $x_{j_1}(k)$, $x_{j_2}(k)$,
$x_{j_3}(k)$. These values are then averaged over all triples of points in the
reference set.

### Properties of a good depth function

Analysis based on depths does not directly rely on probability models, making
it quite distinct from many other methods of statistical data analysis. For
statistical methods based on probability there are standard properties such as
unbiasedness, consistency, accuracy, and efficiency that are used to evaluate
the performance of the approach. Although it is possible to place depth into a
probabilistic framework so that these notions can be applied, several
researchers have attempted to define the geometric properties that a depth
function should exhibit that do not depend on any probability framework. Four
basic such properties are

- *Affine invariance* -- If the data are transformed by the affine orthogonal
  mapping $x\longrightarrow c + Qx$, where $c\in {\mathbb R}^d$ is a fixed
  vector and $Q$ is an orthogonal matrix, then the depths do not change.

- *Maximality at the center* -- If the data are symmetric around zero, i.e. if
  $-x$ is in the dataset whenever $x$ is in the dataset, then the vector $0_d$
  achieves the maximum depth.

- *Monotonicity relative to the deepest point* -- Let $\tilde{x}$ be the
  deepest point and we consider any unit vector $u$, and we then evaluate the
  depth at $\tilde{x} + \lambda u$ for $\lambda \in {\mathbb R}^+$, then the
  depth is a decreasing function of $\lambda$.

- *Vanishing at infinity* -- for any sequence $z_i$ with $\lVert z_i \lVert$
  tending to infinity, the depths of the $z_i$ tend to zero.

### Depth peeling

Data depth can be used in exploratory multivariate analysis to identify the
most central or typical points and then contrast them with the more outlying
points. A systematic way to do this is to stratify the data based on depth and
then inspect the points in each depth stratum. For example, if we stratify the
data into 10 groups based on depth deciles, the first decile consists of the
shallowest 10% of points and the last decile consists of the deepest 10% of
points.

Often (not always) there is little heterogeneity in the deepest decile,
meaning that all of the deepest points are very similar. However there is
nearly always heterogeneity in the shallowest decile, as there are many
different ways to be near the periphery of a collection of points.

## Quantization

A quantization algoithm aims to represent a multivariate distribution through
a small number of representative points. This can be a useful exploratory
technique if the distribution being studied has a complex form that is not
approximately Gaussian or elliptical, and is not well captured through
independent additive factors (as in PCA). The goal of almost any quantization
algorithm is to find a collection of representative points $\\{x_i\\}$
corresponding to the distribution of a random variable $Y$ that are optimal in
some sense. Inspecting the representative points may provide a quick means to
understand the structure of the distribution.

[Support points](https://arxiv.org/abs/1609.01811) are one effective form of
quantization. To understand the idea behind support points, suppose that we
are given a distribution function $G$ on ${\mathbb R}^d$ that we wish to
approximate with a finite set of points. Let $Y$ denote a random draw from
$G$. Now consider an approximating distribution $F$ with random draw $X$.

We are given $G$ and wish to construct $F$ to approximate $G$, so we begin by
defining a distance function that measures how far apart $F$ and $G$ are from
each other. Note that this distance compares two probability distributions (it
is not the simpler notion of a distance between vectors). Distances among
probability distributions play an important role in modern statistics. When
$d=1$, many natural distance measures on probability distributions can be
constructed, but it is harder to construct good distances on probability
distributions when the dimension $d$ is greater than one.

One distance measure on probability distributions that turns out to be very
effective and relatively easy to work with is called the
[energy distance](https://en.wikipedia.org/wiki/Energy_distance), defined as

$$
2E \lVert X-Y \lVert - E \lVert X-X^\prime \lVert - E \lVert Y-Y^\prime \lVert.
$$

Here, $X$, $X^\prime$ are independent draws from $F$, and $Y$, $Y^\prime$ are
independent draws from $G$. It turns out that expression above is equal to
zero if and only if $F \equiv G$. This is a necessary property for a distance
to have.

The interpretation of the expression above is that $F$ is close to $G$ if (i)
a random draw from $F$ tends to be close to a random draw from $G$, (ii) two
independent random draws from $F$ are far from each other, and (iii) two
independent random draws from $G$ are far from each other.

Our goal is to approximate a given distribution $G$ with a simpler
distribution $F$ that we construct. Here "simple" means that $F$ has finite
support, i.e. a finite sample space. Since $G$ is given, the term
$E \lVert Y-Y^\prime \lVert$ in the energy distance is fixed and can be
ignored when constructing $F$. Thus, our goal is to construct $F$ that
minimizes

$$
2E \lVert X - Y \lVert - E \lVert X - X^\prime \lVert.
$$

It is worth considering an alternative approach in which we simply minimize
the first term above, $E \lVert X-Y \lVert$. However doing this always yields
a degenerate solution in which $F$ places all of its probability mass on the
[spatial median](https://en.wikipedia.org/wiki/Geometric_median), which is the
fixed vector $V$ that minimizes $E \lVert Y - V \lVert$. This is the reason
that the "repulsive" term $E \lVert X-X^\prime \lVert$ in the distance measure
is essential.

In practice, we do not observe the distribution $G$ but instead observe a
sample $y_1, \ldots, y_N$. As stated above, the approximating distribution $F$
that we are constructing is supported on a finite set of points
$x_1, \ldots, x_n$. This leads us to the empirical analogue of the distance
function above:

$$
\frac{2}{nN}\sum_{i=1}^n\sum_{j=1}^N \lVert y_j - x_i \lVert - \frac{1}{n^2}\sum_{i=1}^n\sum_{j=1}^n \lVert x_i-x_j \lVert.
$$

Our goal here was to discuss the motivation behind the support point
algorithm. We will not proceed further with discussion of the process of
numerically minimizing this function (see the paper linked above for
computational details).

## Minimum covariance determinant

Much of classical multivariate analysis is based on the
[covariance matrix](https://en.wikipedia.org/wiki/Covariance_matrix) of
several variables. The most common way to estimate this matrix is the
[sample covariance matrix](https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices),
which (but for a minor scaling difference) coincides with the maximum
likelihood estimate for Gaussian data (note however that the data do not need
to follow a Gaussian distribution for the sample covariance matrix to be
informative). The sample covariance matrix has some desirable properties
including being unbiased and elementwise consistent. However for larger
dimensions or smaller sample sizes, the sample covariance matrix can perform
poorly. Also, the sample covariance matrix is not robust to outliers.

To address these limitations, many alternative estimators of the covariance
matrix have been devised. Here we describe one of them, the
[minimum covariance determinant](https://arxiv.org/abs/1709.07045), or MCD. We
will not give full details here, but the main idea is to select a constant
$m \ge n/2$, where $n$ is the sample size, and define the MCD estimate to be
the sample covariance matrix of $m$ observations that has the smallest
determinant among all such subsets.

The determinant of a covariance matrix is a measure of its dispersion. If $x$
follows the multivariate normal distribution $N_p(\mu, \Sigma)$, the squared
Mahalanobis distance $(x-\mu)^T \Sigma^{-1}(x-\mu)$ follows the $\chi^2_p$
distribution, where $x \in {\mathbb R}^p$. It follows that the ellipsoid

$$
\\{(x-\mu)^T \Sigma^{-1}(x-\mu) \le t\\}
$$

covers $F_p(t)$ probability mass of the $N_p(\mu, \Sigma)$ distribution, where
$F_p$ is the CDF of the $\chi^2_p$ distribution. Setting $t = F_p^{-1}(0.95)$,
for example, gives an ellipsoid that covers 95% of the probability. The volume
of this ellipsoid is $t^p\cdot|\Sigma|^{1/2}$ times the volume of the unit
ball, which is $\pi^{p/2}/\Gamma(p/2+1)$.

The rationale for the MCD estimator is that if a few observations greatly
increase the determinant, they may be distorting the covariance estimate. As
with any "robust" estimator, there is an implicit assertion that certain
points may be outliers and should be discounted or excluded from the analysis.
However there is no guarantee that any particular definition of an outlier is
correct. In the case of the MCD, we are asserting that points that inflate the
volume of an estimated _high probability ellipsoid_ may be outliers.

Exact calculation of the MCD estimator of the covariance matrix is
challenging. Enumerating all possible subsets of $m$ out of $n$ observations
is impossible in practice. Therefore a greedy approach called _fast MCD_ is
used in practice.

One use for the covariance matrix is to quantify the "outlyingness" of
individual observations. Specifically, let $x\in{\mathbb R}^d$ denote an
observation, $\mu\in{\mathbb R}^d$ denote the mean, and
$\Sigma\in{\mathbb R}^{d\times d}$ denote the covariance matrix. The squared
[Mahalanobis distance](https://en.wikipedia.org/wiki/Mahalanobis_distance) is
defined to be $(x-\mu)^T \Sigma^{-1}(x-\mu)$. Points with larger Mahalanobis
distance to the center are "less central" or "more outlying". Using the MCD or
another robust estimate of the covariance matrix (and of the mean) can
sometimes reveal more interesting sets of outliers than those that are
identified by the sample covariance matrix.

There is a connection between the MCD and data depth using the "Mahalanobis
depth" (MD), which is defined to be $1 / (1 + (x-\mu)^T \Sigma^{-1}(x-\mu))$.
The conventional MD uses the sample covariance to estimate $\Sigma$, but this
opens the possibility of "masking", where a small group of outliers distorts
the covariance matrix so as to mask their outlyingness. Using the MCD to
estimate $\Sigma$ circumvents this issue.

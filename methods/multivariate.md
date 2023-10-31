# Various methods of multivariate analysis

This document discusses several useful methods for analyzing
multivariate data that are less widely known than the classical
multivariate methods such as PCA, CCA, etc.

## Functional data

A particular type of multivariate data is known as *functional data*,
in which we observe vectors $v \in {\mathbb R}^d$ that arise from evaluating
a function on a grid of points, i.e. $v = [f(t_i)]_{i=1}^d$ for
a grid $t_1 < t_2 < \cdots < t_d$.  If the function $f$ is smooth
then the elements of each $v_i$ will reflect this smoothness.
*Functional Data Analysis* (FDA) encompasses many methods for
analyzing functions as data.  In practice we never
actually observe a function in its entirety, thus
the data we work with in FDA are finite dimensional vectors,
and (superficially) have the same form as other types
of quantitative multivariate data.  But since the data are
considered to arise by evaluating smooth functions, we can
view the population of interest as consisting of functions,
even while the data are vectors. A distinct set of
methods have been developed to take advantage of this property.

## Data Depth

There are various ways to measure the *depth* of a point
$z \in {\mathbb R}^d$ relative to a distribution or collection of points
$\\{x_i\in {\mathbb R}^d; i=1,\ldots,n\\}$.  Formally, a *depth measure* on ${\mathbb R}^d$
is a function from ${\mathbb R}^d\rightarrow{\mathbb R}^+$ that quantifies the
depth of each point as a non-negative real number.

"Deep" points are surrounded in all
directions by many other points, while "shallow" points lie near the
surface or exterior of the point set.  Deep points are often referred to as
having high "centrality" while
shallow points have low centrality or high "outlyingness".
Data depth can be viewed as a multivariate generalization of the
notion of a quantile, with the deepest point in a set being a
type of multivariate median.

Below are several examples of depths.

### Halfspace depth

The original definition of depth was the *halfspace depth* introduced
by John Tukey in 1975.  The definition of the halfspace depth is
simple to describe graphically and a bit more difficult to define
formally.  To calculate the halfspace depth of a single point
$z\in {\cal R}^d$ with respect to a collection of points
$\{x_i; i=1, \ldots, n\}$, with each $x_i \in {\cal R}^d$, let $U$ denote the set
of all unit vectors in ${\cal R}^d$ and define the halfspace depth as

$$
D_H(z; \\{x_i\\}) = {\rm min}_{u\in U} n^{-1}\sum_{i=1}^n {\cal I}(u^\prime (x_i - z) > 0).
$$

$$
D_H(z; \\{x_i\\}) =
$$

$$
{\rm min}_{u\in U} n^{-1}\sum_{i=1}^n {\cal I}(u^\prime (x_i - z) > 0).
$$

What we are doing here is searching for a line passing through $z$
that places the greatest fraction of the $x_i$ on one side of the
line.  If $z$ falls at the geometric center of a collection of
symmetrically distributed points, then $z$ is as deep as possible
and will have halfspace depth approximately equal to 1/2.  At the
other extreme there is a line passing through $z$ such that all of
the $x_i$ are on the same side of this line.  In this case the point
$z$ is as shallow as possible and its halfspace depth is approximately
equal to zero.

The halfspace depth is geometrically natural but expensive to compute
exactly except in two dimensions.

### Spatial depth

The spatial depth has a simple definition that is relatively easy to
compute in high dimensions:

$$
D_S(z; \\{x_i\\}) = 1 - \\|{\rm Avg}_i\\{(x_i-z)/ \\|x_i-z\\|\\}\\|
$$

Note that $(x_i-z)/\\|x_i-z\\|$ is a unit vector pointing in the
direction from $z$ to $x_i$.  If a point $z$ is "shallow" then most of
the unit vectors $(x_i-z)/\\|x_i-z\\|$ point in roughly the same
direction, and therefore their average value will have large
magnitude.  If a point $z$ is "deep" then these unit vectors will
point in many different directions and their average value will have
small magnitude.

### L2 depth

The $L_2$ depth also has a simple definition and is easy to compute:

$$
D_{L_2}(z; \\{x_i\\}) = 1 / (1 + {\rm Avg}_i\{\\|x_i-z\\|\}).
$$

### Band depth

Another way to measure depth that is especially useful for vectors that are
sampled functions is the *band depth*.   In the most common version of the
band depth, we consider all triples of three distinct values in the reference
set $\\{x_i\\}$, say $x_{j_1}$, $x_{j_2}$, $x_{j_3}$ where $j_1 < j_2 < j_3$.  For
each such triple we calculate the proportion of indices $k$ such that $z(k)$
lies between the maximum and minimum of
$x_{j_1}(k)$, $x_{j_2}(k)$, $x_{j_3}(k)$.  These values are then averaged over all
triples of points in the reference set.  A more formal definition is as
follows:

$$
D_B(z; \\{x_i\\}) = {\rm Avg}_{j_1 < j_2 < j_3} {\rm Avg}_k {\cal I}[{\rm min}(x_{j_1}(k), x_{j_2}(k), x_{j_3}(k)) < z(k) < {\rm max}(x_{j_1}(k), x_{j_2}(k), x_{j_3}(k))].
$$

### Properties of a good depth function

Analysis based on depths does not directly rely on probability models,
making it quite distinct from many other methods of statistical data
analysis.  For statistical methods based on probability there are
standard properties such as unbiasedness, consistency, accuracy, and
efficiency that are used to evaluate the performance of the approach.
Although it is possible to place depth into a probabilistic framework so
that these notions can be applied, several researchers have attempted
to define the geometric properties that a depth function should
exhibit that do not depend on any probability framework.  Four basic
such properties are

* *Affine invariance* -- If the data are transformed by the affine
orthogonal mapping $x\longrightarrow c + Qx$, where $c\in {\cal R}^d$
is a fixed vector and $Q$ is an orthogonal matrix, then the depths do
not change.

* *Maximality at the center* -- If the data are symmetric around zero,
i.e. if $-x$ is in the dataset whenever $x$ is in the dataset, then
the vector $0_d$ achieves the maximum depth.

* *Monotonicity relative to the deepest point* -- Let $\tilde{x}$ be
the deepest point and we consider any unit vector $u$, and we then
evaluate the depth at $\tilde{x} + \lambda u$ for
$\lambda \in {\cal R}^+$, then the depth is a decreasing function of $\lambda$.

* *Vanishing at infinity* -- for any sequence $z_i$ with $\\|z_i\\|$
tending to infinity, the depths of the $z_i$ tend to zero.

### Depth peeling

Data depth can be used in exploratory multivariate analysis to identify
the most central or typical points and then contrast them with the more
outlying points.  A systematic way to do this is to stratify the data
based on depth and then inspect the points in each depth stratum.  For
example, if we stratify the data into 10 groups based on depth deciles,
the first decile consists of the shallowest 10% of points and the last
decile consists of the deepest 10% of points.

Often (not always) there is little heterogeneity in the deepest decile,
meaning that all of the deepest points are very similar.  However there
is nearly always heterogeneity in the shallowest decile, as there are
many different ways to be near the periphery of a collection of points.

## Quantization

A quantization algoithm aims to represent a multivariate distribution
through a small number of representative points.
This can be a useful exploratory technique if the distribution being
studied has a complex form that is not approximately Gaussian or elliptical,
and is not well captured through independent additive factors
(as in PCA).  The goal of almost any quantization algorithm is to
find a collection of representative points $\\{x_i\\}$ that are optimal
in some sense - for example we may wish to optimize the $x_i$ so as
to minimize the distance from any observation to its closest representative
point.  Inspecting the representative points may provide a quick
means to understand the structure of the distribution.

A recently developed algorithm constructs
[support points](https://arxiv.org/abs/1609.01811) that are a very
effective form of quantization.  To understand the support point
algorithm, suppose that we are given a distribution function $G$ that
we wish to approximate with a finite set of points.  The sample
space is ${\cal R}^d$, and let $Y$ denote a random draw
from $G$.  Now consider an approximating distribution $F$ with
random draw $X$.

We are given $G$ and wish to construct $F$
to approximate $G$, so we begin by defining a distance function
that measures how far apart $F$ and $G$ are from
each other.  Note that this distance compares two
probability distributions (it is not the familiar notion of a distance
between vectors).  Distances among probability distributions
play an important role in modern statistics.  When $d=1$ many natural
distance measures on probability distributions can be constructed, but it is harder
to construct good distances on probability distributions when the dimension
$d$ is greater than one.

One distance measure on probability distributions that turns out to be
very effective and relatively easy to work with is called the *energy distance* (it has other
names as well), and is defined as

$$
2E\\|X-Y\\| - E\\|X-X^\prime\\| - E\\|Y-Y^\prime\\|.
$$

Here, $X$, $X^\prime$ are independent draws from $F$, and $Y$, $Y^\prime$
are independent draws from $G$.  It turns out that expression above
is equal to zero if and only if $F \equiv G$.  This is a necessary property
for a distance to have.

The interpretation of the expression above is that $F$ is close to $G$
if (i) a random draw from $F$ tends to be close to a random draw from $G$,
(ii) two independent random draws from $F$ are far from each other, and
(iii) two independent random draws from $G$ are far from each other.

Our goal is to approximate a given distribution $G$ with a simpler distribution $F$
that we construct.  Here "simple" means that $F$ has finite support, i.e. a finite
sample space.  Since $G$ is given, the term $E\\|Y-Y^\prime\\|$ in the
energy distance is fixed
and can be ignored when constructing $F$.  Thus, our goal is to construct
$F$ that minimizes

$$
2E\\|X - Y\\| - E\\|X - X^\prime\\|.
$$

It is worth considering an alternative approach in which we simply minimize
the first term above, $E\\|X-Y\\|$.  However doing this always yields a degenerate
solution in which $F$ places all of its probability mass on the
*spatial median*, which is the fixed vector
$V$ that minimizes $E\\|Y - V\\|$.  This is the reason that the "repulsive"
term $E\\|X-X^\prime\\|$ in the distance measure is essential.

In practice, we do not observe the distribution $G$ but instead observe a sample
$y_1, \ldots, y_N$.  As stated above, the approximating distribution $F$ that we are constructing
is supported on a finite set of points $x_1, \ldots, x_n$.  This leads
us to the empirical analogue of the distance function above:

$$
\frac{2}{nN}\sum_{i=1}^n\sum_{j=1}^N\\|y_j - x_i\\| - \frac{1}{n^2}\sum_{i=1}^n\sum_{j=1}^n\\|x_i-x_j\\|.
$$

Our goal here was to discuss the motivation behind the support point algorithm.
We will not proceed further with discussion of the process of numerically minimizing this function (see
the paper linked above for computational details).

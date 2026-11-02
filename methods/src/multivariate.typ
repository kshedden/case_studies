#import "@preview/ilm:2.1.1": *

#show link: set text(fill: blue)

#set text(lang: "en")

#show: ilm.with(
  title: [Multivariate methods],
  authors: "Kerby Shedden",
  figure-index: (enabled: true),
  table-index: (enabled: true),
  listing-index: (enabled: true),
  chapter-pagebreak: false
)

= Introduction

This document discusses several modern methods for analyzing multivariate data. These methods were developed since around 1970, whereas the classical multivariate methods such as PCA, CCA, etc. were developed in the 1930's or earlier. Compared to these classical approaches, the methods discussed here make less use of factor analysis and dimension reduction using the singular value decomposition, and instead make more use of point set geometry and distances on spaces of probability measures. Methods based on covariance matrices work best for #link("https://en.wikipedia.org/wiki/Elliptical_distribution")[elliptically distributed] data, whereas the methods discussed here can in some cases overcome this limitation.

= Functional data

A particular type of multivariate data is known as _functional data_, in which we observe vectors $v in RR^d$ that arise from evaluating a function at a finite number of points, i.e. the $i^"th"$ observation is

$
v_i = [f_(i)(t)]_(t in cal(T)_i).
$

In practice we may only observe $v_i + epsilon_i$, where $epsilon_i$ is random measurement error.  For the functional perspective to be useful, we need to place some constraints on the functions $f_i$, usually involving some form of smoothness.

#link("https://en.wikipedia.org/wiki/Functional_data_analysis")[Functional Data Analysis (FDA)] encompasses many methods for analyzing this type of data. Arguably, this name is misleading since the data themselves are rarely if ever smooth functions in their entirety.  The data we work with in FDA are finite dimensional vectors, and (superficially) have the same form as other types of quantitative multivariate data. But since the data are considered to arise by evaluating smooth functions, we can view the population of interest as consisting of functions, even while the data are vectors.

= Data Depth

The #link("https://arxiv.org/abs/1207.4988")[depth] of a point $z in RR^d$ relative to a distribution or collection of points ${x_i in RR^d; i=1,...,n}$ is a model-free way to understand the geometry of a collection of vectors, viewed as points in a metric space. Formally, a _depth measure_ on $RR^d$ is a function from $RR^d -> RR^+$ that quantifies the depth of each point as a non-negative real number.

"Deep" points (with large depth values) are surrounded in all directions by many other points, whereas "shallow" or "superficial" points lie near the surface or exterior of the point set. Deep points are often referred to as having high "centrality" while shallow points have low centrality or high "outlyingness". Data depth can be viewed as a multivariate generalization of the notion of a quantile, with the deepest point of a distribution being a type of multivariate median.

Below are several specific examples of data depths.

== Halfspace depth

The original definition of depth was the
#link("https://en.wikipedia.org/wiki/Tukey_depth")[halfspace depth] introduced by John Tukey in 1975. The definition of the halfspace depth is simple to describe graphically and a bit more difficult to define formally. To calculate the halfspace depth of a single point $z in RR^d$ with respect to a collection of points ${x_i; i=1, ..., n}$, with each $x_i in RR^d$, let $U$ denote the set of all unit vectors in $RR^d$ and define the halfspace depth as

$
D_("HS")(z; {x_i}) = "min"_(u in U) sum_i cal(I)(u' (x_i - z) > 0).
$

What we are doing here is searching for a hyperplane passing through $z$ that places the greatest fraction of the $x_i$ on one side of the hyperplane ${v: u'(v-z) = 0}$. If $z$ falls at the geometric center of a collection of symmetrically distributed points, then $z$ is as deep as possible and will have halfspace depth approximately equal to 1/2. At the other extreme there is a hyperplane passing through $z$ such that all of the $x_i$ are on the same side of this hyperplane. In this case the point $z$ is as shallow as possible and its halfspace depth is equal to zero.

The halfspace depth is geometrically natural but expensive to compute exactly except in two dimensions.

== Spatial depth

The spatial depth has a simple definition that is relatively easy to compute in high dimensions:

$
D_(S)(z; {x_i}) = 1 - norm("Avg"_i (x_i-z) / norm(x_i-z)).
$

Note that $frac((x_i-z), norm(x_i-z), style: "horizontal")$ is a unit vector pointing in the direction from $z$ to $x_i$. If $z$ is "shallow" then most of the these unit vectors point in roughly the same direction, and therefore their average value will have large magnitude. If $z$ is "deep" then these unit vectors will point in many different directions and their average value will have small magnitude.

== L2 and Mahalanobis depth

The $L_2$ depth also has a simple definition and is easy to compute:

$
D_(L_2)(z; {x_i}) = 1 / (1 + "Avg"_i {norm(x_i-z)}).
$

The Mahalanobis depth incorporates information about the mean $mu$ and covariance $Sigma$ of the pointset ${x_i}$

$
D_(L_2)(z; {x_i}) = 1 / (1 + (z - mu)'Sigma^(-1)(z-mu)).
$

While fast to compute, these two depths have strong connections to the classical methods of multivariate analysis based on covariance matrices, and hence do not effectively achieve our goal of compensating for the limitations of these methods.

== Band depth

Another way to measure depth that is especially useful for vectors that are sampled functions is the _band depth_. In one common version of band depth, to calculate the depth of $z$ relative to the reference set ${x_i}$, we consider all triples of three distinct values in the reference set ${x_i}$, say $x_(j_1)$, $x_(j_2)$, $x_(j_3)$ where $j_1 < j_2 < j_3$. For each such triple we calculate the proportion of indices $k$ such that $z(k)$ lies between the maximum and minimum of $x_(j_1)(k)$, $x_(j_2)(k)$, $x_(j_3)(k)$. These values are then averaged over all triples of points in the reference set.

== Properties of a good depth function

Analysis based on depths does not directly rely on probability models, making it quite distinct from many other methods of statistical data analysis. For statistical methods based on probability there are standard properties such as unbiasedness, consistency, accuracy, and efficiency that are used to evaluate the performance of the approach. Although it is possible to place depth into a probabilistic framework so that these notions can be applied, several researchers have attempted to define the geometric properties that a depth function should exhibit that do not depend on any probability framework. Four
basic such properties are

- _Affine invariance_ -- If the data are transformed by the affine orthogonal mapping $x -> c + Q x$, where $c in RR^d$ is a fixed vector and $Q$ is an orthogonal matrix, then the depths do not change.

- _Maximality at the center_ -- If the data are symmetric around zero, i.e. if $-x$ is in the dataset whenever $x$ is in the dataset, then the vector $0_d$ achieves the maximum depth.

- _Monotonicity relative to the deepest point_ -- Let $tilde(x)$ be the deepest point and consider any unit vector $u$.  If we evaluate the depth at $tilde(x) + lambda u$ for $lambda in RR^+$, then the depth is a decreasing function of $lambda$.

- _Vanishing at infinity_ -- for any sequence $z_i$ with $norm(z_i)$ tending to infinity, the depths of the $z_i$ tend to zero.

== Depth peeling

Data depth can be used in exploratory multivariate analysis to identify the most central or typical points and then contrast them with the more outlying points. A systematic way to do this is to stratify the data based on depth and then inspect the points in each depth stratum. For example, if we stratify the data into 10 groups based on depth deciles, the first decile consists of the shallowest 10% of points and the last decile consists of the deepest 10% of
points.

Often (not always) there is little heterogeneity in the deepest decile, meaning that all of the deepest points are very similar. However there is nearly always heterogeneity in the shallowest decile, as there are many different ways to be near the periphery of a collection of points.

= Metrics on spaces of probability distributions

Many statistical tasks can be defined in terms of comparisons between probability distributions. For example, the "distribution version" of the two-sample problem considers distributions $F_X$ and $F_Y$ corresponding to two populations of interest.  A distance $d(F_X, F_Y)$ can be used as a descriptive statistic to describe the relationship between the two distributions, and we can perform formal hypothesis tests for the null hypothesis $F_X = F_Y$, which is equivalent to $d(F_X, F_Y) = 0$.  A more advanced setting is where we have numerous populations $F_1, ..., F_m$, and we wish to quantify the overall discrepancy among the $F_i$, or test the global null hypothesis that all $F_i$ are equal.

A powerful metric on the space of probability distributions is the #link("https://en.wikipedia.org/wiki/Total_variation_distance_of_probability_measures")[total variation] (TV) metric.  The TV distance between two probability measures $P_A$ and $P_B$ is

$
"sup"_(S in cal(A)) |P_(A)(S) - P_(B)(S)|,
$

where $cal(A)$ consists of all measurable sets.  If $P_A$ and $P_B$ have densities $f_A$ and $f_B$, then the TV distance becomes $integral |frac(f_(A)(x) - f_(B)(x)|, 2, style: "horizontal") "dx"$.  Closely related is the #link("https://en.wikipedia.org/wiki/Hellinger_distance")[Hellinger distance], which is $integral (sqrt(f_(A)(x)) - sqrt(f_(B)(x)))^2 "dx"$.  These metrics are useful but often difficult to work with in practice since estimation usually requires estimation of densities.

Another class of metrics on spaces of probability measures are the Wasserstein distances.  These are closely related to the idea of #link("https://en.wikipedia.org/wiki/Transportation_theory_(mathematics)")[optimal transport], which we discuss below.

A final but very practical class of metrics on spaces of probability distributions is the #link("https://en.wikipedia.org/wiki/Energy_distance")[energy distance].  Given two random vectors $X$ and $Y$, the energy distance between their distributions is

$
2E norm(X-Y) - E norm(X-X') - E norm(Y-Y'),
$

where $X'$ is an independent copy of $X$ and $Y'$ is an independent copy of $Y$.  This idea can be applied whenever the sample spaces for $X$ and $Y$ lie in a common metric space with metric $d$.  We can define the energy distance in this setting as

$
2E d(X,Y) - E d(X, X') - E d(Y, Y').
$

An interpretation of the energy distance is that $F$ is close to $G$ if (i) a random draw from $F$ tends to be close to a random draw from $G$, (ii) two independent random draws from $F$ are far from each other, and (iii) two independent random draws from $G$ are far from each other.  The first of these terms is "attractive" while the latter are "repulsive".  The energy distance is a balance between these attractive and repulsive terms.

There is a straightforward plug-in estimate of the energy distance, which is simply the energy distance between the empirical distributions $n_x^(-1)sum_i bb(1)_(x_i)$ and $n_y^(-1)sum_i bb(1)_(y_i)$.

$
2/(n_x n_y) sum_(i,j) norm(x_i - y_j) - 1/(n_(x)^2) sum_(i,i')norm(x_i - x_(i')) - 1/(n_(y)^2) sum_(i,i')norm(y_i - y_(i')).
$


This is a #link("https://en.wikipedia.org/wiki/V-statistic")[V-statistic], and there is a closely related exactly unbiased estimator that is a #link("https://en.wikipedia.org/wiki/U-statistic")[U-statistic]:

$
2/(n_x n_y) sum_(i,j) norm(x_i - y_j) - 1/(n_(x)(n_x-1)) sum_(i != i')norm(x_i - x_(i')) - 1/(n_(y)(n_y-1)) sum_(i != i')norm(y_i - y_(i')).
$

An advantage of the energy distance over the total variation and Hellinger metrics is that it does not involve densities.  However direct estimation of the energy distance can be expensive as it involves calculating $n^2$ vector norms, where $n$ is the sample size.

= Quantization

A quantization algorithm aims to represent a multivariate distribution through a small number of representative points. This can be a useful exploratory technique if the distribution being studied has a complex form that is not approximately Gaussian or elliptical, and is not well captured through independent additive factors (as in PCA). The goal of almost any quantization algorithm is to find a collection of representative points ${x_i}$ corresponding to the distribution of a random variable $Y$ that are optimal in some sense. Inspecting the representative points may provide a quick means to understand the structure of the distribution.

#link("https://arxiv.org/abs/1609.01811")[Support points] are an effective form of quantization that is based on the energy distance. Suppose we are given a distribution function $G$ on $RR^d$ that we wish to quantize, i.e. to approximate with a finite set of points $x_1, ..., x_m$. Let $Y$ denote a random draw from $G$, and let $X$ denote a draw from $m^(-1) sum_j bb(1)_(x_j)$.  The support points are the values of $x_1, ..., x_m$ that minimize the energy distance between $Y$ and $X$.

Construction of support points involves solving a non-convex optimization problem.  There is a natural majorization/maximization approach, and faster gradient approaches are also used. Due to nonconvexity of the loss function, the support points are not uniquely defined, and it is not generally possible to find more than a local minimizer of the loss. The non-uniqueness is evident in particular when the distribution of $X$ has symmetries.  For example, if the distribution is invariant to rotation, any set of points can be rotated without changing its energy distance to $Y$.

= Distance correlation

A classical approach for quantifying the degree of association between two jointly distributed random vectors is canonical correlation analysis (CCA).  #link("https://en.wikipedia.org/wiki/Distance_correlation")[Distance correlation] addresses this topic in a different way, using ideas that are related to the energy distance.  The goal of distance correlation is to construct a way of measuring the association between two random vectors that is zero if and only if the random vectors are independent.  This can be contrasted with the finding of all canonical correlations being zero, which only implies that $X$ and $Y$ are uncorrelated.

There are several equivalent definitions of the distance covariance, from which the distance correlation is obtained.  Two direct definitions are

$
"dCov"^2(X, Y) = "cov"(norm(X-X'), norm(Y-Y')) - 2"cov"(norm(X-X'), norm(Y-Y''))
$

and

$
"dCov"^2(X, Y) = E[norm(X-X') norm(Y-Y')] + E[norm(X-X')]E[norm(Y-Y')] - \ 2E[norm(X-X')norm(Y-Y'')].
$

In these expressions $(X, Y)$ and $(X', Y')$ are drawn from the joint distribution of $X$ and $Y$, with $(X, Y)$ being independent of $(X', Y')$, and $Y''$ being yet another independent copy of $Y$.

There is yet another characterization of distance covariance defined in terms of characteristic functions.  This definition makes it very clear that the distance covariance is zero if and only if $X$ and $Y$ are independent, but lacks the geometric intuition of the definitions given above.

In direct analogy to covariance and Pearson correlation, the distance variance of $X$ is the distance covariance of $X$ with itself, and the distance correlation between $X$ and $Y$ is their distance covariance divided by the product of their distance standard deviations (the square roots of their distance variances).

= Optimal transport and Wasserstein distance

Given two probability $P_X$ and $P_Y$ with sample spaces $cal(X)$ and $cal(Y)$, respectively, there are infinitely many ways to construct a joint distribution $P_(X,Y)$ on $cal(X) times cal(Y)$ such that $P_X$ and $P_Y$ are the marginal distributions of $P_(X,Y)$.  That is, $P_(X)(A) = P_(X,Y)(A, cal(Y))$ and $P_(Y)(B) = P_(X,Y)(cal(X),B)$. Any such joint distribution is called a _coupling_ of $P_X$ and $P_Y$.

Given a non-negative cost function $c(dot, dot)$, we can seek to minimize the loss function

$
integral c(x, y) "dP"_(X,Y)(x,y)
$

over the space of couplings $P_(X,Y)$.  This is known as _optimal transport_.

There are many applications of this idea in applied statistics that involve interpreting the structure of the optimal coupling.  For example, counterfactual "synthetic controls" can be used in causal analysis.  If ${x_i}$ are a collection of control observations and ${y_i}$ are a collection of treated observations, then we can couple $n_x^(-1)sum_i bb(1)_(x_i)$ with $n_y^(-1)sum_i bb(1)_(y_i)$.  The resulting transport map can be used to infer the matched treated value to each $x_i$ and the matched control value to each $y_i$.

The Wasserstein distance between probability distributions on $cal(S)$ is the solution to an optimal transport problem in which the cost function is a metric on $cal(S)$.  Specifically, the $p^"th"$ Wasserstein distance is

$
"argmin"_(P_(X,Y)) [integral d(x, y)^p "dP"_(X,Y)(x,y)]^(1/p).
$


= Minimum covariance determinant

Much of classical multivariate analysis is based on the
#link("https://en.wikipedia.org/wiki/Covariance_matrix")[covariance matrix] of several variables. The most common way to estimate this matrix is the #link("https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices")[sample covariance matrix],
which (but for a minor scaling difference) coincides with the maximum likelihood estimate for Gaussian data (note however that the data do not need to follow a Gaussian distribution for the sample covariance matrix to be
informative). The sample covariance matrix has some desirable properties including being unbiased and elementwise consistent. However for larger dimensions or smaller sample sizes, the sample covariance matrix can perform poorly. Also, the sample covariance matrix is not robust to outliers.

To address these limitations, many alternative estimators of the covariance matrix have been devised. Here we describe one of them, the #link("https://arxiv.org/abs/1709.07045")[minimum covariance determinant], or MCD. We will not give full details here, but the main idea is to select a constant $m >= n/2$, where $n$ is the sample size, and define the MCD estimate to be the sample covariance matrix of $m$ observations that has the smallest determinant among all such subsets.

The determinant of a covariance matrix is a measure of its dispersion. If $x$ follows the multivariate normal distribution $N_p(mu, Sigma)$, the squared Mahalanobis distance $(x-mu)' Sigma^(-1)(x-mu)$ follows the $chi^2_p$ distribution, where $x in RR^p$. It follows that the ellipsoid

$
{x: (x-mu)' Sigma^(-1)(x-mu) <= t}
$

covers $F_(p)(t)$ probability mass of the $N_p(mu, Sigma)$ distribution, where $F_p$ is the CDF of the $chi^2_p$ distribution. Setting $t = F_p^(-1)(0.95)$, for example, gives an ellipsoid that covers 95% of the probability. The volume of this ellipsoid is $t^p dot.c |Sigma|^(1/2)$ times the volume of the unit ball, which is $frac(pi^(p/2), Gamma(frac(p, 2, style: "horizontal")+1), style: "horizontal")$.

The rationale for the MCD estimator is that if a few observations greatly increase the determinant, they may be distorting the covariance estimate. As with any "robust" estimator, there is an implicit assertion that certain points may be outliers and should be discounted or excluded from the analysis. However there is no guarantee that any particular definition of an outlier is correct. In the case of the MCD, we are asserting that points that inflate the volume of an estimated _high probability ellipsoid_ may be outliers.

Exact calculation of the MCD estimator of the covariance matrix is challenging. Enumerating all possible subsets of $m$ out of $n$ observations is impossible in practice. Therefore a greedy approach called _fast MCD_ is used in practice.

One use for the covariance matrix is to quantify the "outlyingness" of individual observations. Specifically, let $x in RR^d$ denote an observation, $mu in RR^d$ denote the mean, and $Sigma in RR^(d times d)$ denote the covariance matrix. The squared #link("https://en.wikipedia.org/wiki/Mahalanobis_distance")[Mahalanobis distance] is defined to be $(x-mu)' Sigma^(-1)(x-mu)$. Points with larger Mahalanobis distance to the center are "less central" or "more outlying". Using the MCD or another robust estimate of the covariance matrix (and of the mean) can sometimes reveal more interesting sets of outliers than those that are identified by the sample covariance matrix.

There is a connection between the MCD and data depth using the "Mahalanobis depth" (MD), which is defined to be $frac(1, (1 + (x-mu)' Sigma^(-1)(x-mu)), style: "horizontal")$. The conventional MD uses the sample covariance to estimate $Sigma$, but this opens the possibility of "masking", where a small group of outliers distorts the covariance matrix so as to mask their outlyingness. Using the MCD to estimate $Sigma$ circumvents this issue.

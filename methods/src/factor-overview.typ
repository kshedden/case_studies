#import "@preview/ilm:2.1.1": *

#show link: set text(fill: blue)

#set text(lang: "en")

#show: ilm.with(
  title: [Overview of factor analysis and dimension reduction],
  authors: "Kerby Shedden",
  figure-index: (enabled: true),
  table-index: (enabled: true),
  listing-index: (enabled: true),
  chapter-pagebreak: false
)

#let notindep = $#h(0.5em) cancel(perp #h(-6pt) perp) #h(0.5em)$

#let indep = $#h(0.5em) perp #h(-6pt) perp #h(0.5em)$


= Dimension reduction, factor analysis, and embeddings

A large class of powerful statistical methods considers data in which many "objects" ("observations") are measured with respect to multiple variables. These analyses often focus on understanding both the relationships among the variables and the relationships among the observations. The data take the form of a rectangular array, where by convention the rows correspond to objects and the columns correspond to variables.

For example, we may have a #link("https://en.wikipedia.org/wiki/Sampling_(statistics)")[sample] of people (the "objects") with each person measured in terms of their height, weight, age, and sex (the "variables"). In this example, two people are similar if they have similar values on most or all of the variables, and two variables are similar if knowing the value of one of the variables for a particular object allows one to predict the values of other variables for the same object.

The branch of statistics known as #link("https://en.wikipedia.org/wiki/Multivariate_statistics")[multivariate analysis] is concerned with analysis of such datasets. Many methods of classical multivariate statistics operate specifically on rectangular arrays of data, as described above. Modern multivariate statistics provides methods that can operate on much more diverse forms of multivariate data, such as partially observed arrays, infinite dimensional arrays, and
#link("https://en.wikipedia.org/wiki/Tensor")[tensors].

Matrix factorization is a powerful idea from linear algebra, and underlies many of the most important methods in classical multivariate statistics. We will use the term "factor analysis" here to refer in a broad way to any statistical method relying on a linear algebraic factorization of the data matrix (or a transformed version of the data matrix).

To visualize such a data matrix, we can think in terms of _variable space_ and _case space_ (or _object space_). If there are $n$ cases and $p$ variables, then the cases are points in $RR^p$ and the variables are points in $RR^n$.

= Embedding

Many methods of factor analysis can be viewed as a way to obtain an #link("https://en.wikipedia.org/wiki/Embedding")[embedding]. Embedding algorithms take input data vectors $x$ and transform them into output data vectors $z$. Many embedding algorithms take the form of a #link("https://en.wikipedia.org/wiki/Dimensionality_reduction")[dimension reduction], so that $q equiv "dim"(z) < "dim"(x) equiv p$. However some embeddings preserve or even increase the dimension. The latter is related to the notion of #link("https://en.wikipedia.org/wiki/Overcompleteness")[overcompleteness] in linear algebra.

Embeddings can be used for exploratory analysis, especially in visualizations, and can also be used to construct features for prediction (e.g. as in #link("https://en.wikipedia.org/wiki/Principal_component_regression")[Principal Component Regression (PCR)]), as well as being used in formal statistical inference such as hypothesis testing.

An embedding approach is linear if $z = B'x$ for a fixed $p times q$ matrix $B$. Linear embedding embedding algorithms are simpler to devise and characterize than nonlinear embeddings. Many modern embedding algorithms are
#link("https://en.wikipedia.org/wiki/Nonlinear_dimensionality_reduction")[nonlinear], and use this additional flexibility to better capture complex structure.

If $B$ is completely independent of any data then the embedding is _non-adaptive_. For example, the #link("https://en.wikipedia.org/wiki/Discrete_Fourier_transform")[discrete Fourier transfrom (DFT)] projects the data onto a fixed set of orthogonal trigonometric basis vectors. #link("https://en.wikipedia.org/wiki/Wavelet_transform")[Wavelet transforms] can also be seen as a non-adaptive linear dimension reduction in this sense.

If the transform matrix $B$ is based on training data, the embedding is _adaptive_. The purpose of this adaptation is to do a better job capturing the structure in the values that are most likely to be encountered.

Some embedding algorithms embed only the objects while other embedding algorithms embed both the objects and the variables. Embedding the objects provides a reduced feature representation for each object that can be passed on to additional analysis procedures, or interpreted directly. Embedding the variables provides a means to interpret the relationships among the variables. When utilizing an embedding algorithm that embeds both objects and variables we have the opportunity to interpret the results through a #link("https://en.wikipedia.org/wiki/Biplot")[biplot].

Finally, it is important to be aware that some methods of multivariate analysis do not take the form of a dimension reduction at all. A prominent example of this type of approach is
#link("https://en.wikipedia.org/wiki/Cluster_analysis")[clustering], which groups observations into discrete categories such that the observations within a category are similar. Whereas dimension reduction places the observations on a multi-dimensional continuum or spectrum, which is itself situated in a metric space, cluster analysis seeks discrete groups of observations that are only interpretable as unordered categories. Put geometrically, cluster analysis seeks disconnected "islands" of observations, while dimension reduction aims to organize the observations onto a low dimensional continuum. In either case, the goal is to retain as much information as possible while simplifying the representation of the data.

== Centering and standardization

Suppose that $X$ is a $n times p$ matrix whose rows are the observations or objects, and whose columns are the variables. There are many ways to pre-process the data in $X$ prior to performing a factor analysis. The purpose of this preprocessing is to suppress the aspects of the data that are not of primary interest, such as the
#link("https://en.wikipedia.org/wiki/Location_parameter")[locations] and #link("https://en.wikipedia.org/wiki/Scale_parameter")[scales] of the variables. A dimension reduction method that yields the same result upon linear transformation of the variables is called
#link("https://en.wikipedia.org/wiki/Eqiuvariant_map")[equivariant], and this is usually considered to be a desirable property of a factor analysis procedure.

Some factor-type methods depend only on the sample covariance matrix of the
variables, which is a $p times p$
(#link("https://en.wikipedia.org/wiki/Definite_matrix")[positive semidefinite, PSD)] matrix. For these methods, the sample covariance acts as a
#link("https://en.wikipedia.org/wiki/Sufficient_statistic")[sufficient statistic]. Since covariances by definition are derived from mean centered variables, it is common to center the variables (columns) of the data matrix $X$ prior to conducting a factor analysis. Centering usually is accomplished by subtracting the mean from each column, but in some cases the columns may be centered with respect to another measure of location such as the median. After this centering, the case vectors (rows of $X$) become a point cloud that is centered around the origin in $RR^p$, and all subsequent analysis can be interpreted as aiming to understand the pattern of deviations from the mean, or
#link("https://en.wikipedia.org/wiki/Centroid")[centroid] of the original data.

We may wish for our results to be independent of the units in which each variable was measured, and remove any influence of differing dispersions of the variables on our results. This motivates _standardizing_ the columns of $X$ prior to performing a factor analysis. The most common way to standardize data is to divide each column by a #link("https://en.wikipedia.org/wiki/Statistical_dispersion")[measure of scale] such as the standard deviation or inter-quartile range. Usually standardization is done after centering, but in some cases (e.g. for a variable that is naturally constrained to positive values), no centering is performed.

In some datasets there is strong heterogeneity among the observations (rows) and it is desirable to supress this in the analysis. In this case we may choose to _double center_ the data so that every row and every column has mean zero. To achieve this, the following three steps can be performed: (i) center the overall matrix around its grand mean, (ii) center each row of the matrix with respect to the row mean, (iii) center each column of the matrix with respect to the column mean. The resulting matrix will have centered rows and centered columns, and it can be shown that the same result is obtained if steps (ii) and (iii) above are performed in either order.

= Singular Value Decomposition

Many embedding methods make use of a matrix factorization known as the
#link("https://en.wikipedia.org/wiki/Singular_value_decomposition")[Singular Value Decomposition], or SVD. The SVD is defined for any $n times p$ matrix $X$. In most cases we want $n >= p$, and if $n < p$ we can take the SVD of $X'$. When $n >= p$, we decompose $X$ as $X = U S V'$, where $U$ is $n times p$, $S$ is $p times p$, and $V$ is $p times p$. The matrices $U$ and $V$ are orthogonal so that $U' U = V' V = I_p$, and $S$ is diagonal with $S_(1 1) >= S_(2 2) >= dots.h.c >= S_(p p)$. The values on the diagonal of $S$ are the _singular values_ of $S$, and the SVD is unique except when there are ties among the singular values.

One use of the SVD is to obtain a low rank approximation to a matrix $X$. The extreme example of a low rank matrix is a matrix with rank 1, which means that we can write $X_(i j) = r_i c_j$ for vectors $r in RR^n$ and $c in RR^p$, or equivalently $X = r c'$. Our goal will be to decompose $X$ as the (approximate) sum of $k$ rank-1 matrices, and these matrices are ordered by their importance in describing the variation in $X$.

Suppose we truncate the SVD using only the first $k$ components, so that $tilde(U)$ is the $n times k$ matrix consisting of the leading (left-most) $k$ columns of $U$, $tilde(S) = "diag"(S_1, ..., S_k)$, and $V$ is the $p times k$ matrix consisting of the leading $k$ columns of $V$. In this case, the matrix $tilde(X) equiv tilde(U) tilde(S) tilde(V)^T$ is a rank $k$ matrix (it has $k$ non-zero singular values). According to the #link("https://en.wikipedia.org/wiki/Low-rank_approximation")[Eckart-Young] theorem, among all rank $k$ matrices, $tilde(X)$ is the closest matrix to $X$ in the _Frobenius norm_, which is defined as

$
"Frob"(X)^2 = norm(X)_F^2 equiv sum_(i j) X_(i j)^2 = "trace"(X' X).
$

Thus we have

$
tilde(X) = "argmin"_(A: "rank"(A) = k) norm(A - X)_F.
$

Finding a low rank approximation is a type of least square problem, although it is not the linear least squares encountered in regression analysis.

Thinking geometrically, the rows of $tilde(X)$ are spanned by a basis of $q$ elements, and the span of this basis captures most of the variation in the rows of $X$. Also, the columns of $tilde(X)$ are spanned by a basis of $q$ elements, and the span of this basis captures most of the variation in the columns of $X$.

== Analyzing a data matrix using the SVD

Suppose we have a data matrix $X$ and wish to understand its structure. We can begin by double centering the matrix $X$ as discussed above:

$
X_(i j) = m + r_i + c_j + R_(i j).
$

Here, $m in RR$ is the grand mean of $X$, $r in RR^n$ contains the row means of $X-m$, $c in RR^p$ contains the column means of $X-m$, and $R_(i j) equiv X_(i j) - m - r_i - c_j$ are residuals. Next we can take the SVD of $R$, yielding $R = U S V'$, or equivalently

$
R_(i j) = sum_(k=1)^p S_(k k) U_(i k)V_(j k).
$

Although there are $p$ terms in the SVD of $R$, the first few terms may capture most of the structure, so

$
R_(i j) approx sum_(k=1)^q S_(k k) U_(i k)V_(j k)
$

for $q < p$ (this approximation holds better when the _tail singular values_ $S_(q+1,q+1), ..., S_(p,p)$ are small). One important property that results from calculating the SVD for a double-centered matrix is that $U_(dot.c k) = 0$ and $V_(dot.c k) = 0$. That is, the columns of $U$ and $V$ are centered. This column centering means that the SVD captures "deviations from the mean" represented by the additive model $m + r_i + c_j$.

=== Assessing dimensionality

The rate at which the singular values decay contains important information about the intrinsic dimensionality of the population under study. Let $s_j$ denote the $j^"th"$ singular value, where $s_1 >= s_2 dots.h.c$. Two canonical patterns of decay that may be found are an exponential decay $s_j approx a exp(-b dot.c j)$ and a power-law decay $s_j approx frac(a, j^b, style: "horizontal")$. To assess the decay of the eigenvalues, we can consider plots of $log(s_j)$ against $j$, which is linear in the case of exponential decay, and plots of $log(s_j)$ against $log(j)$, which is linear in the case of power-law decay.

== Spectral decomposition of a covariance matrix

Given a collection of $p$ jointly distributed random variables $X_1, ..., X_p$ with sample space $RR$ (i.e. the random variables are _quantitative_ and _real-valued_), the covariance matrix $Sigma$ is a $p times p$ matrix such that $Sigma_(i j) = "Cov"(X_i, X_j)$.

The covariance matrix is symmetric, since $"Cov"(X_i, X_j) = "Cov"(X_j, X_i)$. The diagonal elements of a covariance matrix are the variances of the variables, since $"Cov"(X_i, X_i) = "Var"(X_i)$. The covariance matrix is positive semi-definite (PSD), which means that for any vector $v in RR^p$, $v' Sigma v >= 0$. This inequality is expected since $v' Sigma v = "Var"(v' x)$, where $X = (X_1, \ldots, X_p)'$, but the PSD property originates in the study of _quadratic forms_ (multivariate quadratic functions), so a PSD matrix defines a quadratic form that cannot take on negative values.

The existence of _eigenvalues_ and a complementary orthogonal basis of _eigenvectors_ is a complicated question in general. But symmetric matrices always have a basis of orthogonal eigenvectors and corresponding real eigenvalues. That is, we can write any $p times p$ symmetric matrix $S$ in the form $S = V Lambda V'$ where $V$ is $p times p$ and orthogonal ($V'V = I$), and $Lambda$ is a diagonal matrix. If in addition we know that $S$ is PSD, then the diagonal elements of $Lambda$ are non-negative. If the diagonal elements of $Lambda$ are strictly positive, then $S$ is _positive definite_ (PD). The factorization $S = V Lambda V'$ is known as the #link("https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix")[spectral decomposition].

In practice, we take the spectral decomposition of the sample covariance matrix $hat(Sigma)$. If $tilde(X)$ is the column-centered data matrix, then $tilde(X)' frac(tilde(X), n, style: "horizontal")$ is the sample covariance matrix (dividing by $n-1$ yields an unbiased estimate of the population covariance matrix but this only matters for small samples). There is a close connection between the SVD of $tilde(X)$ and the spectral decomposition of the sample covariance matrix. Let $tilde(X) = U S V'$ denote the SVD of the column-centered data matrix, so that the sample covariance matrix is $hat(Sigma) = tilde(X)' frac(tilde(X), n, style: "horizontal") = V S^2 frac(V', n, style: "horizontal")$. Then if $hat(Sigma) = V Lambda V'$ is the spectral decomposition of $hat(Sigma)$, we see that the $V$ terms in these two decompositions are the same, and
$Lambda = S^2$.

A small caveat to this discussion is when there are tied eigenvalues (i.e. several elements of $S$ are identical). In this case, neither the SVD nor the spectral decomposition are unique. The sample covariance matrix will not exhibit such ties (except in corner cases where the data are perfectly degenerate). However it is difficult to use the sample covariance matrix to assess whether the population covariance matrix has tied eigenvalues. Moreover, when the eigenvalues of $Sigma$ are very similar, even if they are not identical, we will still experience the consequences of the near non-uniqueness of the spectral decomposition. A special case of tied eigenvalues is when all of the eigenvalues are equal, a condition known as _sphericity_.

= Dimension reduction of grouped data

Suppose we have multivariate grouped data, such as vectors of individual characteristics for people who have been partitioned into $G$ groups, The data for each individual is a $p$-dimensional vector of characteristics, so we let $y_(i j) in RR^p$ denote the $j^"th"$ observation in the $i^"th"$ group. There are $n_i$ observations in the $i^"th"$ group, so $j=1,... ,n_i$. Let $m_i in RR^p$ denote the mean (centroid) of group $i$, and let $r_(i j) = y_(i j) - m_i in RR^p$ denote the residual vector for observation $y_(i j)$. The overall _scatter matrix_ is $E = sum_(i j) r_(i j)r'_(i j) in RR^(p times p)$. Letting $n=sum_i n_i$, one could view $frac(E, n, style: "horizontal")$ as the _residual covariance matrix_. Let $m$ be the centroid of the entire collection of observations, and define the _hypothesis matrix_ as $H = sum_i n_i (m_i - m)(m_i - m)' in RR^p$. We can interpret $frac(H, n, style: "horizontal")$ as the "between group covariance" matrix.

In classical multivariate analysis, a technique known as
#link("https://en.wikipedia.org/wiki/Multivariate_analysis_of_variance")[multivariate analysis of variance] (MANOVA) allows us to understand the within and between group structure captured by the matrices $H$ and $E$, and to reduce the dimension of the data while preserving as much of the between-group information as possible. A closely related
technique is
#link("https://en.wikipedia.org/wiki/Linear_discriminant_analysis")[linear discriminant analysis (LDA)].

The classical approaches to LDA and MANOVA are all based on the eigenvalues and eigenvectors of $E^(-1)H$, which are also the generalized eigenvalues of $H$ relative to $E$. Let $lambda_1 >= lambda_2 >= dots.h.c$ represent these eigenvalues in descending order, and let $eta_1, eta_2, ...$ denote the corresponding eigenvectors. The natural interpretation of $lambda_k$ is that it represents the ratio of between-group to within-group variance of the data projected in the direction $eta_k$, i.e. of the scalars $eta_k^' y_(i j)$. The transformed eigenvalue statistic $frac(lambda_k, (1 + lambda_k), style: "horizontal")$ is the proportion of total variance of the projected data that is between-group, which is essentially an $R^2$ statistic for the projected data.

= Principal Components Analysis

Suppose that $x$ is a $p$-dimensional random vector with mean $0$ and
covariance matrix $Sigma$ (our focus here is not the mean, so if $x$ does not have mean zero we can replace it with $x - mu$, where $mu = E[x]$). Principal Components Analysis (PCA) seeks a linear embedding of $x$ into a lower dimensional space of dimension $q < p$. The standard approach to PCA gives us a $p times q$ orthogonal matrix $B$ of _loadings_, which can be used to produce _scores_ $Q = B' x$, which are $q$-dimensional vectors.

PCA can be viewed in terms of linear compression and decompression via dimension reduction of the distribution of $x$. Let

$
x -> B'_(:,1:q) x equiv Q(x)
$

denote the compression that reduces the data from $p$ dimensions to $q$ dimensions, where $B_(:,1:q)$ is the $p times q$ matrix consisting of the leading $q$ columns of $B$. We can now decompress the data as follows:

$
Q(x) -> B_(:,1:q)Q(x) equiv hat(x).
$

This represents a two-step process of first reducing the dimension, then predicting the original data using the reduced data. Among all possible loading matrices $B$, the PCA loading matrix loses the least information in that it minimizes the expected value of $norm(x - hat(x))^2$.

The loading matrix $B$ used in PCA is a truncated eigenvector matrix of $Sigma = "Cov"(x)$. Specifically, we can write $Sigma = B Lambda B'$, where $B$ is an orthogonal matrix and $Lambda$ is a diagonal matrix with $Lambda_(1 1) >= Lambda_(2 2) >= dots.h.c >= Lambda_(p p) > 0$. This is the spectral decomposition of $Sigma$. The columns of $B$ are orthogonal in both the Euclidean metric and in the metric of $Sigma$, that is, $B'B = I_p$ and $B' Sigma B = I_p$. As a result, the scores $Q equiv B'x$ are uncorrelated, $"Cov"(Q) = I_q$.

Next we consider how PCA can be carried out with a sample of data, rather than in a population. Given a $n times p$ matrix of data $X$ whose rows are iid copies of the random vector $x$, we can estimate the covariance matrix $Sigma$ using the method of moments.  First, column center $X$ to produce $tilde(X) equiv X - 1_n overline(x)'$, where $overline(x) in RR^p$ is the vector of column-wise means of $X$, Then set $hat(Sigma) =  tilde(X)'frac(tilde(X), n, style: "horizontal")$. Letting $B$ denote the eigenvectors of $hat(Sigma)$, the scores have the form $Q = tilde(X)B$.

Since the eigenvalues $Lambda_(i i)$ are non-increasing, the leading columns of $Q$ contain the greatest fraction of information about the distribution of $x$. Thus, visualizations (e.g. scatterplots) of the first two columns of $Q$ best reflect the relationships in $X$ (compared to any other scatterplot formed from linear scores).

PCA can also be carried out using the SVD of the column-centered data matrix $tilde(X)$. Write $tilde(X) = U S V'$ as discussed above and note that

$
hat(Sigma) = tilde(X)' frac(tilde(X), n, style: "horizontal") = V S^2 frac(V', n, style: "horizontal").
$

Thus, $hat(Sigma)V = V frac(S^2, n, style: "horizontal")$, so $V$ contains the eigenvectors of $hat(Sigma)$ and the corresponding eigenvalues are in $frac("diag"(S^2), n, style: "horizontal")$.

== Karhunen-Loeve decomposition

PCA can be seen as a way to estimate the terms of the
#link("https://en.wikipedia.org/wiki/Kosambi%E2%80%93Karhunen%E2%80%93Lo%C3%A8ve_theorem")[Karhunen-Loeve (K-L) decomposition]
of a multivariate probability distribution. If $Y in RR^p$ is a random vector, then the K-L decomposition expresses $Y$ in the form

$
Y = mu + sum_(j=1)^p eta_j V_j,
$

where $mu = E[Y] in RR^p$ is a fixed vector containing the element-wise mean of $Y$, the $eta_j$ are a collection of mutually uncorrelated scalar random variables, and the $V_j$ are a collection of mutually orthogonal fixed vectors of dimension $p$. The expected values of the $eta_j$ are all zero and their variances are non-increasing in $j$. That is, $E[eta_j] = 0$, $"Var"(eta_1) >= "Var"(eta_2) >= dots.h.c$, and $"Cor"(eta_j, eta_k) = 0$ if $j != k$. If $Sigma$ is the $p times p$ covariance matrix of $Y$ and $Sigma = V Lambda V'$ is the spectral decomposition of $Sigma$, then the variance of $eta_j$ in the KL decomposition is $Lambda_(j j)$.

Each component $eta_j V_j$ represents random variation in the direction of the unit vector $V_j$. The leading term $eta_1 V_1$ captures the greatest variation of any such one-dimensional component, and the subsequent components $eta_j V_j$ ($j=2, 3, ...$) capture progressively smaller fractions of the total variance. Since the $eta_j$ are uncorrelated, these "axes of variation" capture distinct and unrelated contributions to the overall variance of $Y$.

== Biplots

Suppose that $tilde(X)$ is the column-centered or double-centered data (where the rows are observations and the columns are variables). Let $tilde(X) = U S V'$ denote the SVD. A _biplot_ is a plot that displays both the variables and the observations in a way that conveys (i) how the variables are related to each other, (ii) how the observations are related to each other, and (iii) how the variables are related to the observations. In other words, the biplot shows which variables load on each of the two displayed components ($j$ and $k$), and simultaneously shows which of the observations most exhibit the characteristics summarized by each component.

To construct the biplot, observation $i$ is plotted as the point $(S_(j j)^alpha U_(i j), S_(k k)^alpha U_(i k))$ (the _object scores_) and variable $ell$ is plotted as the point $(S_(j j)^(1-alpha) V_(ell j), S_(k k)^(1-alpha) V_(ell k))$ (the _variable scores_). Here, $j, k$ are chosen components, usually $j=1$ and $k=2$ (to show the dominant factors). The parameter $0 <= alpha <= 1$ is usually set to either $0$, $frac(1, 2, style: "horizontal")$, or $1$.

If we set $alpha=1$, the object scores approximate Euclidean distances in the
data. For example, let $d = (1, -1, 0, 0, ...)'$, and note that

$
norm(d' tilde(X))^2 = d' tilde(X) tilde(X)' d = d' U S^2 U' d = norm(d' U S)^2.
$

This shows that the Euclidean distance between two rows of $tilde(X)$ is equal to the Euclidean distance between two rows of $U S$ (which are the object scores if $alpha=1$).

If we set $alpha=0$, then $n^(-1) tilde(X)' tilde(X)$, which is the covariance matrix among the variables, can be written (up to a proportionality constant) as

$
tilde(X)' tilde(X) = V S^2 V'.
$

This shows that the magnitudes of the variable scores are equal to the variances of the variables, and the cosine of the angles between any two variables scores is equal to the correlation between two variables.

The interpretations of the object and variable scores given above are exact if we consider the non-truncated vectors of object/variable scores. In practice we truncate to the dominant two components, to enable two-dimensional plotting. In this case, the statements above are approximations. In particular, the magnitudes of the variable scores reflect the variance of a variable that is reflected in the biplot, rather than its overall variance.

Another pragmatic issue in biplots is that theoretically we should set $alpha=1$ if we want to focus on the objects and $alpha=0$ if we want to focus on the variables. In practice we may set $alpha=frac(1, 2, style: "horizontal")$ to approximate both sets of relationships in the same plot.

== Principal Components Regression

Principal Components Analysis can be used to understand relationships in multivariate data when there is no single variable that can be seen as the "response" relative to the other variables. PCA can also be used to create covariates for use in a regression analysis, which is called "Principal Components Regression" (PCR). The motivation behind PCR is that if we have a large number of covariates and do not wish to explicitly include all of them in a regression model, we can use PCA to reduce the variables to a smaller vector of scores that capture most of the variation in the original variables, and then use the reduced variables (scores) as covariates in our regression.

More formally, suppose that $X in RR^(n times p)$ is our regression design matrix, for a regression with $n$ observations (cases) and $p$ covariates. The response variable for this regression is a vector $Y in RR^n$, but this is not used in the first few steps of PCR. First, the mean should be subtracted from each column of $X$, and we write $X = U S 'T$ using the SVD. Next we truncate this representation so that $tilde(U) in RR^(n times q)$, contains the first $q$ columns of $U$. As discussed above, these are the columns of $U$ that explain the most variation in $X$.

We use the columns of $tilde(U)$ instead of the columns of $X$ as covariates in our regression, obtaining a linear predictor $tilde(U) hat(beta)$, where $hat(beta)$ are coefficients obtained using a regression procedure, e.g. least squares or a GLM. To relate the coefficients $hat(beta)$ back to the original covariates, note that $U = X V S^{-1}$, and $tilde(U) = X(V S^(-1))_(:,1:q)$. Therefore

$
tilde(U) hat(beta) = X(V S^(-1))_(:,1:q) hat(beta),
$

and we can define

$
hat(gamma) = (V S^(-1))_(:,1:q)hat(beta).
$

The coefficients $hat(gamma)$ align with the original covariates, while the coefficients $hat(beta)$ align with the principal component scores. Both determine the same fitted values and linear predictor.

PCR can be very effective, but may or may not peform well depending on the circumstances. In PCR, instead of allowing the coefficients $gamma$ to take on any value in $RR^p$, we are constraining $gamma$ to lie in the columnspace of $(V S^(-1))_(:,1:q)$. PCR works well if the projections of $X$ that explain the most variation in $X$ also explain the most variation in $Y$. This is often, but not always the case.

= Canonical Correlation Analysis

Suppose we have two vectors (of measurements or variables) and a collection of objects that are assessed with respect to all variables in both vectors. Canonical Correlation Analysis (CCA) aims to understand the relationships between the two vectors. This goal differs from the goal of PCA, which is used to understand the relationships within each vector of variables separately. Specifically, in CCA we can search for linear combinations of the components of each vector that are maximally correlated with each other.

Let $X in RR^(n times p)$ and $Y in RR^(n times q)$ denote the data for the two collections of variables. In these data, we have measurements on $n$ objects (the rows of $X$ and $Y$), a vector of $p$ variables in the columns of $X$, and a vector of $q$ variables in the columns of $Y$. Row $i$ of $X$ and row $i$ of $Y$ correspond to the same object.

Given coefficient vectors $a in RR^p$ and $b in RR^q$, we can form linear predictors $X a$ and $Y b$, both of which are $n$-dimensional vectors. The goal of CCA is to find $a$ and $b$ to maximize the correlation coefficient between $X a$ and $Y b$. Inspecting plots of the canonical loadings $a$ and $b$ can yield insight into the relationship between the variables in $X$ and the variables in $Y$.

Specifically, taking $X$ and $Y$ to be column-centered, the objective function is

$
(a^T S_(x y) b) / ((a'S_(x x)a)^(1/2) dot.c (b' S_(y y)b)^(1/2)),
$

where $S_(x y) = X' frac(Y, n, style: "horizontal")$, $S_(x x) = X' frac(X, n, style: "horizontal")$, and $S_(y y) = Y' frac(Y, n, style: "horizontal")$ are the cross-correlation and correlation matrices.

Once we have found the leading canonical loadings, say $a^((1))$ and $b^((1))$, we can maximize the same objective function subject to the constraints $a' S_(x x)a^((1)) = 0$ and $b' S_(y y)b^((1)) = 0$. This can be repeated up to $"min"(p, q)$ times to produce a series of canonical variate loading pairs

If the dimension is high, CCA can overfit the data. A procedure analogous to principal components regression can be used to address this. We first reduce $X$ and $Y$ to lower dimensions using PCA, then we fit the CCA using these reduced data. The coefficients from the reduced PCA can then be re-expressed in the original coordinates for interpretation.

= Dimension Reduction Regression

Dimension reduction regression (DRR) is a flexible, nonparametric approach to regression analysis that can be seen as extending PCA to the setting of regression. At the population level in a DRR analysis, we observe copies of $(Y, X)$, where both $Y in RR$ and $X in RR^p$ are random.  The goal is to find a matrix $B in RR^(p times q)$ such that

$
Y indep X | B'X.
$

This is a very general and nonparametric way to express that all of the information in the predictors $X$ about the outcome $Y$ is contained in the reduced variates $B'X$.  Only the columnspace $"col"(B)$ is identified, and we seek a $B$ of minimal column rank that achieves this independence condition.  This minimal subspace is known as the _sufficient dimension reduction_ (SDR) subspace.

In an empirical DRR analysis with data, we have a matrix $X in RR^(n times p)$ containing data on the explanatory (independent) variables, and a vector $Y in RR^n$ containing the values of a response (dependent) variable. Like PCA, the goal is to find linear _factors_ (also known as _components_ or _variates_) among the $X$ variables, but in the case of DRR the goal is for these factors to predict $Y$, not to predict $X$ itself. One way to view DRR is as a means to "steer" PCA toward variates that explain the variation in $Y$, rather than seeking variates that explain the variation in $X$.

In some cases, the variates in $X$ that explain $X$ and the variates in $X$ that explain $Y$ can be quite similar, but in other cases they dramatically differ. This is a common critique of principal component regression (PCR), which uses the PCs as explanatory variables and therefore can only find regression relationships that happen to coincide with the dominant modes of variation of $X$ (the PCs).

We can also view DRR \as fitting the multi-index regression model

$
E[Y | X=x] = f(b_1^T X, ..., b_q^T X),
$

where the $b_j in RR^p$ are coefficient vectors defining the "indices" in $X$ that predict $Y$. The link function $f$ must be non-linear, since in a linear $f$, the multiple indices would collapse to a single index.

An appealing feature of the DR approach is that it incorporates a link function $f$ but this function does not need to be known. Thus, DR provides a means to capture a wide range of nonlinear and non-additive regression relationships. The parameter $q$ above determines the dimension of the dimension reduction subspace, and is selected based on the data. If $q=1$ we have a single-index model, like in a GLM, except that here the link function does not need to be pre-specified. As $q$ grows, the model becomes more complex. DR is most effective when relatively small values of $q$ can be used, thereby compressing the regression structure into a few variates.

It is advantageous that the link function $f$ does not need to be known while estimating the coefficients $b_j$. However later in the analysis we may wish to estimate $f$, and commonly a nonparametric regression method like
#link("https://en.wikipedia.org/wiki/Local_regression")[loess] can be used for this purpose. Since $f$ operates on the dimensionally reduced variates, rather than on the full set of candidate covariates, the _curse of dimensionality_ is partially overcome, potentially allowing basic nonparametric methods to be employed for estimation of $f$.

One of the simplest and most widely-used approaches to dimension reduction regression is Sliced Inverse Regression (SIR). Recall that the loading vectors of PCA are the eigenvectors of $Sigma_x = "Cov"(X)$. In SIR, we wish to steer the PCA loadings toward directions that are relevant for predicting $Y$. To do this we consider $M_(x|y) equiv "Cov"[E[X|Y]]$, which is a covariance matrix for $X$ that suppresses the variation in $X$ that is irrelevant for $Y$. Specifically, this is accomplished by replacing $X$ in $"Cov"[X]$ with $E[X|Y]$, which suppresses the variation in $X$ that occurs when $Y$ is fixed.

The relationship between $Sigma_x$ and $M_(x|y)$ can be further elucidated by
employing the
#link("https://en.wikipedia.org/wiki/Law_of_total_covariance")[law of total covariance],
which asserts that

$
"Cov"(X) = "Cov" E[X|Y] + E "Cov"(X|Y).
$

Equivalently, $Sigma_x = M_(x|y) + R$, where $R = E["Cov"(X|Y)]$. The matrices $Sigma_x = "Cov"[X]$, $M_(x|y) = "Cov"[E[x|y]]$, and $R = E["Cov"[X|Y]]$ are all
#link("https://en.wikipedia.org/wiki/Positive_semidefinite_matrix")[positive semidefinite]. Therefore, the matrix $M_(x|y)$ is smaller than $Sigma_x$ in the sense of
#link("https://en.wikipedia.org/wiki/Loewner_order")[Loewner ordering].

The solutions to the generalized eigenvalue problem $Sigma_x b = lambda M_(x|y)b$ correspond to the solutions of the Rayleigh quotient optimization

$
"argmax"_b (b^T M_(x|y) b) / (b^T Sigma_x b),
$

which in turn is equivalent to the solutions of

$
"argmax"_b "Var"[E[b' X|Y]] / "Var"[b'X].
$

This last expression provides the best intuition for how SIR works -- we are seeking direction vectors such that if we project the data in the given direction, the explained variance $"var"[E[b^T X|Y]]$ is maximized relative to the total variance $"Var"[b'x]$.

The matrix $M_(x|y)$ is estimated by sorting the data ${(x_i, y_i)}$ so that $y$ is non-decreasing, dividing this sorted sequence into slices (blocks), and averaging the values of $x_i$ within each slice. Let $u_k in RR^d$ denote the mean of slice $k$. We then estimate $M_(x|y)$ as the covariance matrix of the $u_k$ (or as the weighted covariance with weights proportional to the number of observations in each slice).

= Correspondence Analysis

Correspondence analysis (CA) is an embedding approach that aims to represent _chi-square distances_ in the data space as Euclidean distances for visualization. The chi-square metric is defined as follows: if $X, Y in RR^p$ are random vectors with mean $E[X] = E[Y] = mu$ and covariance $"Cov"(X) = "Cov"(Y) = "diag"(mu)$, then the (squared) chi-square distance from $X$ to the mean is $(X - mu)'"diag"(mu)^(-1)(X-mu)$, and the squared chi-square distance between $X$ and $Y$ is $(X-Y)'"diag"(mu)^(-1)(X-Y)$. As discussed further below, correspondence analysis arises naturally when working with nominal (categorical) variables, but there are other situations where it makes sense to apply CA when the data are not nominal.

The motivation for embedding using the chi-square metric this is that in many settings chi-square distances may best represent information in the data, while Euclidean distances and angles are arguably the best approach for producing visualizations for human interpretation. As discussed further below, the chi-square metric is especially appropriate for understanding relationships among the objects or observations (the rows of the data matrix), which often can be viewed as an independent and identically distributed collection of draws from a multivariate distribution.

== Mean/variance relationships

To better understand the motivation for using chi-square distances to measure distances in the data space, let $X in RR^p$ be a random vector with mean $mu >= 0$ and covariance matrix $Sigma$. In some cases, $mu$ and $Sigma$ are unrelated (i.e. knowing $mu$ places no constraints on $Sigma$, and vice-versa). On the other hand, in many settings it is plausible that $mu$ and $Sigma$ are related via a _mean-variance relationship_. A common form of mean-variance relationship is $"diag"(Sigma) prop mu$. A particular example is the Poisson distribution where $"diag"(Sigma) = mu$, so the constant of proportionality is $1$. In a broader class of settings we may have quasi-Poisson _over-dispersion_ or _under-dispersion_, meaning that $Sigma_(i i) = c dot.c mu_i$, where $c > 1$ or $c < 1$ for over and under-dispersion, respectively.

When comparing two random vectors $X$, $Y$ sharing covariance matrix $Sigma$, it is reasonable to use the #link("https://en.wikipedia.org/wiki/Mahalanobis_distance")[Mahalanobis distance].

$
d(X, Y)^2 = (X - Y)^T Sigma^(-1) (X - Y).
$

The chi-square distance is a special case of the Mahalanobis distance in which $Sigma = "diag"(mu)$.

== Multiple Correspondence Analysis (MCA)

Suppose we have $n$ observations on $p$ variables, and the data are represented in an $n times p$ matrix $X$ whose rows are the cases (observations) and columns are the variables. Correspondence analysis can be applied when each $X_(i j) >= 0$, and where it makes sense to compare any two rows or any two columns of $X$ using chi-square distance.

Let $P equiv frac(X, N, style: "horizontal")$, where $N = sum_(i j) X_(i j)$. We introduce the following notation: let $P_(i,:)$ denote row $i$ of $P$, let $r equiv P dot.c 1_p$ denote the row sums of $P$, and let $c = P' dot.c 1_n$ denote the column sums of $P$. Also, let $W_r = "diag"(r) in RR^(n times n)$ and $W_c = "diag"(c) in RR^(p times p)$. Then $P^r equiv W_r^(-1) dot.c P$ are the _row profiles_ of $P$, which are simply the rows of $P$ (or of $X$) normalized by their sum. Analogously, $P^c equiv P dot.c W_c^(-1)$ are the _column profiles_ of $P$ (or of $X$).

The primary goal of MCA is to embed the row profiles $P^r$ into _row scores_ $F$ and the column profiles $P^c$ into _column scores_ $G$, where $F$ is an $n times p$ array and $G$ is a $p times p$ array, and both embeddings approximate chi-square distances. For any $j=1, 2, ..., p$, the $j^"th"$ columns of $F$ and $G$ comprise an MCA factor.

In many cases, the row sums will be approximately constant, so $W_r prop 1_n$ and this matrix can effectively be ignored. However it is very unlikely that $W_c prop 1_p$, and in fact the inhomogeneity in $c$ is the key feature in the data that motivates use of MCA.

Our goals are as follows:

- For any $1 <= i, j <= n$, the Euclidean distance from $F_(i,:)$ to
  $F_(j,:)$ is equal to the chi-square distance from $P^r_(i,:)$ to
  $P^r_(j,:)$. Also, for any $1 <= i,j <= p$ the Euclidean distance from
  $G_(i,:)$ to $G_(j,:)$ is equal to the chi-square distance from $P^c_(:,i)$ to $P^c_(:,j)$. Thus, $F$ provides an embedding of the rows of $P^r$ and $G$ provides an embedding of the columns of $P^c$. These embeddings map chi-square distances in the data space (specifically, the space of row profiles or column profiles) to Euclidean distances in the embedded space.

- The columns of $F$ and $G$ are ordered in terms of importance. Specifically, if we select $1 <= q <= p$ then the Euclidean distance from $F_(i,1:q)$ to $F_(j,1:q)$ is the best possible $q$-dimensional approximation to the chi-square distance from $P_(i,:)$ to $P_(j,:)$. Note that if $q=p$ then the approximation becomes exact, but for $q < p$ the approximation is inexact.

== Derivation of the algorithm

The array $P - r c'$ is a doubly-centered array of residuals, since $1'_n (P - r c') = 0_p$, and $(P - r c')1_p = 0_n$. We can further standardize these residuals by scaling the rows and columns by their respective standard deviations (taking the variance to be proportional to the mean). This gives us the doubly-standardized array of residuals

$
W_r^(-1/2)(P - r c')W_c^(-1/2).
$

We next take the singular value decomposition of these doubly-standardized residuals:

$
W_r^(-1/2)(P - r c')W_c^(-1/2) = U S V'.
$

Now let $F = W_r^(-1/2) U S$ and $G = W_c^(-1/2) V S$. We show that this specification of $F$ and $G$ satisfies the conditions stated above. Specifically, we will show that the rows of $F$ embed the objects (the rows of $X$), converting chi-square distances among the rows of $X$ to Euclidean distances among the rows of $F$. If we make a scatterplot of the points $(F_(i 1), F_(i 2))$, then the Euclidean distances among these points approximate the chi-square distances among the rows of $X$, providing a meaningful Euclidean embeding of the rows of $X$.

First, note that since $V$ is orthogonal

$
norm(F_(i,:) - F_(j,:)) = norm(F_(i,:) - F_(j,:))V'.
$

Therefore,

$
norm(F_(i,:) - F_(j,:))^2 =&
norm(r_i^(-1/2)U_(i,:)S - r_j^(-1/2)U_(j,:)S)^2\
=& norm(r_i^(-1)(P_(i,:) - r_i c')W_c^(-1/2) - r_j^(-1)(P_(j,:) - r_j c^T)W_c^(-1/2))^2\
=& norm(r_i^(-1)P_(i,:)W_c^(-1/2) - r_j^(-1)P_(j,:)W_c^(-1/2))^2\
=& (P_(i,:)/r_i - P_(j,:)/r_j)' W_c^(-1)(P_(i,:)/r_i - P_(j,:)/r_j)\
=& (P^r_i - P^r_j)' W_c^(-1) (P^r_i - P^r_j)\
=& n^(-1)(P^r_i - P^r_j)' (W_c/n)^(-1) (P^r_i - P^r_j).
$

Suppose for simplicity that $W_r prop I_n$. Since $frac(W_c, n, style: "horizontal") = "diag"(frac(c, n, style: "horizontal"))$, where $frac(c, n, style: "horizontal")$ is the mean of the rows of $P^r$, it follows that $norm(F_(i,:) - F_(j,:))$ is $n^(-1)$ times the squared chi-square distance between $P^r_(i,:)$ and $P^r_(j,:)$. Thus, the rows of $F$ embed the rows of $P^r$ as desired. Applying the same argument to $X^T$ shows that the rows of $G$ embed the columns of $P^c$, also reflecting chi-square distances.

== Correspondence analysis and Multiple Correspondence analysis for nominal data

One common application of correspondence analysis (CA) arises when analyzing datasets in which all variables are nominal. First, suppose we have a single nominal variable and code it using an _indicator matrix_. That is, we define $X$ to be a matrix whose values are entirely $0$ and $1$, such that $X_(i j) = 1$ if and only if the value of the nominal variable for case $i$ is equal to level $j$. Correspondence analysis as defined above can be used to analyze this indicator matrix, revealing how the objects and categories are related. Note that in this case $r = 1_n$ and $W_r = I_n$.

An important extension of CA is _Multiple Correspondence Analysis_, in which we have several nominal variables. In this case, we recode each nominal variable with its own indicator matrix, and then concatenate these matrices horizontally. If there are $p_j$ levels for variable $j$, and we set $p = sum_j p_j$, then the concatenated indicator matrix (the _Burt matrix_) is $n times p$. We then apply CA to this concatenated indicator matrix, yielding insights into the relationships among the objects, and relationships among levels of different variables.

Note that in the special case where MCA is applied to a collection of nominal variables, the Burt matrix $X$ has the property that every row sums to $q$ (where $q$ is the number of nominal variables). Thus every row of $P = frac(X, (n q), style: "horizontal")$ sums to $frac(1, n, style: "horizontal")$ and the column means of $P^c$ are also all equal to $frac(1, n, style: "horizontal")$, since

$
overline(P)^c = n^(-1)1'_n P Q_c^(-1) = n^(-1)1'_p.
$

== Angles and magnitudes of category scores

In many applications of CA and MCA, the rows are viewed as an IID sample from a distribution, and the variables of this distribution exhibit a Poisson mean/variance relationship. This motivates embedding the row profiles $P^r$ derived from the data matrix $X$ using chi-square distances. However the argument for embedding the columns (variables) using the chi-square metric is less clear (although it is true that MCA achieves such an embedding). A more interpretable embedding of the column profiles is based on covariances and correlations rather than on chi-square distances. Here we show that MCA also achieves a column embedding that can be interpreted in terms of covariances and correlations.

The angle between two embedded category points is an easily intepreted visual feature of the plot of category embeddings (in practice usually restricted to the dominant two components). The angle $theta$ between two vectors $v$, $w$ satisfies

$
lr(chevron.l v, w chevron.r) = norm(v) dot.c norm(w) dot.c cos(theta).
$

Therefore, to understand the angles among the embedded categories, we should consider the Gram matrix $G G'$ containing inner products of every pair of embedded category points.

The matrix $G G'$ has the following relationship to the data in $P$:

$
G G' = W_c^(-1/2) V S S V' W_c^(-1/2) = W_c^(-1)(P - r c')'W_r^(-1)(P - r c')W_c^(-1).
$

In most applications of MCA, $W_r = q I_n$, where $q$ is the number of
variables. Further, as noted above $overline(P)^c_(:,i) = 1/n$ for each $i$, and $r = n^(-1)1_n$. Therefore, a single element of the Gram matrix has the form

$
[G G']_(i j) = q^(-1)(P^c_(:,i) - r)'(P^c_(:,j) - r) = q^(-1)n dot.c "Cov"(P^c_(:,i), P^c_(:,j)).
$

The diagonal elements of the Gram matrix are proportional to the variances of the column profiles:

$
[G G']_(i i) = q^(-1)n dot.c "Var"(P^c_(:,i)).
$

It follows that the angle between two category profiles, say in columns $i$ and $j$, is proportional to the correlation coefficient between the column profiles $P^c_i$ and $P^c_j$.

In addition, a category point lies further from the origin if the column profiles are not close to being constant vectors. These are the column profiles that most strongly influence the MCA fit.

Combining these observations, two column profiles corresponding to categories of different variables that are both far from the origin, and that have a small angle, correspond to substantially correlated indicator vectors.

Note that the interpretation of the category points in MCA usually focuses on the relationships between pairs of categories for different variables. Comparing two categories of the same variable is rarely interesting since these indicators are by definition mutually exclusive.

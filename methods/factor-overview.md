# Factor analyses and embeddings

A large class of powerful statistical methods considers data in which many
"objects" ("observations") are measured with respect to multiple variables.
These analyses often focus on understanding both the relationships among the
variables and the relationships among the observations. The data takes the
form of a rectangular array, where by convention the rows correspond to
objects and the columns correspond to variables.

For example, we may have a sample of people (the "objects") with each person
measured in terms of their height, weight, age, and sex. In this example, two
people are similar if they have similar values on most or all of the
variables, and two variables are similar if knowing the value of just one of
the variables for a particular object allows one to predict the value of the
other variable for that same object.

Matrix factorization is a powerful idea from linear algebra, and underlies
many of the most important methods in statistics. We will use the term "factor
analysis" here to refer in a broad way to any statistical method relying on a
linear algebraic factorization of the data matrix (or a transformed version of
the data matrix).

To visualize such a data matrix, we can think in terms of _variable space_ and
_case space_ (or _object space_). If there are $n$ cases and $p$ variables,
then the cases are points in ${\mathbb R}^p$ and the variables are points in
${\mathbb R}^n$.

## Embedding

Many methods of factor analysis can be viewed as a way to obtain an
*embedding*. Embedding algorithms take input data vectors $x$ and transform
them into output data vectors $z$. Many embedding algorithms take the form of
a "dimension reduction", so that
$q \equiv {\rm dim}(z) < {\rm dim}(x) \equiv p$. However some embeddings
preserve or even increase the dimension.

Embeddings can be used for exploratory analysis, especially in visualizations,
and can also be used to construct features for prediction, as well as being
used in formal statistical inference, e.g. in hypothesis testing.

An embedding approach is linear if $z = B^Tx$ for a fixed $p\times q$ matrix
$B$. The matrix $B$ may be based on training data, in which case the embedding
is *adaptive*. In other settings $B$ is completely independent of any data and
the embedding is *non-adaptive*. Linear embedding algorithms are simpler to
devise and characterize than nonlinear embeddings. Many modern embedding
algorithms are nonlinear, and use this additional flexibility to better
capture complex structure.

Some embedding algorithms embed only the objects while other embedding
algorithms embed both the objects and the variables. Embedding the objects
provides a reduced feature representation for each object that can be passed
on to additional analysis procedures, or interpreted directly. Embedding the
variables provides a means to interpret the relationships among the variables.
When utilizing an embedding algorithm that embeds both objects and variables
we have the opportunity to interpret the results through a
[biplot](https://en.wikipedia.org/wiki/Biplot).

## Centering and standardization

Suppose that $X$ is a $n\times p$ matrix whose rows are the observations or
objects, and whose columns are the variables. There are many ways to
pre-process the data in $X$ prior to performing a factor analysis.

Some factor-type methods work with the covariance matrix of the variables,
which is a $p\times p$
[positive semidefinite](https://en.wikipedia.org/wiki/Definite_matrix) (PSD)
matrix. Since covariances by definition are derived from mean centered
variables, it is common to center the variables (columns) of the data matrix
$X$ prior to conducting a factor analysis. Centering usually involves
subtracting the mean from each column, but in some cases the columns may be
centered with respect to another measure of location such as the median. After
this centering, the case vectors become a point cloud that is centered around
the origin in $R^p$.

We may wish for our results to be independent of the units in which each
variable was measured, and remove any influence of differing dispersions of
the variables on our results. This motivates standardizing the columns of $X$
prior to performing a factor analysis. This *standardization* involves first
centering each column as discussed above, then dividing each column by a
measure of scale such as the standard deviation or inter-quartile range.

In some datasets there is strong heterogeneity among the observations (rows)
and it is desirable to supress this in the analysis. In this case we may
choose to *double center* the data so that every row and every column has mean
zero. To achieve this, the following three steps can be performed: (i) center
the overall matrix around its grand mean, (ii) center each row of the matrix
with respect to the row mean, (iii) center each column of the matrix with
respect to the column mean.

## Singular Value Decoposition

Many embedding methods make use of a matrix factorization known as the
*Singular Value Decomposition* (SVD). The SVD is defined for any $n\times p$
matrix $X$. In most cases we want $n \ge p$, and if $n < p$ we can take the
SVD of $X^T$. When $n\ge p$, we decompose $X$ as $X = USV^T$, where $U$ is
$n\times p$, $S$ is $p\times p$, and $V$ is $p\times p$. The matrices $U$ and
$V$ are orthogonal so that $U^TU = I_p$, $V^TV = I_p$, and $S$ is diagonal
with $S_{11} \ge S_{22} \ge \cdots \ge S_{pp}$. The values on the diagonal of
$S$ are the *singular values* of $S$, and the SVD is unique except when there
are ties among the singular values or if one or more of the singular values
are zero.

One use of the SVD is to obtain a low rank approximation to a matrix $X$. The
extreme example of a low rank matrix is a matrix with rank 1, which means that
we can write $X_{ij} = r_i\cdot c_j$ for vectors $r\in {\mathbb R}^n$ and
$c\in {\mathbb c}^p$. Our goal will be to decompose $X$ as the (approximate)
sum of $k$ rank-1 matrices, and these matrices are ordered by their importance
in describing the variation in $X$.

Suppose we truncate the SVD using only the first $k$ components, so that
$\tilde{U}$ is the $n\times k$ matrix consisting of the leading (left-most)
$k$ columns of $U$, $\tilde{S} = (S_1, \ldots, S_k)$, and $V$ is the
$p\times k$ matrix consisting of the leading $k$ columns of $V$. In this case,
the matrix
$\tilde{X} \equiv \tilde{U}\cdot {\rm diag}(\tilde{S})\cdot \tilde{V}^T$ is a
rank $k$ matrix (it has $k$ non-zero singular values). According to the
[Eckart-Young](https://en.wikipedia.org/wiki/Low-rank_approximation) theorem,
among all rank $k$ matrices, $\tilde{X}$ is the closest matrix to $X$ in the
*Frobenius norm*, which is defined as

$$
{\rm Frob}(X)^2 = \lVert X \rVert_F^2 \equiv \sum\_{ij}X\_{ij}^2 = {\rm trace}(X^\prime X).
$$

Thus we have

$$
\tilde{X} = {\rm argmin}_{A: {\rm rank}(A) = k} \lVert A - X \rVert_F.
$$

Finding a low rank approximation is a type of least square problem, although
it is not the linear least squares used in regression analysis.

### Analyzing a data matrix using the SVD

Suppose we have a data matrix $X$ and wish to understand its structure. We can
begin my double centering the matrix $X$ as discussed above:

$$
X_{ij} = m + r_i + c_j + R_{ij}.
$$

Here, $m$ is the grand mean of $X$, $r$ contains the row means of $X-m$, $c$
contains the column means of $X-m$, and $R_{ij} \equiv X_{ij} - m - r_i - c_j$
are residuals. Next we can take the SVD of $R$, yielding
$R = U\cdot {\rm diag}(S)\cdot V^T$, or equivalently

$$
R_{ij} = \sum_{k=1}^p S_k U_{ik}V_{jk}.
$$

Although there are $p$ terms in the SVD of $R$, the first few terms may
capture most of the structure, so

$$
R_{ij} \approx \sum_{k=1}^q S_k U_{ik}V_{jk}
$$

for $q < p$ (this approximation holds better when the *tail singular values*
$S_{q+1}, \ldots, S_p$ are small). One important property that results from
calculating the SVD for a double-centered matrix is that $U_{\cdot k} = 0$ and
$V_{\cdot k} = 0$. That is, the columns of $U$ and $V$ are centered. This
column centering means that the SVD captures "deviations from the mean"
represented by the additive model $m + r_i + c_j$.

#### Assessing dimensionality

The rate at which the singular values decay contains important information
about the intrinsic dimensionality of the population under study. Let $s_j$
denote the $j^{\rm th}$ singular value, where $s_1 \ge s_2 \cdots$. Two
canonical patterns of decay that may be found are an exponential decay
$s_j \approx a\cdot \exp(-b\cdot j)$ and a power-law decay
$s_j \approx a/j^b$. To assess the decay of the eigenvalues, we can consider
plots of $\log s_j$ against $j$, which is linear in the case of exponential
decay, and plots of $\log s_j$ against $\log j$, which is linear in the case
of power-law decay.

## Principal Components Analysis

Suppose that $x$ is a $p$-dimensional random vector with mean $0$ and
covariance matrix $\Sigma$ (our focus here is not the mean, so if $x$ does not
have mean zero we can replace it with $x-\mu$, where $\mu=E[x]$). Principal
Components Analysis (PCA) seeks a linear embedding of $x$ into a lower
dimensional space of dimension $q < p$. The standard approach to PCA gives us
a $p\times q$ orthogonal matrix $B$ of *loadings*, which can be used to
produce *scores* $Q = B^Tx$, which are $q-$dimensional vectors.

For a single vector $x$, the scores are obtained via the mapping
$Q(x) = B^Tx$. For a data matrix $X$ whose rows are independent and
identically distributed (IID) copies of the random vector $x$, the scores can
be obtained via the mapping $Q = XB$, where each row of $Q$ contains the
scores for the corresponding row of $X$.

PCA can be viewed in terms of linear compression and decompression of the
variables in $X$. Let

$$
X \rightarrow B_{:,1:q}^TX \equiv Q(X)
$$

denote the compression that reduces the data from $p$ dimensions to $q$
dimensions, where $B_{:,1:q}$ is the $p\times q$ matrix consisting of the
leading $q$ columns of $B$. We can now decompress the data as follows:

$$ Q \rightarrow B_{:,1:q}Q \equiv \hat{X}.  $$

This represents a two-step process of first reducing the dimension, then
predicting the original data using the reduced data. Among all possible
loading matrices $B$, the PCA loading matrix loses the least information in
that it minimizes the expected value of $\lVert X - \hat{X} \rVert$.

The loading matrix $B$ used in PCA is a truncated eigenvector matrix of
$\Sigma = {\rm cov}(x)$. Specifically, we can write $\Sigma = B\Lambda B^T$,
where $B$ is an orthogonal matrix and $\Lambda$ is a diagonal matrix with
$\Lambda_{11} \ge \Lambda_{22} \ge \cdots \ge \Lambda_{pp} > 0$. This is the
spectral decomposition of $\Sigma$. The columns of $B$ are orthogonal in both
the Euclidean metric and in the metric of $\Sigma$, that is, $B^TB = I_p$ and
$B^T\Sigma B = I_p$. As a result, the scores $Q \equiv B^TX$ are uncorrelated,
${\rm cov}(Q) = I_q$.

Next we consider how PCA can be carried out with a sample of data, rather than
in a population. Given a $n\times p$ matrix of data $X$ whose rows are iid
copies of the random vector $x$, we can estimate the covariance matrix
$\Sigma$ by column centering $X$ to produce $X_c \equiv X - 1_n\bar{x}^T$,
where $\bar{x} \in {\mathbb R}^p$ is the vector of column-wise means of $X$,
Then set $\hat{\Sigma} = X_c^TX_c/n$. Letting $B$ denote the eigenvectors of
$\hat{\Sigma}$, the scores have the form $Q = X_cB$.

Since the eigenvalues $\Lambda_{ii}$ are non-increasing, the leading columns
of $Q$ contain the greatest fraction of information about the distributiuon of
$x$. Thus, visualizations (e.g. scatterplots) of the first two columns of $Q$
best reflect the relationships in $X$ (compared to any other scatterplot
formed from linear scores).

PCA can also be carried out using the SVD of the column-centered data matrix
$X_c$. Write $X_c = USV^T$ as discussed above and note that

$$
\hat{\Sigma} = X_c^T X_c/n = VS^2V^T/n.
$$

Thus, $\hat{\Sigma}V = VS^2/n$, so $V$ contains the eigenvectors of
$\hat{\Sigma}$ and the corresponding eigenvalues are in ${\rm diag}(S^2)/n$.

### Karhunen-Loeve decomposition

PCA can be seen as a way to estimate the terms of the
[Karhunen-Loeve (K-L) decomposition](https://en.wikipedia.org/wiki/Kosambi%E2%80%93Karhunen%E2%80%93Lo%C3%A8ve_theorem)
of a multivariate distribution (a random vector). If $Y$ is a $p$-dimensional
random vector, then the K-L decomposition expresses $Y$ in the form

$$
Y = \mu + \sum_{j=1}^p \eta_j V_j,
$$

where $\mu = E[Y]$ is a fixed $p$-dimensional vector containing the
element-wise mean of $Y$, the $\eta_j$ are a collection of mutually
uncorrelated random scalar variables, and the $V_j$ are a collection of
mutually orthogonal fixed vectors of dimension $p$. The means of the $\eta_j$
are all zero and the variances of the $\eta_j$ are non-increasing in $j$. That
is, $E[\eta_j] = 0$, ${\rm Var}(\eta_1) \ge {\rm Var}(\eta_2) \ge \cdots$, and
${\rm cor}(\eta_j, \eta_k) = 0$ if $j \ne k$.

Each component $\eta_j V_j$ represents random variation in the direction of
the unit vector $V_j$. The leading term $\eta_1 V_1$ captures the greatest
variation of any one-dimensional component, and the subsequent components
$\eta_j V_j$ ($j=2, 3, \ldots$) capture progressively less of the variance.
Since the $\eta_j$ are uncorrelated, these "axes of variation" capture
distinct and unrelated contributions to the overall variance of $Y$.

### Biplots

Suppose that $X_c$ is the column-centered or double-centered data (where the
rows are observations and the columns are variables). Let $X_c = USV^T$ denote
the SVD. A *biplot* is a plot that displays both the variables and the
observations in a way that conveys (i) how the variables are related to each
other, (ii) how the observations are related to each other, and (iii) how the
variables are related to the observations. In other words, the biplot shows
which variables load on each of the two displayed components ($j$ and $k$),
and simultaneously shows which of the observations most exhibit the
characteristics summarized by each component.

To construct the biplot, observation $i$ is plotted as the point
$(S_{jj}^\alpha U_{ij}, S_{kk}^\alpha U_{ik})$ (the *object scores*) and
variable $\ell$ is plotted as the point
$(S_{jj}^{1-\alpha} V_{\ell j}, S_{kk}^{1-\alpha} V_{\ell k})$ (the *variable
scores*). Here, $j, k$ are chosen components, usually $j=1$ and $k=2$ (to show
the dominant factors). The parameter $0 \le \alpha \le 1$ is usually set to
either $0$, $1/2$, or $1$.

If we set $\alpha=1$, the object scores approximate Euclidean distances in the
data. For example, let $d = (1, -1, 0, 0, \ldots)^\prime$, and note that

$$
\lVert d^\prime X_c \rVert^2 = d^\prime X_cX_c^\prime d =
d^\prime US^2U^\prime d = \lVert d^\prime US \rVert^2.
$$

This shows that the Euclidean distance between two rows of $X_c$ is equal to
the Euclidean distance between two rows of $US$ (which are the object scores
if $\alpha=1$).

If we set $\alpha=0$, then $n^{-1}X_c^\prime X_c$, which is the covariance
matrix among the variables, can be written (up to a proportionality constant)
as

$$
X_c^\prime X_c = VS^2V^\prime.
$$

This shows that the magnitudes of the variable scores are equal to the
variances of the variables, and the cosine of the angles between any two
variables scores is equal to the correlation between two variables.

The interpretations of the object and variable scores given above are exact if
we consider the non-truncated vectors of object/variable scores. In practice
we truncate to the dominant two components, to enable two-dimensional
plotting. In this case, the statements above are approximations. In
particular, the magnitudes of the variable scores reflect the variance of a
variable that is reflected in the biplot, rather than its overall variance.

Another pragmatic issue in biplots is that theoretically we should set
$\alpha=1$ if we want to focus on the objects and $\alpha=0$ if we want to
focus on the variables. In practice we may set $\alpha=1/2$ to approximate
both sets of relationships in the same plot.

### Principal Components Regression

Principal Components Analysis is a method for multivariate data analysis that
can be used to understand relationships in multivariate data when there is no
single variable that can be seen as the "response" relative to the other
variables. PCA can also be used to create covariates for use in a regression
analysis. This method is called "Principal Components Regression" (PCR). The
motivation behind PCR is that if we have a large number of covariates and do
not wish to explicitly include all of them in a regression model, we can use
PCA to reduce the variables to a smaller vector of scores that capture most of
the variation in the original variables, and then use the reduced variables
(scores) as covariates in our regression.

More formally, suppose that $X \in {\mathbb R}^{n\times p}$ is our regression
design matrix, for a regression with $n$ observations (cases) and $p$
covariates. The response variable for this regression is a vector
$Y \in {\mathbb R}^n$, but this is not used in the first few steps of PCR.
First, the mean should be subtracted from each column of $X$, and we write
$X = USV^T$ using the SVD. Next we truncate this representation so that
$\tilde{U} \in {\mathbb R}^{n\times q}$, contains the first $q$ columns of
$U$. As discussed above, these are the most important columns of $U$ for
explaining the variation in $X$.

We use $\tilde{U}$ instead of $X$ as covariates in our regression, obtaining a
linear predictor $\tilde{U}\hat{\beta}$, where $\hat{\beta}$ are coefficients
obtained using a regression procedure, e.g. least squares or a GLM. To relate
the coefficients $\hat{\beta}$ back to the original covariates, note that
$U = XVS^{-1}$, and $\tilde{U} = X(VS^{-1})_{:,1:q}$. Therefore

$$
\tilde{U}\hat{\beta} = X(VS^{-1})_{:,1:q}\hat{\beta},
$$

and we can define

$$
\hat{\gamma} = (VS^{-1})_{:,1:q}\hat{\beta}.
$$

The coefficients $\hat{\gamma}$ align with the original covariates, while the
coefficients $\hat{\beta}$ align with the principal component scores. Both
determine the same fitted values and linear predictor.

PCR can be very effective, but may or may not peform well depending on the
circumstances. In PCR, instead of allowing the coefficients $\gamma$ to take
on any value in ${\mathbb R}^p$, we are constraining $\gamma$ to lie in the
columnspace of $(VS^{-1})_{:,1:q}$. PCR works well if the projections of $X$
that explain the most variation in $X$ also explain the most variation in $Y$.
This is often, but not always the case.

## Canonical Correlation Analysis

Suppose we have two vectors (of measurements or variables) and a collection of
objects that are assessed with respect to all variables in both vectors.
Canonical Correlation Analysis (CCA) aims to understand the relationships
between the two vectors. This goal differs from the goal of PCA, which is used
to understand the relationships within each vector of variables separately.
Specifically, in CCA we can search for linear combinations of the components
of each vector that are maximally correlated with each other.

Let $X \in {\mathbb R}^{n\times p}$ and $Y \in {\mathbb R}^{n\times q}$ denote
the data for the two collections of variables. In these data, we have
measurements on $n$ objects (the rows of $X$ and $Y$), a vector of $p$
variables in the columns of $X$, and a vector of $q$ variables in the columns
of $Y$. Row $i$ of $X$ and row $i$ of $Y$ correspond to the same object.

Given coefficient vectors $a\in {\mathbb R}^p$ and $b\in {\mathbb R}^q$, we
can form linear predictors $Xa$ and $Yb$, both of which are $n$-dimensional
vectors. The goal of CCA is to find $a$ and $b$ to maximize the correlation
coefficient between $Xa$ and $Yb$. Inspecting plots of the canonical loadings
$a$ and $b$ can yield insight into the relationship between the variables in
$X$ and the variables in $Y$.

Specifically, taking $X$ and $Y$ to be column-centered, the objective function
is

$$
\frac{a^\prime S_{xy} b}{(a'S_{xx}a)^{1/2}\cdot (b^\prime S_{yy}b)^{1/2}},
$$

where $S_{xy} = X^\prime Y/n$, $S_{xx} = X^\prime X/n$, and
$S_{yy} = Y^\prime Y/n$ are the cross-correlation and correlation matrices.

Once we have found the leading canonical loadings, say $a^{(1)}$ and
$b^{(1)}$, we can maximize the same objective function subject to the
constraints $a^\prime S_{xx}a^{(1)} = 0$ and $b^\prime S_{yy}b^{(1)} = 0$.
This can be repeated up to ${\rm min}(p, q)$ times to produce a series of
canonical variate loading pairs

If the dimension is high, CCA can overfit the data. A procedure analogous to
principal components regression can be used to address this. We first reduce
$X$ and $Y$ to lower dimensions using PCA, then we fit the CCA using these
reduced data. The coefficients from the reduced PCA can then be re-expressed
in the original coordinates for interpretation.

## Dimension reduction regression

Dimension reduction regression (DR) is a flexible approach to regression
analysis that can be seen as extending PCA to the setting of regression. In a
DR analysis, we have a matrix $X \in {\mathbb R}^{n\times p}$ containing data
on the explanatory (independent) variables, and we also have a vector
$Y \in {\mathbb R}^n$ containing values of a response (dependent) variable.
Like PCA, the goal is to find *factors*, *components*, or *variates* among the
$X$ variables, but in the case of DR the goal is for these factors to predict
$Y$, not to predict $X$ itself. One way to view DR is as a means to "steer"
PCA toward variates that explain the variation in $Y$, rather than seeking
variates that explain the variation in $X$.

In some cases, the variates in $X$ that explain $X$ and the variates in $X$
that explain $Y$ can be quite similar, but in other cases they dramatically
differ. This is a critique of principal component regression (PCR), which uses
the PCs as explanatory variables and therefore can only find regression
relationships that happen to coincide with the PCs.

DR can be viewed as fitting the multi-index regression model

$$
E[Y | X=x] = f(b_1^\prime x, \ldots, b_q^\prime x),
$$

where the $b_j \in {\mathbb R}^p$ are coefficient vectors defining the
"indices" in $X$ that predict $Y$. The span of $b_1, \ldots, b_q$ is known as
the "dimension reduction" subpsace for the regression function $E[Y|X=x]$.

An appealing feature of the DR approach is that it incorporates a link
function $f$ but this function does not need to be known. Thus, DR provides a
means to capture a wide range of nonlinear and non-additive regression
relationships. The parameter $q$ above determines the dimension of the
dimension reduction subspace, and is selected based on the data. If $q=1$ we
have a single-index model, like in a GLM, except that here the link function
does not need to be pre-specified. As $q$ grows, the model becomes more
complex and therefore may underperform due to high variance and overfitting.
DR is most effective when relatively small values of $q$ can be used, thereby
compressing the regression structure into a few variates.

It is advantageous that the link function $f$ does not need to be known while
estimating the coefficients $b_j$. However later in the analysis we will
probably want to estimate $f$, and commonly a nonparametric regression method
like loess can be used for this purpose.

One of the simplest and most widely-used approaches to dimension reduction
regression is known as Sliced Inverse Regression (SIR). Recall that the
loading vectors of PCA are the eigenvectors of ${\rm cov}(X)$. In SIR, we wish
to steer the PCA loadings toward directions that are relevant for predicting
$Y$. To do this we consider $M_{xy} \equiv {\rm cov}E[X|Y]$, which is also a
type of covariance matrix for $X$, but which suppresses the variation in $X$
that is irrelevant for $Y$. This is accomplished by replacing $X$ in
${\rm cov}(X)$ with $E[X|Y]$, which suppresses the variation in $X$ that
occurs when $Y$ is fixed. The matrix $M_{xy}$ is estimated by sorting the data
$\{(x_i, y_i)\}$ by increasing values of $y$, dividing this sorted sequence
into "slices" (blocks), and averaging the values of $x_i$ within each slice.
Let $u_k \in {\mathbb R}^d$ denote the mean of slice $k$. We then estimate
$M_{xy}$ as the covariance matrix of the $u_k$.

In SIR, we estimate the coefficient vectors $b_j$ by solving the generalized
eigenvalue problem (GEP)

$$
M_{xy}b = \lambda S_xb,
$$

where $S_x = {\rm cov(X)}$ is the covariance matrix of $X$ (ignoring $Y$).
Solving the GEP identifies vectors $b$ such that $b^\prime E[X|y]$ has high
variance relative to $b^\prime X$.

## Correspondence Analysis

Correspondence analysis (CA) is an embedding approach that aims to represent
*chi-square distances* in the data space as Euclidean distances for
visualization. The chi-square metric is defined as follows: if
$X, Y \in {\mathbb R}^p$ are random vectors with mean $EX = EY = \mu$ and
covariance ${\rm cov}(X) = {\rm cov}(Y) = {\rm diag}(\mu)$, then the (squared)
chi-square distance from $X$ to the mean is
$(X-\mu)^T{\rm diag}(\mu)^{-1}(X-\mu)$, and the squared chi-square distance
between $X$ and $Y$ is $(X-Y)^T{\rm diag}(\mu)^{-1}(X-Y)$. As discussed
further below, correspondence analysis arises most naturally when working with
nominal (categorical) variables, but there are some situations where it makes
sense to apply CA when the data are not nominal.

The motivation for embedding using the chi-square metric this is that in many
settings chi-square distances may best represent information in the data,
while Euclidean distances and angles are arguably the best approach for
producing visualizations for human interpretation. As discussed further below,
the chi-square metric is especially appropriate for understanding
relationships among the objects or observations (the rows of the data matrix),
which often can be viewed as an independent and identically distributed
collection of draws from a multivariate distribution.

### Mean/variance relationships

To better understand the motivation for using chi-square distances to measure
distances in the data space, let $X \in {\mathbb R}^p$ be a random vector with
mean $\mu \ge 0$ and covariance matrix $\Sigma$. In some cases, $\mu$ and
$\Sigma$ are unrelated (i.e. knowing $\mu$ places no constraints on $\Sigma$,
and vice-versa). On the other hand, in many settings it is plausible that
$\mu$ and $\Sigma$ are related via a *mean-variance relationship*. A common
form of mean-variance relationship is ${\rm diag}(\Sigma) \propto \mu$. A
particular example is the Poisson distribution where
${\rm diag}(\Sigma) = \mu$, so the constant of proportionality is $1$. In a
broader class of settings we may have *over-dispersion* or *under-dispersion*,
meaning that $\Sigma_{ii} = c\cdot \mu_i$, where $c>1$ or $c<1$ for over and
under-dispersion, respectively.

When comparing two random vectors $X$, $Y$ sharing covariance matrix $\Sigma$,
it is reasonable to use the
[Mahalanobis distance](https://en.wikipedia.org/wiki/Mahalanobis_distance)

$$
d(X, Y)^2 = (X-Y)^T\Sigma^{-1}(X-Y).
$$

The chi-square distance is a special case of the Mahalanobis distance in which
$\Sigma = {\rm diag}(\mu)$.

### Goals of MCA

Suppose we have $n$ observations on $p$ variables, and the data are
represented in an $n\times p$ matrix $X$ whose rows are the cases
(observations) and columns are the variables. Correspondence analysis can be
applied when each $X_{ij} \ge 0$, and where it makes sense to compare any two
rows or any two columns of $X$ using chi-square distance. Let $P \equiv X/N$,
where $N = \sum_{ij} X_{ij}$. We introduce the following notation: let
$P_{i,:}$ denote row $i$ of $P$, let $r \equiv P\cdot 1_p$ denote the row sums
of $P$, and let $c = P^T\cdot 1_n$ denote the column sums of $P$. Also, let
$W_r = {\rm diag}(r)\in {\mathbb R}^{n\times n}$ and
$W_c = {\rm diag}(c) \in {\mathbb R}^{p\times p}$. Then
$P^r \equiv W_r^{-1}\cdot P$ are the *row profiles* of $P$, which are simply
the rows of $P$ (or of $X$) normalized by their sum. Analogously,
$P^c \equiv P\cdot W_c^{-1}$ are the *column profiles* of $P$ (or of $X$).

The primary goal of MCA is to embed the row profiles $P^r$ into *row scores*
$F$ and the column profiles $P^c$ into *column scores* $G$, where $F$ is an
$n\times p$ array and $G$ is a $p\times p$ array, and both embeddings respect
chi-square distances.

In many cases, the row sums will be approximately constant, so
$W_r \propto 1_n$ and this matrix can effectively be ignored. However it is
very unlikely that $W_c \propto 1_p$, and in fact the inhomogeneity in $c$ is
the key feature in the data that motivates use of MCA.

Our goals are as follows:

- For any $1 \le i, j \le n$, the Euclidean distance from $F_{i,:}$ to
  $F_{j,:}$ is equal to the chi-square distance from $P^r_{i,:}$ to
  $P^r_{j,:}$. Also, for any $1 \le i,j \le p$ the Euclidean distance from
  $G_{i,:}$ to $G_{j,:}$ is equal to the chi-square distance from $P^c_{:,i}$
  to $P^c_{:,j}$. Thus, $F$ provides an embedding of the rows of $P^r$ and $G$
  provides an embedding of the columns of $P^c$. These embeddings map
  chi-square distances in the data space (specifically, the space of row
  profiles or column profiles) to Euclidean distances in the embedded space.

- The columns of $F$ and $G$ are ordered in terms of importance. Specifically,
  if we select $1 \le q \le p$ then the Euclidean distance from $F_{i,1:q}$ to
  $F_{j,1:q}$ is the best possible $q$-dimensional approximation to the
  chi-square distance from $P_{i,:}$ to $P_{j,:}$. Note that if $q=p$ then the
  approximation becomes exact, but for $q < p$ the approximation is inexact.

### Derivation of the algorithm

The array $P - rc^T$ is a doubly-centered array of residuals, since
$1_n^T(P - rc^T) = 0_p$, and $(P - rc^T)1_p = 0_n$. We can further standardize
these residuals by scaling the rows and columns by their respective standard
deviations (taking the variance to be proportional to the mean). This gives us
the doubly-standardized array of residuals

$$
W_r^{-1/2}(P - rc^T)W_c^{-1/2}.
$$

We next take the singular value decomposition of these doubly-standardized
residuals:

$$
W_r^{-1/2}(P - rc^T)W_c^{-1/2} = USV^T.
$$

Now let $F = W_r^{-1/2}US$ and $G = W_c^{-1/2}VS$. We show that this
specification of $F$ and $G$ satisfies the conditions stated above.
Specifically, we will show that the rows of $F$ embed the objects (the rows of
$X$), converting chi-square distances among the rows of $X$ to Euclidean
distances among the rows of $F$. If we make a scatterplot of the points
$(F_{i1}, F_{i2})$, then the Euclidean distances among these points
approximate the chi-square distances among the rows of $X$, providing a
meaningful Euclidean embeding of the rows of $X$.

First, note that since $V$ is orthogonal

$$
\lVert F_{i,:} - F_{j,:} \rVert = \lVert (F_{i,:} - F_{j,:})V^T \rVert.
$$

Therefore,

$$
\lVert F_{i,:} - F_{j,:} \rVert^2 =
\lVert r_i^{-1/2}U_{i,:}S - r_j^{-1/2}U_{j,:}S \rVert^2 =
\lVert r_i^{-1}(P_{i,:} - r_ic^T)W_c^{-1/2} - r_j^{-1}(P_{j,:} - r_jc^T)W_c^{-1/2}\rVert^2 =
$$

$$
\lVert r_i^{-1}P_{i,:}W_c^{-1/2} - r_j^{-1}P_{j,:}W_c^{-1/2} \rVert^2 =
(P_{i,:}/r_i - P_{j,:}/r_j)^TW_c^{-1}(P_{i,:}/r_i - P_{j,:}/r_j) =
$$

$$
(P^r_i - P^r_j)^TW_c^{-1}(P^r_i - P^r_j) = n^{-1}(P^r_i - P^r_j)^T(W_c/n)^{-1}(P^r_i - P^r_j).
$$

Suppose for simplicity that $W_r \propto I_n$. Since
$W_c/n = {\rm diag}(c/n)$, where $c/n$ is the mean of the rows of $P^r$, it
follows that $\lVert F_{i,:} - F_{j,:} \rVert$ is $n^{-1}$ times the squared
chi-square distance between $P^r_{i,:}$ and $P^r_{j,:}$. Thus, the rows of $F$
embed the rows of $P^r$ as desired. Applying the same argument to $X^T$ shows
that the rows of $G$ embed the columns of $P^c$, also reflecting chi-square
distances.

### Correspondence analysis and Multiple Correspondence analysis for nominal data

One common application of correspondence analysis (CA) arises when analyzing
datasets in which all variables are nominal. First, suppose we have a single
nominal variable and code it using an *indicator matrix*. That is, we define
$X$ to be a matrix whose values are entirely $0$ and $1$, such that $X_{ij}=1$
if and only if the value of the nominal variable for case $i$ is equal to
level $j$. Correspondence analysis as defined above can be used to analyze
this indicator matrix, revealing how the objects and categories are related.
Note that in this case $r = 1_n$ and $W_r = I_n$.

An important extension of CA is *Multiple Correspondence Analysis*, in which
we have several nominal variables. In this case, we recode each nominal
variable with its own indicator matrix, and then concatenate these matrices
horizontally. If there are $p_j$ levels for variable $j$, and we set
$p = \sum_j p_j$, then the concatenated indicator matrix (the *Burt matrix*)
is $n\times p$. We then apply CA to this concatenated indicator matrix,
yielding insights into the relationships among the objects, and relationships
among levels of different variables.

Note that in the special case where MCA is applied to a collection of nominal
variables, the Burt matrix $X$ has the property that every row sums to $q$
(where $q$ is the number of nominal variables). Thus every row of $P = X/(nq)$
sums to $1/n$ and the column means of $P^c$ are also all equal to $1/n$, since

$$
\bar{P}^c = n^{-1}1_n^\prime PQ_c^{-1} = n^{-1}1_p^T.
$$

### Angles and magnitudes of category scores

In many applications of CA and MCA, the rows are viewed as an IID sample from
a distribution, and the variables of this distribution exhibit a Poisson
mean/variance relationship. This motivates embedding the row profiles $P^r$
derived from the data matrix $X$ using chi-square distances. However the
argument for embedding the columns (variables) using the chi-square metric is
less clear (although it is true that MCA achieves such an embedding). A more
interpretable embedding of the column profiles is based on covariances and
correlations rather than on chi-square distances. Here we show that MCA also
achieves a column embedding that can be interpreted in terms of covariances
and correlations.

The angle between two embedded category points is an easily intepreted visual
feature of the plot of category embeddings (in practice usually restricted to
the dominant two components). The angle $\theta$ between two vectors $v$, $w$
satisfies

$$
\langle v, w\rangle = \lVert v \rVert \cdot \lVert w \rVert \cdot \cos(\theta).
$$

Therefore, to understand the angles among the embedded categories, we should
consider the Gram matrix $GG^T$ containing inner products of every pair of
embedded category points.

The matrix $GG^T$ has the following relationship to the data in $P$:

$$
GG^T = W_c^{-1/2}VSSV^TW_c^{-1/2} = W_c^{-1}(P - rc^T)^TW_r^{-1}(P - rc^T)W_c^{-1}.
$$

In most applications of MCA, $W_r = qI_n$, where $q$ is the number of
variables. Further, as noted above $\bar{P}^c_{:,i} = 1/n$ for each $i$, and
$r = n^{-1}1_n$. Therefore, a single element of the Gram matrix has the form

$$
[GG^T]\_{ij} = q^{-1}(P^c_{:,i} - r)^T(P^c_{:,j} - r) = q^{-1}n\cdot {\rm cov}(P^c_{:,i}, P^c_{:,j}).
$$

The diagonal elements of the Gram matrix are proportional to the variances of
the column profiles:

$$
[GG^T]\_{ii} = q^{-1}n\cdot {\rm var}(P^c_{:,i}).
$$

It follows that the angle between two category profiles, say in columns $i$
and $j$, is proportional to the correlation coefficient between the column
profiles $P^c_i$ and $P^c_j$.

In addition, a category point lies further from the origin if the column
profiles are not close to being constant vectors. These are the column
profiles that most strongly influence the MCA fit.

Combining these observations, two column profiles corresponding to categories
of different variables that are both far from the origin, and that have a
small angle, correspond to substantially correlated indicator vectors.

Note that the interpretation of the category points in MCA usually focuses on
the relationships between pairs of categories for different variables.
Comparing two categories of the same variable is rarely interesting since
these indicators are by definition mutually exclusive.

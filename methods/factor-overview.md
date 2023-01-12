# Factor analyses and embeddings

A large class of powerful statistical methods considers data in which
many "objects" ("observations") are measured with respect to multiple variables.  At a
high level, these analyses usually aim to simultaneously understand
the relationships among the objects and the relationships among the
variables.

For example, we may have a sample of people (the "objects") with each
person measured in terms of their height, weight, age, and sex.  In
this example, two people are "similar" if they have similar values on
most or all of the variables, and two variables are "similar" if
knowing the value of just one of the variables for a particular object
allows one to predict the value of the other variable for that same
object.

## Embedding

Many methods of factor analysis can be viewed as a means to obtain
an *embedding*.  Embedding algorithms take input data vectors $X$ and transform them
into output data vectors $Z$.  Many embedding algorithms take the form
of a "dimension reduction", so that
$q \equiv {\rm dim}(Z) < {\rm dim}(X) \equiv p$.
However some embeddings preserve or even increase the dimension.

Embeddings can be used for exploratory analysis, especially in
visualizations, and can also be used to construct features for
prediction, as well as being used in formal statistical inference,
e.g. in hypothesis testing.

An embedding approach is linear if $Z = BX$ for a fixed but
data-dependent $q\times p$ matrix $B$.  Linear embedding algorithms
are simpler to devise and characterize than nonlinear embeddings.
Many modern embedding algorithms are nonlinear, exploiting the
potential to better capture complex structure.

Some embedding algorithms embed only the objects while other embedding
algorithms embed both the objects and the variables.  Embedding the
objects provides a reduced feature representation for each object that
can be passed on to additional analysis procedures, or interpreted
directly.  Embedding the variables provides a means to interpret the
relationships among the variables.

## Centering and standardization

There are many ways to pre-process the data prior to performing
a factor analysis.  Suppose that $X$ is a $n\times p$ matrix whose rows
are the observations or objects, and whose columns are the variables.

Some factor-type methods work with the covariance matrix of the variables,
which is a $p\times p$ positive semidefinite matrix.  Since covariances
by definition are derived from mean centered variables, it is common to
mean center the variables (columns) of $X$ prior to running a factor analysis.

Furthermore, in many cases we do not wish our results to depend on the
units in which each variable was measured, or we wish to explicitly
remove any influence of differing dispersions of the variables on our results.
This motivates standardizing the columns of $X$ prior to performing
a factor analysis, where *standardization* involves first mean
centering each column and then dividing each column by its standard
deviation.

While standardization is commonly performed, in some cases it is desirable for the variables
with more dispersion to have more influence on the results of a factor
analysis, and in these cases one may choose not to standardize the
variables.

In some datasets there is no clear distinction between
an observation and a variable, and it may be desirable to standardize
both the rows and the columns.  To achieve this, the following three
steps can be performed: (i) center the overall matrix around its
grand mean, (ii) center each row of the matrix with respect to the
row mean, (iii) center each column of the matrix with respect to
the column mean.  After these three steps, each row and each column
of the matrix will have mean zero and we can refer to the matrix as
having been *double centered*.

## Singular Value Decoposition

Many embedding methods make use of a matrix factorization known as the
*Singular Value Decomposition* (SVD).  The SVD is defined for any
$n\times p$ matrix $X$.  In most cases we want $n \ge p$, and if $n < p$
we would take the SVD of $X^T$ instead of $X$.  When $n\ge p$, we
decompose $X$ as $X = USV^T$, where $U$ is $n\times p$, $S$ is
$p\times p$, and $V$ is $p\times p$.  The matrices $U$ and $V$ are
orthogonal so that $U^TU = I_p$, $V^TV = I_p$, and $S$ is diagonal
with $S_{11} \ge S_{22} \ge \cdots \ge S_{pp}$.  The values on the
diagonal of $S$ are the *singular values* of $S$, and the SVD is
unique except when there are ties among the singular values.

One use of the SVD is to obtain a low rank approximation to a matrix
$X$.  Suppose we truncate the SVD using only the first $k$ components,
so that $\tilde{U}$ is the $n\times k$ matrix consisting of the
leading (left-most) $k$ columns of $U$, $\tilde{S}$ is the upper-left
$k\times k$ block of $S$, and $V$ is the $p\times k$ matrix consisting
of the leading $k$ columns of $V$.  In this case, the matrix
$\tilde{X} \equiv \tilde{U}\tilde{S}\tilde{V}^T$ is a rank $k$ matrix
(it has $k$ non-zero singular values).  Among all rank $k$ matrices,
$\tilde{X}$ is the closest matrix to $X$ in the *Frobenius norm*,
which is defined as

$$
{\rm Frob}(X)^2 = \|X\|_F^2 \equiv \sum_{ij}X_{ij}^2 = {\rm trace}(X^\prime X).
$$

Thus we have

$$
\tilde{X} = {\rm argmin}_{A: {\rm rank}(A) = k} \|A - X\|_F.
$$

### Analyzing a data matrix using the SVD

Suppose we have a data matrix $X$ and wish to understand its structure.  We
can begin my double centering the matrix $X$ as discussed above:

$$
X_{ij} = m + r_i + c_j + S_{ij}.
$$

Here, $S_{ij} \equiv X_{ij} - m - r_i - c_j$ are residuals.  Next we can
take the SVD of $S$, yielding $S = USV^T$, or equivalently

$$
S_{ij} = \sum_{k=1}^p S_{kk} U_{ik}V_{jk}.
$$

Although there are $p$ terms in the SVD of $S$, the first few terms may capture
most of the structure.  One important property that results from calculating
the SVD for a double-centered matrix is that $U_{\cdot k} = 0$ and $V_{\cdot k} = 0$.
That is, the columns of $U$ and $V$ are centered.  This column centering means
that the SVD captures "deviations from the mean" represented by the additive
model $m + r_i + c_j$.

## Principal Components Analysis

Suppose that $X$ is a $p$-dimensional random vector with mean $0$ and
covariance matrix $\Sigma$ (our focus here is not the mean, so if $X$
does not have mean zero we can replace it with $X-\mu$, where
$\mu=EX$).  Principal Components Analysis (PCA) seeks a linear
embedding of $X$ into a lower dimensional space of dimension $q < p$.
The standard PCA approach gives us an orthogonal matrix $B$ of
*loadings*, which can be used to produce *scores* denoted $Q$.

For a single vector $X$, the scores are obtained via the mapping $Q(X)
= B^TX$.  For a data matrix $Z$ whose rows are independent and
identically distributed (IID) copies of the random vector $X$, the
scores can be obtained via the mapping $Q = ZB$, where each row of $Q$
contains the scores for the corresponding row of $Z$.

PCA can be viewed in terms of linear compression and decompression of
the variables in $X$.  Let

$$
X \rightarrow B_{:,1:q}^TX \equiv Q(X)
$$

denote the compression that reduces the data from $p$ dimensions to
$q$ dimensions, where $B_{:,1:q}$ is the $p\times q$ matrix consisting
of the leading $q$ columns of $B$.  We can now decompress the data as
follows:

$$ Q \rightarrow B_{:,1:q}Q \equiv \hat{X}.  $$

This represents a two-step process of first reducing the dimension,
then predicting the original data using the reduced data.  Among all
possible loading matrices $B$, the PCA loading matrix loses the least
information in that it minimizes the expected value of $\|X -
\hat{X}\|$.

The loading matrix $B$ used in PCA is the eigenvector matrix of
$\Sigma$.  Specifically, we can write $\Sigma = B\Lambda B^T$, where
$B$ is an orthogonal matrix and $\Lambda$ is a diagonal matrix with
$\Lambda_{11} \ge \Lambda_{22} \ge \cdots \ge \Lambda_{pp} > 0$.  This
is the spectral decomposition of $\Sigma$.  The columns of $B$ are
orthogonal in both the Euclidean metric and in the metric of $\Sigma$,
that is, $B^TB = I_p$ and $B^T\Sigma B = I_p$.  As a result, the
scores $Q \equiv B^TX$ are uncorrelated, ${\rm cov}(Q) = I_q$.

Next we consider how PCA can be carried out with a sample of data,
rather than in a population.  Given a $n\times p$ matrix of data $Z$
whose rows are iid copies of the random vector $X$, we can estimate
the covariance matrix $\Sigma$ by column centering $Z$ to produce
$Z_c \equiv X - 1_n\bar{X}^T$, where $\bar{X} \in {\cal R}^p$ is the vector
of column-wise means of $X$, Then set
$\hat{\Sigma} = Z_c^TZ_c/n$. Letting $B$ denote the eigenvectors of $\hat{\Sigma}$,
the scores have the form $Q = Z_cB$.

Since the eigenvalues $\Lambda_{ii}$ are non-increasing, the leading
columns of $Q$ contain the greatest fraction of information about $Z$.
Thus, visualizations (e.g. scatterplots) of the first two columns of
$Z$ best reflect the relationships among the rows of $Z$ (compared to
any other scatterplot formed from linear scores).

PCA can also be carried out using the SVD of the column-centered data
matrix $Z_c$.  Write $Z_c = USV^T$ as discussed above and note that

$$
\hat{\Sigma} = Z_c^T/Z_c/n = VS^2V^T/n.
$$

Thus, $\hat{\Sigma}V = VS^2/n$, so $V$ contains the eigenvectors of
$\hat{\Sigma}$ and the corresponding eigenvalues are in ${\rm
diag}(S^2)/n$.

### Principal Components Regression

Principal Components Analysis is a method for multivariate data analysis
that can be used to understand relationships in multivariate data when
there is no single variable that can be seen as the "response" relative
to the other variables.  PCA can also be used to create covariates
for use in a regression analysis.  This method is called "Principal
Components Regression" (PCR).  The motivation behind PCR is that if we
have a large number of covariates and do not wish to explicitly
include all of them in a regression model, we can use PCA to
reduce the variables to a smaller set of scores that capture most
of the variation in the original variables, and then use the
reduced variables (scores) as covariates in our regression.

More formally, suppose that $X \in {\cal R}^{n\times p}$ is our
regression design matrix, for a regression with $n$ observations
(cases) and $p$ covariates.  The response variable for this
regression is a vector $Y \in {\cal R}^n$, but this is not used
in the first few steps of PCR.  First, the mean should be subtracted
from each column of $X$, and we write $X = USV^T$ using the SVD.
Next we truncate this representation so that $\tilde{U} \in {\cal R}^{n\times q}$,
contains the first $q$ columns of $U$.  As discussed above, these
are the most important columns of $U$ for explaining the variation
in $X$.

We use $\tilde{U}$ instead of $X$ as covariates in our regression, obtaining a
linear predictor $\tilde{U}\hat{\beta}$, where $\hat{\beta}$
are coefficients obtained using a regression procedure,
e.g. least squares or a GLM.  To relate the coefficients
$\hat{\beta}$ back to the original covariates, note that
$U = XVS^{-1}$, and $\tilde{U} = X(VS^{-1})_{:,1:q}$.
Therefore

$$
\tilde{U}\hat{\beta} = X(VS^{-1})_{:,1:q}\hat{\beta},
$$

and we can define

$$
\hat{\gamma} = (VS^{-1})_{:,1:q}\hat{\beta}.
$$

The coefficients $\hat{\gamma}$ align with the original covariates,
while the coefficients $\hat{\beta}$ align with the principal
component scores.  Both determine the same fitted values and
linear predictor.

PCR can be very effective, but may or may not peform well
depending on the circumstances.  In PCR, instead
of allowing the coefficients $\gamma$ to take on any
value in ${\cal R}^p$, we are constraining $\gamma$
to lie in the columnspace of $(VS^{-1})_{:,1:q}$. PCR
works well if the projections of $X$ that explain the
most variation in $X$ also explain the most variation
in $Y$.  This is often, but not always the case.

## Canonical Correlation Analysis

Suppose we have two sets of variables and a collection of objects that
are assessed with respect to all variables in both sets. CCA aims to understand
the relationships between the two sets of variables.  This goal differs from
the goal of PCA, which could be used to understand the
relationships within each set of variables separately.
Specifically, we can search for linear combinations of each set of variables
that are maximally correlated with each other.

Let $X \in {\cal R}^{n\times p}$ and $Y \in {\cal R}^{n\times q}$ denote
the data for the two sets of variables.  In these data, we have measurements
on $n$ objects (the rows of $X$ and $Y$), a set of $p$ variables in the columns of $X$,
and a set of $q$ variables in the columns of $Y$.  Row $i$ of $X$ and row $i$ of $Y$
correspond to the same object.

Given coefficient vectors $a\in {\cal R}^p$ and $b\in {\cal R}^q$, we can
form linear predictors $Xa$ and $Yb$, both of which are $n$-dimensional
vectors.  The goal of CCA is to find $a$ and $b$ to maximize the correlation
coefficient between $Xa$ and $Yb$.  Inspecting plots of these "canonical coefficients"
can yield insight into the relationship between the variables in $X$ and the variables
in $Y$.

If the dimension is high, CCA can overfit the data.  A procedure analogous to
principal components regression can be used to address this.  We first reduce $X$ and $Y$ to lower
dimensions using PCA, then we fit the CCA using these reduced data.  The coefficients
from the reduced PCA can then be re-expressed in the original coordinates for interpretation.

## Dimension reduction regression

Dimension reduction regression (DR) is a flexible approach to regression analysis that
can be seen as extending PCA to the setting of regression.  In a DR analysis, we have
a matrix $X \in {\cal R}^{n\times p}$ containing data on the explanatory variables,
and we also have a vector $Y \in {\cal R}^n$
containing values of a response variable.  Like PCA, the goal is to find *factors* or *components*
among the $X$ variables, but in the case of DR the goal is for these factors to predict $Y$,
not to predict $X$ itself.  One way to view DR is as a means to "steer" PCA toward variates
that explain the variation in $Y$, rather than seeking variates that explain the variation
in $X$.

In some cases, the variates in $X$ that explain $X$ and the variates in $X$ that
explain $Y$ can be quite similar, but in other cases they dramatically differ.  This is
one of the main critiques of principal component regression (PCR), which uses the PCs as
explanatory variables and therefore can only find regression relationships that happen
to coincide with the PCs.

DR can be viewed as fitting the multi-index regression model

$$
E[Y|X=x] = f(b_1^\prime x, \ldots, b_q^\prime x) + \epsilon,
$$

where the $b_j \in {\cal R}^p$ are coefficient vectors defining the "indices" in
$X$ that predict $Y$.  The span of $b_1, \ldots, b_q$ is known as the "dimension reduction"
subpsace for the regression function $E[Y|X=x]$.

An appealing feature of the DR approach is that it incorporates
a link function $f$ but this function does not need to be known.  Thus, DR provides
a means to capture a wide range of nonlinear and non-additive regression relationships.
The parameter $q$ above determines the dimension of the dimension reduction subspace, and
is selected based on the data.  If $q=1$ we have a single-index model, like in a
GLM, except that here the link function does not need to be pre-specified.  As $q$ grows,
the model becomes more complex and therefore may underperform due to high variance
and overfitting.  DR is most effective when relatively small values of $q$ can
be used, thereby compressing the regression structure into a few variates.

It is advantageous that the link function $f$ does not need to be known while estimating
the coefficients $b_j$.  However later in the analysis we will probably want to estimate
$f$, and commonly a nonparametric regression method like loess can be used for
this purpose.

One of the simplest and most widely-used approaches to dimension reduction regression
is known as Sliced Inverse Regression (SIR).  Recall that the loading vectors of PCA
are the eigenvectors of ${\rm cov}(X)$.  In SIR, we wish to steer the PCA loadings
toward directions that are relevant for predicting $Y$.  To do this we consider
$M_{xy} \equiv {\rm cov}E[X|Y]$, which is also a type of covariance matrix for $X$, but which
suppresses the variation in $X$ that is irrelevant for $Y$.  This is accomplished
by replacing $X$ in ${\rm cov}(X)$ with $E[X|Y]$, which suppresses the variation
in $X$ that occurs when $Y$ is fixed. The matrix $M_{xy}$ is estimated by sorting
the data $\{(x_i, y_i)\}$ by increasing values of $y$, dividing this sorted sequence
into "slices" (blocks), and averaging the values of $x_i$ within each slice.  Let
$u_k \in {\cal R}^d$
denote the mean of slice $k$.  We then estimate $M_{xy}$ as the covariance matrix of the
$u_k$.

In SIR, we estimate the coefficient vectors $b_j$ by solving the generalized eigenvalue problem
(GEP)

$$
M_{xy}b = \lambda S_xb,
$$

where $S_x = {\rm cov(X)}$ is the covariance matrix of $X$ (ignoring $Y$).  Solving the GEP
identifies vectors $b$ such that $b^\prime E[X|y]$ has high variance relative to $b^\prime X$.

## Correspondence Analysis

Correspondence analysis is an embedding approach that aims to
represent *chi-square distances* in the data space as Euclidean
distances for visualization.  The motivation for doing this is that in
many settings chi-square distances may be the best approach for
summarizing the information in the data, while Euclidean distances are
arguably the best approach for producing visualizations for human
interpretation.

Let $X \in {\cal R}^p$ be a random vector with mean $\mu \ge 0$ and
covariance matrix $\Sigma$.  In some cases, it is reasonable to view
$\mu$ and $\Sigma$ as unrelated (i.e. knowing $\mu$ places no
constraints on $\Sigma$, and vice-versa).  On the other hand, in many
settings it is plausible that $\mu$ and $\Sigma$ are related in that
${\rm diag}(\Sigma) \propto \mu$.  Specifically in a Poisson
distribution ${\rm diag}(\Sigma) = \mu$, but in a broader class of
settings we may have over-dispersion or under-dispersion, meaning that
$\Sigma_{ii} = c\cdot \mu_i$, where $c>1$ or $c<1$ for over and
under-dispersion, respectively.

### Mean/variance relationships

In any setting where the variance is proportional to the mean, it is
reasonable to compare vectors using chi-square distances.
Specifically the (squared) chi-square distance from $X$ to the mean is
$(X-\mu)^T{\rm diag}(\mu)^{-1}(X-\mu)$, and the squared chi-square
distance between two random vectors $X$ and $Y$ having the same mean
$\mu$ is $(X-Y)^T{\rm diag}(\mu)^{-1}(X-Y)$.  If $\Sigma$ is a
diagonal matrix, the chi-square distance is also the *Mahalanobis
distance*, which is arguably the proper way to measure distances among
vectors whose components have differing variances.

Suppose we have $n$ observations on $p$ variables, and the data are
represented in an $n\times p$ matrix $X$ whose rows are the cases
(observations) and columns are the variables.  Correspondence analysis
can be applied when each $X_{ij} \ge 0$, and where it makes sense to
compare any two rows or any two columns of $X$ using chi-square
distance.  Let $P \equiv X/N$, where $N = \sum_{ij} X_{ij}$.  The goal
is to transform $P$ into *row scores* $F$ and *column scores* $G$,
where $F$ is an $n\times p$ array and $G$ is a $p\times p$ array.

### Goals of MCA

We introduce the following notation: let $P_{i,:}$, $F_{i,:}$, and
$G_{i,:}$ denote row $i$ of the arrays $P$, $F$, and $G$ respectively,
let $r \equiv P\cdot 1_p$ (the row sums of $P$) and $c = P^T\cdot 1_n$
(the column sums of $P$), and let $W_r = {\rm diag}(r)\in {\cal
R}^{n\times n}$ and $W_c = {\rm diag}(c) \in {\cal R}^{p\times p}$.
Then let $P^r \equiv W_r^{-1}\cdot P$ denote the *row profiles* of
$P$, which are simply the rows of $P$ (or of $X$) normalized by their
sum.  Analogously, let $P^c \equiv P\cdot W_c^{-1}$ denote the *column
profiles* of $P$ (or of $X$).

Our goals are as follows:

* For any $1 \le i, j \le n$, the Euclidean distance from $F_{i,:}$ to
$F_{j,:}$ is equal to the chi-square distance from $P^r_{i,:}$ to
$P^r_{j,:}$.  Also, for any $1 \le i,j \le p$ the Euclidean distance
from $G_{i,:}$ to $G_{j,:}$ is equal to the chi-square distance from
$P^c_{:,i}$ to $P^c_{:,j}$.  Thus, $F$ provides an embedding of the
rows of $P^r$ and $G$ provides an embedding of the columns of $P^c$.

* The columns of $F$ and $G$ are ordered in terms of importance.
Specifically, if we select $1 \le q \le p$ then the Euclidean distance
from $F_{i,1:q}$ to $F_{j,1:q}$ is approximately equal to the
chi-square distance from $P_{i,:}$ to $P_{j,:}$.  Note that if $q=p$
then the approximation becomes exact, but for $q < p$ the approximation
is inexact.

### Derivation of the algorithm

The array $P - rc^T$ is a doubly-centered array of residuals,
in that $1_n^T(P - rc^T) = 0_p$, and $(P - rc^T)1_p = 0_n$.
We can further standardize these residuals by scaling the rows
and columns by their respective standard deviations (taking the
variance to be proportional to the mean).  This gives us the
doubly-standardized array of residuals

$$
W_r^{-1/2}(P - rc^T)W_c^{-1/2}.
$$

We next take the singular value decomposition of these doubly-standardized
residuals:

$$
W_r^{-1/2}(P - rc^T)W_c^{-1/2} = USV^T.
$$

Now let $F = W_r^{-1/2}US$ and $G = W_c^{-1/2}VS$.  We now show that
this specification of $F$ and $G$ satisfies the conditions stated
above.

First, note that since $V$ is orthogonal

$$
\|F_{i,:} - F_{j,:}\| = \|(F_{i,:} - F_{j,:})V^T\|.
$$

Therefore,

$$
\|F_{i,:} - F_{j,:}\|^2 =
\|r_i^{-1/2}U_{i,:}S - r_j^{-1/2}U_{j,:}S\|^2 =
\|r_i^{-1}(P_{i,:} - r_ic^T)W_c^{-1/2} - r_j^{-1}(P_{j,:} - r_jc^T)W_c^{-1/2}\|^2 =
$$

$$
\|r_i^{-1}P_{i,:}W_c^{-1/2} - r_j^{-1}P_{j,:}W_c^{-1/2}\|^2 =
(P_{i,:}/r_i - P_{j,:}/r_j)^TW_c^{-1}(P_{i,:}/r_i - P_{j,:}/r_j).
$$

Since $W_c = {\rm diag}(\hat{\mu})$, where $\hat{\mu}$ is an estimate
of $\mu$, it follows that $\\|F_{i,:} - F_{j,:} \\|$ is an estimate of
the chi-square distance between $P^r_{i,:}$ and $P^r_{j,:}$.  Thus,
the rows of $F$ embed the rows of $P^r$ as desired.  Applying the same
argument to $X^T$ shows that the rows of $G$ embed the columns of
$P^c$.

### Correspondence analysis and Multiple Correspondence analysis for nominal data

One common application of correspondence analysis (CA) arises when
analyzing datasets in which all variables are nominal.  First, suppose
we have a single nominal variable and code it using an *indicator
matrix*.  That is, we define $X$ to be a matrix whose values are
entirely $0$ and $1$, such that $X_{ij}=1$ if and only if the value of
the nominal variable for case $i$ is equal to level $j$.
Correspondence analysis as defined above can be used to analyze this
indicator matrix, revealing how the objects and categories are
related.

An important extension of CA is *Multiple Correspondence Analysis*, in
which we have several nominal variables.  In this case, we recode each
nominal variable with its own indicator matrix, and then concatenate
these matrices horizontally.  If there are $p_j$ levels for variable
$j$, and we set $p = \sum_j p_j$, then the concatenated indicator
matrix is $n\times p$.  We then apply CA to this concatenated
indicator matrix, yielding insights into the relationships among the
objects, the relationships between different levels of a single
variable, and relationships among levels of different variables.

Suppose we are mainly interested in the relationships among the
variables, and we use the scatterplot of variable scores to address
this.  Every category of every variable will have a point in this
scatterplot.  The main interest is in learning about how a category of
one variable relates to a category of another (different) variable.
As noted above, the distances between these points can be interpreted
in terms of chi-square distances, but it is also informative to
consider the magnitudes of, and angles among the vectors from the
origin to each point defined by the category scores.

### Angles and magnitudes of category scores

Above we discussed how the Euclidean distances between the rows of $G$
correspond to the chi-square distances between the columns of $P^c$.
Another insight into the geometry of MCA can be obtained my
considering the angles among the category scores (or specifically, the
angles between the vectors from the origin to each embedded category
point).  The angle between two vectors is the dot product between the
vectors, divided by the product of the norms of the vectors.
Therefore, to understand the angles among the embedded categories, we
should consider the Gram matrix $GG^T$ containing inner products of
every category score with every other category score.

The matrix $GG^T$ has the following relationship to the data in $P$:

$$
GG^T = W_c^{-1/2}VSSV^TW_c^{-1/2} = W_c^{-1}(P - rc^T)^TW_r^{-1}(P - rc^T)W_c^{-1}.
$$

In most applications of MCA, $W_r = n^{-1}I_n$, so a single element of
this matrix has the simpler form

$$
[GG^T]_{ij} = n(P^c_{:,i} - r)^T(P^c_{:,j} - r).
$$

Note that $\bar{P}^c_{:,i} = 1/n$ for each $i$, and $r \equiv 1/n$,
therefore, $P^c_{:,i} - r$ has average value zero and can be
interpreted as a vector of deviations (residuals) of each variable
from its mean.  The dot product between two such vectors is a type of
covariance (up to a scale factor) that captures whether two indicators
tend to co-occur within observations.

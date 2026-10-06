#import "@preview/ilm:2.1.1": *

#show link: set text(fill: blue)

#set text(lang: "en")

#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

#show: ilm.with(
  title: [Basic regression analysis],
  authors: "Kerby Shedden",
  figure-index: (enabled: true),
  table-index: (enabled: true),
  listing-index: (enabled: true),
  chapter-pagebreak: false
)


= Introduction

This document covers several topics that are relevant for many kinds of regression analysis where the focus is on the conditional location (e.g. mean or median) of a univariate response.  We do not focus here on more advanced topics including mean/variance relationships (best handled using generalized linear models), non-independent samples (handled using mixed effects regression or generalized estimating equations regression), or multivariate responses.

Throughout the discussion below, we mostly avoid speaking in terms of _generative models_.  This means that you will usually not see expressions such as $y = beta' x + epsilon$.  Instead, we will limit ourselves to explaining particular numerical characteristics of interest through the covariates.  For example, we may be interested in models for the conditional mean, conditional variance, or conditional median, which can be expressed $E[Y|X=x] = g(x)$, $"Var"[Y|X=x] = g(x)$, or $Q_(0.5)[Y|X=x] = g(x)$, respectively.

If we specify a regression model, e.g. for the conditional mean, in the form $E[Y|X=x] = g(x)$, then some of the main questions are: how to specify the range of possibilities for $g$, how to estimate $g$, and how to interpret the estimated function $hat(g)$?  Placing meaningful constraints on the range of possible values of $g$ is an important strategy.  These constraints often involve the notions of _additivity_ and _linearity_.  These terms will be discussed in more detail below, but the basic idea is that $g$ is additive with respect to $x_j$ and $x_k$ ($j != k$) if the mixed partial derivatives $frac(partial^2 g, partial x_j partial x_k, style: "horizontal")$ are identically 0, and $g$ is linear in $x_j$ when $frac(partial^2 g, partial x_j^2, style: "horizontal") equiv 0$.

= Basis functions

Many approaches to regression analysis relate the expected value of the response variable $y$ to a _linear predictor_ $beta' x = beta_1 x_1 + dots.h.c + beta_p x_p$, formed from the covariates $x in RR^p$ using coefficients $beta in RR^p$. Examples include the linear mean structure model $E[y|x] = x' beta$ and the _single index model_ with _link function_ $g$, $E[y|x] = g^(-1)(x' beta)$.

The linear mean structure model is linear in two senses -- the conditional mean of $y$ given $x$ is linear in $x$ for fixed $beta$, and it is linear in $beta$ for fixed $x$. These two forms of linearity have very different implications.

- Linearity in $beta$ for fixed $x$ makes it much easier to characterize the theoretical properties of the estimation process. Linearity in $beta$ generally also makes it easier to develop algorithms to compute the estimates.

- Linearity in $x$ for fixed $beta$ is sometimes cited as a weakness of this type of model. People incorrectly argue that models with this property are only suitable for describing systems that behave linearly, and since most natural and social processes are not linear, such models are sometimes claimed to have limited utility.

The (apparent) linearity of the mean structure model in the covariate vector $x$ is easily overcome. For a quantitative covariate $x$, it is possible to include both $x$ and $x^2$ as covariates in a "linear" model, leading to a linear predictor of the form $beta_1 x + beta_2 x^2$. This retains the benefits of linear estimation, while allowing the model for the conditional mean function to be non-linear in the covariates.

Including powers of covariates (like $x^2$) as regressors is a method known as _polynomial regression_. It may be the earliest example of the general technique of utilizing _basis functions_ to incorporate nonlinearity into regression analyses. A family of (univariate) basis functions is a collection of functions $g_1, g_2, ...$, each from $RR -> RR$, such that we can include $g_1(x), g_2(x), ...$ as covariates in a model in place of $x$. This allows the fitted mean function to take on any form that can be represented as a linear combination

$
beta_1 g_1(x) + beta_2 g_2(x) + dots.h.c + beta_p g_(p)(x).
$

The parameters $beta_j$ can be estimated using least squares, or other approaches such as penalized least squares (lasso/ridge/elastic net),
#link("https://en.wikipedia.org/wiki/Least_absolute_deviations")[least absolute deviations], or other forms of #link("https://en.wikipedia.org/wiki/M-estimator")[M-estimation]. Using a large collection of basis functions allows a wide range of non-linear forms to be represented. Basis functions thus allow non-linear mean structures to be fit to data using linear estimation techniques.  As discussed below, basis functions alone do not provide an optimal solution for this type of nonlinear regression analysis, but when combined with regularization they arguably do.

While the basis function approach is very powerful, some families of basis functions have undesirable properties. For example, polynomial basis functions can exhibit poor scaling and high colinearity (depending on the range of the data). Also, polynomial basis functions are not "local", meaning that when using polynomial basis functions, the fitted value at a point $x$ can depend
on data values $(x_i, y_i)$ where $x_i$ is far from $x$.

_Local basis functions_ are families of functions such that each element in the family has limited support (that is, the functions are exactly zero outside of a fairly small interval). One of the most useful forms of local basis functions is #link("https://en.wikipedia.org/wiki/Spline_(mathematics)")[splines], specifically, _polynomial splines_. Roughly speaking, a polynomial spline is a continuous and somewhat smooth function that has bounded support, and is a piecewise polynomial. There are several different constructions, and we will not derive the mathematical form of a polynomial spline here in detail. Polynomial splines are "somewhat smooth" in that they have a given finite number of continuous derivatives, e.g. the function, and its first and second derivative are continuous, but the third and higher-order derivative are not.

Other useful families of basis functions are #link("https://en.wikipedia.org/wiki/Wavelet")[wavelets], #link("https://en.wikipedia.org/wiki/Fourier_series")[Fourier series], #link("https://en.wikipedia.org/wiki/Radial_basis_function")[radial basis functions], and higher-dimensional basis functions formed via tensor products of univariate basis functions such as splines.

When working with basis functions, it is important to remember that terms derived from the same parent variable cannot vary independently of each other. For example, age determines $"age"^2$ and vice-versa. This means that in general, the coefficients of variables in a regression model that uses basis functions may be difficult to interpret (e.g. the "effect" of age is represented through the coefficients for $"age"$ and for $"age"^2$). There are many ways to resolve this, especially using plots. For example, if we have a model relating BMI to age and sex, and we use basis functions to capture a non-linear role for age, we can make a plot showing the fitted values of $E["BMI" | "age", "sex"]$ plotted against age, for each sex.

Using splines or other families of basis functions is a very powerful technique because it allows familiar estimation methods to be used in a much broader range of settings, simply by augmenting the regression design matrix with additional columns.

== Generalized Additive Models

Additive models are a powerful class of regression methods that combines the use of basis functions with smoothness penalties. Many of these methods broadly can be considered forms of _generalized additive modeling_ (GAM). The basic idea of a GAM is that we begin with the mean structure model

$
g^(-1)(E[y|x]) = beta_0 + beta_1 g_1(x) + beta_2 g_2(x) + dots.h.c + beta_p g_(p)(x).
$

This is a type of _semi-parametric model_, since the $g_j$ are functions that lie in large function spaces whose dimension can grow with the sample size.

For the moment suppose that there is only one covariate $x$. In a GAM, we impose a smoothing penalty, often based on the second derivative of the fitted regression function, such as the integrated squared second derivative

$
sum_i (beta_1 g_1^('')(x_i) + dots.h.c + beta_p g_p^('')(x_i))^2.
$

The quantity above is larger when the fitted regression function $sum_j beta_j g_(j)(x)$ is less smooth (and is zero when the fitted regression is linear in $x$). Note that the penalty is quadratic in $beta$ and therefore when using this expression to penalize the usual regression sum of squares, the calculations of the estimator and standard errors are straightforward.

Considering now the setting with more than one covariate, a true generalized additive model is additive in the sense that we model the conditional mean in the form

$
g^(-1)(E[y|x_1, ..., x_p]) = sum_j sum_k beta_(j k) g_(j k)(x_j),
$

and use a smoothing penalty of the form

$
sum_j (sum_k beta_(j k) g_(j k)^('')(x_j))^2.
$

A true GAM is additive in the sense that $g^(-1)(E[y | x_1, ..., x_p]) = sum_j h_(j)(x_j)$, where each $h_j$ is represented in terms of basis functions. Such a model does not permit any interactions among the covariates.  The above penalty can be expressed as a quadratic form $beta' P_g beta$, where $P_g$ is a $p times p$ matrix that does not depend on $beta$ or on $y$.

In practice, strict additivity is often too limiting, so most GAM software supports the inclusion of selected pairwise or higher order interactions, in which case the model is no longer additive, and the "automatic" nature of a pure GAM is lost. Defining tractable smoothing penalties for non-additive models is challenging and remains an area of research.

Another form of model that is often encountered is a _partial linear model_.  In this approach we partition the covariates into $x in RR^p$ and $z in RR^q$, and the model has the form

$
g^(-1)(E[y|x_1, ..., x_p, z]) = sum_(j=1)^p h_(j)(x_j) + theta' z.
$

= Transformations

A _transformation_ in statistics refers to any setting in which a function is applied to a variable being analyzed. It is easy to apply transformations, and there are many principled reasons for transforming data.

In regression analysis, transformations can be applied to the dependent variable $y$, to one or more of the independent variables $x_j$, or to both the independent and dependent variables simultaneously. One reason for transforming the data is to induce it to fit into a given regression framework. For example, ordinary least squares (OLS) is most efficient when the conditional mean $E[y|x]$ is linear in $x$, and the conditional variance $"Var"[y|x]$ is constant ("homoscedasticity"). A transformation that achieves the latter is called a _variance stabilizing transformation_.  Sometimes, applying a transformation such as replacing $y$ with $log(y)$ will induce linearity of the mean structure and homoscedasticity of the variance structure. However achieving linearity of $E[y|x]$ and achieving homoscedasticity (constant conditional variance) cannot always be achieved with the same transformation.

Methods for automating the process of selecting a transformation in linear regression have been proposed, the most well-known being the Box-Cox method. However this only automates the process of selecting a transformation for the dependent variable $y$. Transforming covariates is usually a manual process of trial and error.

It is sometimes mistakenly believed that the dependent variable $y$ in a linear model should marginally follow a symmetric or (even stronger) a Gaussian distribution. In general however, the marginal distribution of $y$ is irrelevant in regression analysis. In a linear model, we might like the "unexplained variation" $y - E[y|x]$ (the "errors") to be approximately symmetrically-distributed. This can be assessed with a histogram of the residuals $y - hat(y)$. But the marginal distribution of $y$, such as is assessed with a histogram of $y$, is largely irrelevant in a linear regression or in a generalized linear model (GLM). Similarly, the marginal distributions of the covariates $x_j$ in a regression analysis are usually not relevant.

GLMs are often a useful alternative to transforming variables. A GLM models the mean function $E[y|x]$ using a link function g, so that $g(E[y|x]) = beta' x$. This is not a transformation in the sense discussed here, as the data themselves are not transformed. A GLM has distinct mean and variance functions, providing great flexibility in specifying the model. The conditional mean and conditional variance can be specified with care to best fit a particular population.  We will have a separate document providing a thorough overview of GLMs, so do not discuss them further here.

Another reason for transforming variables is to make the results more interpretable. Most commonly, log transformations are used for this purpose. Log transforms convert multiplication to addition. Many physical, biological, and social processes are better described by multiplicative relationships than additive relationships. Thus, log transforming the independent and/or dependent variables in a regression analysis may produce a fitted model that is both more interpretable, and that may provide a better fit to the data.

A special case of using log transformations is a "log/log" regression, in which both the dependent and one (or more) of the independent variables are log transformed. In this case, the coefficients can be interpreted in terms of the percent change in the mean of the dependent variable corresponding to a given percent change in an independent variable.

To see this, suppose that we have a non-stochastic relationship $log(y) = beta_0 + beta_1 log(x)$, and let $x_1 = x(1 + q)$, so that $x_1$ is $100 times q$ percent different from $x$. Then

$
log(y_1) = beta_0 + beta_1 log(x_1) = beta_0 + beta_1 log(x) + beta_1 log(1+q),
$

and so

$
log(y_1) - log(y) = beta_1 log(1+q).
$

By linearization, $log(y_1) - log(y) = log(1 + (y_1-y)/y) approx (y_1 - y)/y$, and $log(1+q) approx q$. Therefore we have $(y_1 - y)/y approx q beta_1$.

= Categorical variables

Categorical variables can be _nominal_ or _ordinal_, with nominal variables having no ordering or metric information whatsoever, while ordinal variables have an ordering, but there is no precise quantitative meaning to the levels of the variable beyond the ordering. For example, country of birth (US, China, Canada) is a nominal variable, whereas if someone is asked to state their views regarding a policy as being "negative", "neutral", or "positive", then this is ordinal.

In a regression analysis, quantitative, semi-quantitative, and ordinal variables can be modeled directly. But a nominal variable cannot be included directly in a regression, as it must first be "coded". The usual way of doing this is to select one level of the variable as the _reference level_, and then create "dummy" or "indicator" variables for each of the other levels. For example, if a nominal variable $x$ can take on values "A", "B", "C", and we choose level "A" to be the reference level, then we create two indicators, $z_1 = cal(I)(x=B)$ and $z_2 = cal(I)(x=C)$. We cannot also include $cal(I)(x=A)$, in the same regression, since these three indicators sum to 1 and therefore are colinear with the intercept (we could include all three indicators and omit the intercept, but then if there were another categorical variable in the model, we would need to omit one of its categories as a reference level).

There are other ways to code nominal variables in a regression, but the "reference category" approach described above is by far the most common, and is the default in most software. In fact, all standard coding schemes are equivalent via linear change of variables, so we are fitting the same model regardless of which coding scheme is chosen.

Regression coefficients for dummy variables must be interpreted in light of the coding scheme. If the standard reference category scheme is used, then the coefficients are interpreted as contrasts between one non-reference category and the reference category. For example, in the example given above, the regression coefficient for $z_1$ captures the difference in mean values for a case with $x=B$ relative to a case with $x=A$, when all other covariates in the model are equal.

There are some exceptional cases where indicators for all levels of a categorical variable can be included in a model, despite being perfectly collinear.  One such situation would occur when fitting a linear model using methods that do not require a non-singular design matrix (e.g. by employing the pseudo-inverse).  In this case, we get a fitted coefficient vector $hat(beta)$ and its estimated inverse variance/covariance matrix $hat(Psi)^(-1)$ is singular.  However many contrasts of interest may be well-defined in spite of the vector $beta$ not being identified.  Another setting where it is not essential to exclude a reference category is when using a penalized fitting method such as ridge regression or the lasso.  For example, the lasso would automatically drop any redundant covariates early in the solution path.

= Moderation and interactions

An "additive regression" is one in which the expected value of the response variable (possibly after a transformation) is expressed additively in terms of the covariates. A linear mean structure is additive, since

$
E[y | x_1, ..., x_p] = beta_0 + beta_1 x_1 + dots.h.c + beta_p x_p.
$

A more general additive model is:

$
E[y | x_1, ..., x_p] = g_(1)(x_1) + dots.h.c + g_(p)(x_p),
$

where the $g_j$ are functions $RR -> RR$. Models of the second form given above can be estimated using a framework called "GAM" (Generalized Additive Models).  If the mean function is differentiable, then additivity is equivalent to all of the mixed partial derivatives $frac(partial^2 E[y | x_1, ..., x_p], partial x_j partial x_k, style: "horizontal")$ being identically zero when $j != k$.

In an additive model, the change in the mean $E[y]$ associated with changing one covariate by a fixed amount does not depend on the values of the other covariates. For example, in the GAM, if we observe $x_1$ to change from $a$ to $b$, then the expected value of $y$ changes by $g_1(b) - g_1(a)$. This change is universal in the sense that its value does not depend on the values of the other covariates $x_2, ... x_p$.

_Effect modification_ arises when the difference of means resulting from a change in one covariate is not invariant to the values of the other covariates. Put another way, we can say that the "effect" of one covariate is "modified" or "moderated" by the value of another covariate. In practice, effect modification is usually modeled by including _interactions_ in the model, where an interaction can be defined as a product of two or more covariates or derived terms.  For example, we may have the
mean structure

$
E[y | x_1, ..., x_p] = beta_1 x_1 + beta_2 x_2 + beta_3 x_1 x_2.
$

In this model, the parameters $beta_1$ and $beta_2$ are the _main effects_ of $x_1$ and $x_2$, respectively. If we observe $x_1$ to change from 0 to 1, then $E[y | x_1, ..., x_p]$ changes by $beta_1 + beta_3x_2$ units. Note that in this case, the change in $E[y | x_1, ..., x_p]$ corresponding to a specific change in $x_1$ depends on the value of $x_2$, so is not universal in the sense described above.

Including products of covariates in a statistical model is the most common way to model an interaction. But the notion of an interaction, as defined above, is much more general than what can be expressed just by
including products of covariates in the linear predictor.

Focusing on interactions of the product type, a regression model with interactions can be represented by including products of two, three, or more variables, or by including products of transformed variables. For example $log(x_1) dot.c sqrt(x_2-2)$ is an interaction between $x_1$ and $x_2$. If basis functions or categorical variables are present, things can get
complicated:

- If $x_1$ is categorical and $x_2$ is quantitative, then $x_1$ will be represented in the model through dummy variables $z_1, ..., z_q$. The interaction of $x_1$ and $x_2$ is the set of products $x_2z_1, x_2z_2, ..., x_2z_q$.

- If $x_1$ and $x_2$ are both categorical, and we represent $x_1$ with dummy variables $w_1, ..., w_q$, and we represent $x_2$ with dummy variables $z_1, ..., z_(q')$, then the interaction is the set of all $q dot.c q'$ products $w_1 dot.c z_1, w_1 dot.c z_2, ..., w_2 dot.c z_1, w_2 dot.c z_2, ..., w_q dot.c z_(q').$

- If $x_1$ is represented using three basis functions $f_1$, $f_2$, and $f_3$, then the interaction of $x_1$ with another quantitative variable $x_2$ is represented by the terms $f_1(x_1) dot.c x_2, f_2(x_1) dot.c x_2, f_3(x_1) dot.c x_2.$

One challenge that arises when working with interactions is that people struggle to interpret the regression parameters (slopes) of the fitted models. This problem can be reduced by centering all the covariates (or at least by centering the covariates that are present in interactions).

If the covariates are centered, and we work with the mean structure $E[y | x_1, ..., x_p] = beta_1 x_1 + beta_2 x_2 + beta_3 x_1 x_2$, then $beta_1$ is the rate at which $E[y| x_1, ..., x_p]$ changes as $x_1$ changes, as long as $x_2 approx 0$. Similarly, $beta_2$ is the rate at which $E[y| x_1, ..., x_p]$ changes as $x_2$ changes, as long as $x_1 approx 0$. Roughly speaking, when $x_1$ and $x_2$ are close to their means (which are both zero due to centering), then $beta_1$ and $beta_2$ can be interpreted like main effects in a model without interactions. As we move away from the mean, we need to consider the interaction, so the change in $E[y | x_1, ..., x_p]$ corresponding to a unit change in $x_1$ is $beta_1 + beta_3 x_2$, and the change in $E[y | x_1, ..., x_p]$ corresponding to a unit change in $x_2$ is $beta_2 + beta_3 x_1$.

There is a connection between interactions and derivatives. The "regression effect" of $x_j$ can be defined in very general terms as the derivative $frac(d E[y], "dx"_j, style: "horizontal")$. In an additive model, $frac(d E[y], "dx"_j, style: "horizontal")$ is a constant, i.e. it does not depend on the value of $x_k$ for $k != j$. If an interaction between $x_j$ and $x_k$ is present, then $frac(d E[y], "dx"_j, style: "horizontal")$ will depend on $x_k$.

There are two main reasons why it is often a good idea to center covariates that are to be included in interactions:

- If the covariates are centered, then the main effects in a model with interactions have clear interpretations as the rate of change of $E[y]$ corresponding to a unit change in one explanatory variable, when the other explanatory variables are close to their means.

- When the covariates are not centered, variables formed as products, e.g. $x_1x_2$, have complex colinearity properties with other variables, especially with $x_1$ and $x_2$. This can lead to very large standard errors for the main effects, or to settings where models converge slowly or not at all. Often these convergence problems can be easily resolved by centering variables.

It is important to note that main effects have no meaningful interpretation if interactions are present and the covariates are not centered. For example, suppose that $y$ is blood pressure, $x_1$ is body mass index (BMI), and $x_2$ equals 1 for females and 0 for males. We then fit the working model

$
E[y] = beta_0 + beta_1 x_1 + beta_2 x_2 + beta_3 x_1 x_2.
$

In this case, if we do not center the covariates, then the main effect $beta_2$ would mathematically represent the expected difference in blood pressure between a female with BMI equal to zero and a male with BMI equal to zero. Since it is not possible to have BMI equal to zero, this interpretation is meaningless. On the other hand, if we were to center the covariates (including the binary covariate indicating female sex), then $beta_1$ would be equal to the weighted average of the rate of change in $E[y | x_1, ..., x_p]$ per unit change in $x_2$ for females and for males, when $x_1$ is near its mean value, weighted by the proportions of females and males.

We can work through the above example in more detail. If $overline(f)$ is the proportion of females, then after centering the sex variable, the coding for $x_2$ becomes $x_2 = 1-overline(f)$ for females, and $x_2 = -overline(f)$ for males. The regression equation for females can be rearranged to

$
E[y | x_1, ..., x_p] = beta_0 + (1 - overline(f)) beta_2 + (beta_1 + beta_3 (1-overline(f))) x_1
$

and similarly for males we get

$
E[y | x_1, ..., x_p] = beta_0 - overline(f) beta_2 + (beta_1 - beta_3 overline(f))x_1
$

The weighted averages of the BMI slopes for females and males are

$
overline(f) (beta_1 + beta_3 (1 - overline(f))) + (1 - overline(f))(beta_1 - beta_3 overline(f)) = beta_3.
$

Another important thing to note is that the interpretation of the interaction coefficient itself is completely unrelated to how the variables are centered. As shown below, regardless of how we center $x_1$ and $x_2$, $beta_3$ is always the coefficient of $x_1 x_2$.

$
E[y | x_1, ..., x_p] =& beta_1(x_1-c_1) + beta_2(x_2-c_2) + beta_3(x_1-c_1)(x_2-c_2)\
=& beta_3c_1c_2 -beta_1c_1 - beta_2c_2 + (beta_1 - c_2 beta_3)x_1 + (beta_2-c_1beta_3)x_2 + beta_3x_1x_2.
$

Another debate that comes up when working with interactions is whether it is necessary to include all nested "lower order terms" when including an interaction term in a model. For example, if $x_1x_2$ is included in a model, must we also include $x_1$ and $x_2$ as main effects? There are different points of view on this. One argument is that $x_1$, $x_2$, and $x_1x_2$ are just three covariates, and can be selected or excluded from a model independently. However, many variable selection procedures enforce a _hereditary constraint_ in which main effects cannot be dropped in a model selection process if their interaction is included.

= Mediation analysis

A conventional regression model focuses exclusively on how the covariates predict the outcome, not on how the covariates predict each other.  _Mediation analysis_ posits that the relationship between _exposures_ $X$ and an outcome $Y$ can flow through _mediators_ $M$, giving rise to a causal diagram (a directed acyclic graph) $X -> M -> Y$.  In this setting, $X$, $M$, and $Y$ are all observed, and ideally $X$ is assigned through randomization.  It is possible to conduct a mediation analysis with fully observational data, but it is important to remember that unobserved confounders can obscure the true mediation relationship, or create a false one even when no mediation is present.

Here we focus on model-based mediation analysis employing regression, as illustrated #link("https://imai.fas.harvard.edu/research/files/BaronKenny.pdf")[here].  This approach can be understood through the use of _potential outcomes_.  Suppose that for each subject $i$, there are potential outcomes for the mediator $M$ as a function of the exposure $X$.  Denote these mediator potential outcomes as $M_(i)(X=x)$.  We get to observe $M_(i)(X=X_i)$, all the other points on the $M_(i)(dot.c)$ function are "counterfactual".  Similarly, there are potential outcomes for $Y$ as a function of the mediator and the exposure, $Y_(i)(M=m, X=x)$.  We observe $Y_(i)(M=M_i, X=X_i)$ with all other values of $Y_(i)(dot.c, dot.c)$ being counterfactual.

The main goals of mediation analysis are to identify the _direct_ and _indirect_ effects of $X$ on $Y$.  The indirect effects of $X$ on $Y$ are posited to be mediated by $M$, while the direct effects are not.  Let $x_0$, $x_1$ denote two values of the exposure that we wish to compare.  The direct effect is defined as

$
E[Y_(i)(M=M(x_0), X=x_1) - Y_(i)(M=M(x_0), X=x_0)],
$

the indirect effect is defined as

$
E[Y_(i)(M=M(x_1), X=x_0) - Y_(i)(M=M(x_0), X=x_0)],
$

and the total effect is defined as

$
E[Y_(i)(M=M(x_1), X=x_1) - Y_(i)(M=M(x_0), X=x_0)].
$

Intuitively, the direct effect blocks the effect that changing $X$ from $x_0$ to $x_1$ would have on the mediator, and the indirect effect blocks all effects of changing $X$ from $x_0$ to $x_1$ except those that go through the mediator.

Estimation proceeds by building models for $E[Y|X, M]$ and $P(M|X)$ using the observed data.  We can then impute any counterfactual values needed when forming estimates of the direct and indirect effects.  For example, to impute $Y_(i)(M=M(x_1), X=x_0)$, we replace $M(x_1)$ with a sample from the fitted model $hat(P)(M|X=x_0)$.  To account for the fact that $hat(P)(M|X)$ is an estimate of the true model, we should usually either bootstrap the data for each imputation (which is expensive), or, if $hat(P)$ is parameterized by parameters $theta$ for which the estimates $hat(theta)$ are approximately unbiased with estimated variance/covariance matrix $hat(Psi)$, we can replace $hat(theta)$ with $tilde(theta) = hat(theta) + hat(Psi)^(-frac(1, 2, style: "horizontal"))eta$, where $eta$ is an iid standard normal random vector with the same dimension as $theta$.  This is sometimes described as a "quasi-Bayes" approach.

The models for $P(Y|M, X)$ and $P(M|X)$ can be extended to include additional covariates $Z$, in which case we would have $P(Y|M, X, Z)$ and $P(M|X, Z)$.  These could be possible confounders or precision variables.  Furthermore, $Z$ could interact with $X$ in the model for $M$ or in the model for $Y$, giving rise to _moderated mediation_.  This allows the strength of mediation $X -> M -> Y$ to vary with the value of $Z$.

= Causal roles of covariates

Covariates can be play various causal roles, including being: _exposures_, _treatment variables_, _confounders_, _control variables_, _moderators_, _colliders_, and _mediators_, among other roles. These terms can refer to unobserved variables as well as to variables that are available to include in an analysis. It is rarely possible to identify the causal role of every covariate in a proposed analysis. In most cases, analysis of the data cannot fully resolve these roles.  External knowledge or substantive theory are the main basis for identifying how variables are causally related.

The diagram below shows an _exposure_ $X$ for an _outcome_ $Y$, along with a _confounder_ $Z$. A confounder is a _common cause_ of the exposure and the outcome, and should normally be included in the regression analysis to reduce bias.

#align(center + horizon)[
#diagram(
  spacing: (20mm, 15mm),
  node-outset: 3pt,
  node-corner-radius: 5pt,
  node((0, 0), [$X$], name: <a>, fill: blue.lighten(70%)),
  node((2, 0), [$Y$], name: <b>, fill: blue.lighten(70%)),
  node((1, 1), [$Z$], name: <c>, fill: blue.lighten(70%)),
  edge(<a>, "->", <b>),
  edge(<c>, "->", <a>),
  edge(<c>, "->", <b>),
)]

A _precision variable_ $Z$ explains some of the variation in an outcome $Y$, and is unrelated to the exposure $X$.  Including a precision variable in an analysis generally increases power/precision but has no impact on bias.

#align(center + horizon)[
#diagram(
  spacing: (20mm, 15mm),
  node-outset: 3pt,
  node-corner-radius: 5pt,
  node((0, 0), [$X$], name: <a>, fill: blue.lighten(70%)),
  node((1, 0), [$Y$], name: <b>, fill: blue.lighten(70%)),
  node((1, 1), [$Z$], name: <c>, fill: blue.lighten(70%)),
  edge(<a>, "->", <b>),
  edge(<c>, "->", <b>),
)]

A _mediator_ $Z$ lies on the causal pathway between an exposure $X$ and an outcome $Y$.  Including mediators may mask the role of the exposure $X$, but also may explain its mechanism.

#align(center + horizon)[
#diagram(
  spacing: (20mm, 15mm),
  node-outset: 3pt,
  node-corner-radius: 5pt,
  node((0, 0), [$X$], name: <a>, fill: blue.lighten(70%)),
  node((1, 0), [$Z$], name: <b>, fill: blue.lighten(70%)),
  node((2, 0), [$Y$], name: <c>, fill: blue.lighten(70%)),
  edge(<a>, "->", <b>),
  edge(<b>, "->", <c>),
  edge(<a.north>, "->", <c.north>, bend: +40deg)
)]

A _moderator_ ($Z$), also called an _effect modifier_, is a variable that changes the relationship between the exposure $X$ and the outcome $Y$.  Including moderators in analysis can reveal _effect heterogeneity_.

#align(center + horizon)[
#diagram(
  spacing: (20mm, 15mm),
  node-outset: 3pt,
  node-corner-radius: 5pt,
  node((0, 0), [$X$], name: <a>, fill: blue.lighten(70%)),
  node((2, 0), [$Y$], name: <b>, fill: blue.lighten(70%)),
  node((1, 1), [$Z$], name: <c>, fill: blue.lighten(70%)),
  edge(<a>, "->", <b>),
  edge(<c>, "->", (1,0)),
)]

A _collider_ ($Z$) is a variable that is caused by the exposure $X$ and the outcome $Y$.  Including a collider in the model introduces bias.

#align(center + horizon)[
#diagram(
  spacing: (20mm, 15mm),
  node-outset: 3pt,
  node-corner-radius: 5pt,
  node((0, 0), [$X$], name: <a>, fill: blue.lighten(70%)),
  node((2, 0), [$Y$], name: <b>, fill: blue.lighten(70%)),
  node((1, 1), [$Z$], name: <c>, fill: blue.lighten(70%)),
  edge(<a>, "->", <b>),
  edge(<c>, "<-", <a>),
  edge(<c>, "<-", <b>),
)]

= Sensitivity analysis for unmeasured confounding

When we fit a regression model using observational data, the associations identified by the model may not represent causal effects due to the possibility of unmeasured confounding.  It is sometimes useful to quantify how strong an unmeasured confounder would need to be in order to strongly or completely attenuate an effect of interest.  This analysis is relatively easy to perform in the setting of ordinary least squares, which we develop here.

Suppose our model based on observed data produces estimates of the form $hat(beta) = M_(x x)^(-1)M_(x y)$, where $M_(x x) = frac(X'X, n, style: "horizontal")$, and $M_(x y) = frac(X'y, n, style: "horizontal")$.  Now suppose that there is an unmeasured confounder $z$, and we augment $M_(x x)$ and $M_(x y)$ to accommodate $z$.  Specifically,

$
tilde(M)_(x x) =& mat(M_(x x), v; v', 1)\
tilde(M)_(x y) =& mat(M_(x y); s),
$

where $v = frac(X'z, n, style: "horizontal")$ and $s = frac(y'z, n, style: "horizontal")$, which are unknown.  A necessary and sufficient condition for $tilde(M)_(x x)$ to be PSD is that $v' M_(x x)^(-1)v <= 1$.  Also, the partial $R^2$ of $z$ with respect to $X$ is $v' M_(x x)^(-1)v$.  Without loss of generality, we take $z$ to have unit variance (hence $tilde(M)_(x x)[p+1, p+1] = 1$) and zero mean, so $v[1] = 0$ where the first covariate is the intercept.

Now suppose that there is a subset $iota subset {1, .., p}$, so that $beta[iota]$ are the coefficients corresponding to the effects of interest.  That is, we wish to assess whether $beta[iota]$ would be attenuated or even become zero in absence of confounding from $z$.  For any given $v$, the contribution of the effects of interest can be represented by $L(v, s) = norm(X[:, iota]beta_(v s)[iota])^2$, where $beta_(v s)$ are the coefficients when the moments are augmented with $v$ and $s$.  Further, $L$ can be analytically minimized with respect to $s$ for fixed $v$, since $L(v, s)$ is a quadratic polynomial in $v$.  This allows rapid exploration of the space of possible values of $v$, so as to assess whether $L(v, s)$ can be reduced to a very small value or even to zero, while still having a moderate $R^2$ between $z$ and $x$ (given by $v' M_(x x)^(-1)v$), and a moderate level of correlation between $z$ and $y$ (given by $s = frac(y'z, n, style: "horizontal")$ when $y$ is standardized).

= Automating model specification

Many modern regression methods, including most methods from "machine learning", can be seen as aiming to automate the process of model building, specifically with regard to non-linear and non-additive structure.  Here we briefly review two such methods.

== MARS/EARTH

MARS (multivariate adaptive regression splines), also known as EARTH (extended additive regression through hinges) is a method for adaptively constructing multivariate basis functions. We will only describe the approach at a high level here. A _hinge function_ is a function of a single variable of the form
$h(x) = "max"(x-a, 0)$ or $h(x) = "min"(x-a, 0)$. In EARTH, multivariate regression functions are constructed by multiplying hinges for different variables, and nonlinearity can be obtained by summing and/or taking products of hinge functions of a single variable.

EARTH is a greedy algorithm that sequentially searches through the space of basis functions derived as products of hinges. It can capture additive and non-additive relationships. While EARTH remains useful, some drawbacks of this approach have been noted. One drawback is that as a greedy algorithm, it typically cannot achieve the statistical performance of methods that use ensembles or regularization to more efficiently manage the bias/variance tradeoff. A second weakness of EARTH is that there is no rigorous way to perform statistical inference on models fitted using EARTH-like methods.

== Kernel regression

Since the 1990's new approaches to regression based on _kernels_ have become increasingly widely used. These approaches automatically incorporate non-additivity and non-linearity into the fitted regression function, and use regularization to optimize the fitted values based on the bias/variance tradeoff. They share some properties with regression splines, but do not require specification of an explicit model formula and may perform better in
higher dimensions due to their use of regularization.

As an aside, note that the term "kernel" in statistics can have different meanings, and there is a different approach to nonparametric regression based on using kernel weights to localize a regression procedure. That is a different and unrelated use of the term "kernel" to what we are discussing here.

In the present setting, a kernel is a bivariate function $K(dot.c, dot.c)$ mapping $RR^p times RR^p -> RR$. The arguments to the kernel function are covariate vectors, so $p$ is the number of explanatory variables in the model. The kernel must be _positive semi-definite_ meaning that $K(x, x) >= 0$ for all $x in RR^p$.

Two common choices for kernel functions are the _squared exponential kernel_ with bandwith $omega$

$
K(u, v) = exp(-norm(u - v)^2 / (2 omega^2))
$

and the _polynomial kernel_ of degree $d$

$
K(u, v) = (1 + u' v)^d.
$

Each of these kernel functions has a tuning parameter: $omega > 0$ in the first case and $d in 1, 2, ...$ in the second case.

Two ways to use a kernel to build a regression model are as follows. Let $bold(K)$ denote the $n times n$ kernel matrix defined by

$
bold(K)_(i j) = K(x_i, x_j),
$

where $x_i$ and $x_j$ are the covariate vectors for the $i^"th"$ and $j^"th"$ covariates. Note that this is a very large matrix and is expensive to produce and store (usually in regression analysis we construct $n times p$ and $p times p$ matrices but avoid constructing any $n times n$ matrices).

_Kernel ridge regression_ (KRR) estimates coefficients using the ridge-like estimator

$
hat(alpha) = (bold(K)^2 + lambda I)^(-1) bold(K) y,
$

which minimizes the criterion

$
||y - bold(K) alpha||^2 + lambda alpha' bold(K) alpha.
$

It turns out that $alpha' bold(K) alpha$ is a form of regularization in that it shrinks the fitted values toward the nullspace of the Reproducing Kernel Hilbert Space (RKHS) corresponding to the kernel K. More concretely, $alpha' bold(K) alpha$ measures and penalizes the non-smoothness of the function $x -> sum_i alpha_i K(x,x_i)$.

A variant of kernel ridge regression is _Kernel Principal Components Regression_ (KPCR), which finds a limited number of leading eigenvectors of $bold(K)$ and uses them as covariates in an OLS or ridge regression. This would be computationally expensive if done directly, but there is an efficient class of algorithms (including the Lanczos method) for finding a limited number of leading eigenvectors of a large symmetric matrix that are much faster than calculating all of the eigenvectors.

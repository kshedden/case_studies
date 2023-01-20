# Additional regression topics

This document covers several topics that are relevant for many kinds
of regression analysis.

## Basis functions

Many approaches to regression analysis relate the expected value of
the response variable $y$ to a linear predictor $x^\prime \beta$,
formed from the covariates $x \in {\cal R}^p$ using coefficients
$\beta \in {\cal R}^p$.  Examples include the linear mean structure
model $E[y|x] = x^\prime\beta$ and the single index model with link
function $g$, $E[y|x] = g^{-1}(x^\prime\beta)$.

The linear mean structure model is linear in two senses -- the
conditional mean of $y$ given $x$ is linear in $x$ for fixed $\beta$,
and it is linear in $\beta$ for fixed $x$.  These two forms of
linearity have very different implications.

Linearity in $\beta$ for fixed $x$ makes it much easier to
characterize the theoretical properties of the estimation process.
Linearity in $\beta$ generally also makes it easier to develop
algorithms to compute the estimates.

Linearity in $x$ for fixed $\beta$ is sometimes cited as a weakness of
this type of model.  People incorrectly argue that models with this
property are only suitable for systems that behave linearly, and since
most natural and social processes are not linear, such models are
sometimes claimed to have limited utility.

However the (apparent) linearity of the mean structure model in the
covariate vector $x$ is easily overcome.  Given a covariate $x$, say a
person's age, it is possible to include both $x$ and $x^2$ as
covariates in a "linear" model.  This retains the benefits of linear
estimation, while allowing the model for the conditional mean function
to be non-linear in the covariates.

Including powers of covariates (like $x^2$) as regressors is a method
known as *polynomial regression* may be the earliest example
of a general technique utilizing *basis functions* to incorporate nonlinearity
into regression analyses.  A family of basis
functions is a collection of functions $g_1, g_2, \ldots$, each from
${\cal R}\rightarrow {\cal R}$, such that we can include
$g_1(x), g_2(x), \ldots$ as covariates in a model in place of $x$.  This allows the
fitted mean function to take on any form that can be represented as a
linear combination $\beta_1 g_1(x) + \beta_2 g_2(x) + \cdots$.  Using
a large collection of basis functions allows a wide range of
non-linear forms to be represented.  Basis functions thus
allow non-linear mean structures to be fit to data using linear
estimation techniques.

While the basis function approach is very powerful, some families of
basis functions have less desirable properties.
For example, polynomial basis functions ($x^2$, etc.) can be badly scaled
(e.g. this happens when raising a large number to a high power).
Also, polynomial basis functions can be highly colinear with each
other, depending on the range of the data.  Finally, polynomial basis
functions are not "local", meaning that when using polynomial basis
functions, the fitted value at a point $x$ can depend on data values
$(x_i, y_i)$ where $x_i$ is far from $x$.

_Local basis functions_ are families of functions such that each
element in the family has limited support (that is, the functions are
exactly zero outside of a fairly small interval).  One of the most popular forms
of local basis functions is _splines_, specifically, _polynomial
splines_.  Roughly speaking, a polynomial spline is a continuous and
somewhat smooth function that has bounded support, and is a piecewise
polynomial.  We will not derive the mathematical form of a polynomial
spline here in detail.  Polynomial splines are "somewhat smooth" in
that they have a finite number of continuous derivatives (often the
function, and its first and second derivative are continuous, but the
third and higher-order derivative are not continuous).

Other popular families of basis functions are
wavelets, Fourier series, radial basis functions, and
higher-dimensional basis functions formed via tensor products of
univariate basis functions such as splines.

When working with basis functions, it is important to remember that
terms derived from the same parent variable cannot vary independently
of each other.  For example, age determines ${\rm age}^2$ and
vice-versa.  This means that in general, the coefficients of variables
in a regression model that uses basis functions may be difficult to
interpret (i.e. the "effect" of age is represented through the
coefficiencts for ${\rm age}$ and for ${\rm age}^2$).  There are many
ways to resolve this using plots.  For example, if we have a model
relating BMI to age and sex, and we use basis functions to capture a
non-linear role for age, we can make a plot showing the fitted values
of $E[{\rm BMI} \, | \, {\rm age}, {\rm sex}]$ plotted against age,
for each sex.

Using splines or other families of basis functions is a very powerful
technique because it allows familiar methods to be used in a much
broader range of settings, simply by augmenting the regression design
matrix with additional columns.  However in recent years even more
powerful approaches to accommodating non-linearities in regression
modeling have been devised that combine the use of basis functions
with smoothing penalties.  These approaches can do a better job of
capturing non-linearity while also controlling the degree of
smoothness.  However they cannot use standard algorithms for fitting,
or standard methods for inference.  A class of models known as
_generalized additive models_ provides a practical framework for
regression analysis using penalized splines, and addresses both the
computational and inferential aspects of performing this type of
analysis.


## Transformations

A *transformation* in statistics refers to any setting in which a
function $f$ is applied to a variable being analyzed.  It is easy to
apply transformations, and there are many principled reasons for
transforming data.  But it is difficult to come up with a generally
applicable justification to use when applying transformations in a
broad range of settings.

In regression analysis, transformations can be applied to the
dependent variable $y$, to one or more of the independent variables
$x_j$, or to both the independent and dependent variables
simultaneously.  One reason for transforming the data is to induce it
to fit into a pre-existing regression framework.  For example,
ordinary least squares (OLS) is most efficient when the conditional
mean $E[y|x]$ is linear in $x$, and the conditional variance ${\rm
Var}[y|x]$ is constant.  Sometimes, applying a transformation such as
replacing $y$ with $\log(y)$ will induce linearity of the mean
structure and homoscedasticity of the variance structure.  However
achieving linearity of $E[y|x]$ and achieving homoscedasticity
(constant conditional variance) cannot always be achieved with the
same transformation.

Methods for automating the process of selecting a transformation in
linear regression have been proposed, the most well-known being the
Box-Cox method.  However this only automates the process of selecting
a transformation for the dependent variable $y$.  In general,
transforming variables is a manual process of trial and error.

It is sometimes mistakenly believed that the dependent variable $y$ in
a linear model should marginally follow a symmetric or (even stronger)
a Gaussian distribution.  In general however, the marginal
distribution of $y$ is irrelevant in regression analysis.  In a linear
model, we might like the "errors" $y - E[y|x]$ to be approximately
symmetrically-distributed.  This can be assessed with a histogram of
the residuals, or with a plot of residuals on fitted values.  But the
marginal distribution of $y$, e.g. as assessed with a histogram of
$y$, is largely irrelevant in a linear regression or a generalized
linear model (GLM).  Similarly, the marginal distributions of the
covariates $x_j$ in a regression analysis are usually not relevant.

GLM's are often a useful alternative to transforming variables.  A GLM
models the mean function $E[y|x]$ using a link function g, so that
$g(E[y|x]) = \beta^\prime x$.  This is not a transformation, as the
data themselves are not transformed.  A GLM has distinct mean and
variance functions, providing much more flexibility in specifying the
model.  The conditional mean and conditional variance can be specified
independently to best fit a particular population.

Another reason for transforming variables is to make the results more
interpretable.  Most commonly, log transformations are used for this
purpose.  Log transforms convert multiplication to addition.  Many
physical, biological, and social processes are better described by
multiplicative relationships than additive relationships.  Thus, log
transforming the independent and/or dependent variables in a
regression analysis may produce a fitted model that is both more
interpretable, and that may provide a better fit to the data.

A special case of using log transformations is a "log/log" regression,
in which both the dependent and one (or more) of the independent
variables are log transformed.  In this case, the coefficients can be
interpreted in terms of the percent change in the mean of the
dependent variable corresponding to a given percent change in an
independent variable.

To see this, suppose that we have a non-stochastic relationship
$\log(y) = \beta_0 + \beta_1 \log(x)$, and let $x_1 = x(1 +q)$, so
that $x_1$ is $100\times q$ percent different from $x$.  Then

$$
\log(y_1) = \beta_0 +\beta_1 \log(x_1) = \beta_0 +\beta_1\log(x) +\beta_1\log(1+q),
$$

and so

$$
\log(y_1) - \log(y) = \beta_1\log(1+q).
$$

By linearization,
$\log(y_1) - \log(y) = \log(1 + (y_1-y)/y) \approx (y_1 - y)/y$,
and $\log(1+q) \approx q$.  Therefore we have
$(y_1 - y)/y \approx q\beta_1$.  If the relationship is non-stochastic, as in
a regression model, then the we would start with

$$
E[\log(y) | x] = \beta_0 +\beta_1\log(x),
$$

and using the first order approximation

$$
E[\log(y) |x] \approx \log E[y|x],
$$

the analogous result

$$
(E[\log(y) |x_1] - E[\log(y)|x])/E[\log(y)|x] \approx q\beta
$$

is obtained.

## Categorical variables

Categorical variables can be nominal or ordinal, with nominal
variables having no ordering or metric information whatsoever, while
ordinal variables have an ordering, but there is no precise
quantitative meaning to the levels of the variable beyond the
ordering.  For example, country of birth (US, China, Canada) is a
nominal variable, whereas if someone is asked to state their views
regarding a policy as "negative", "neutral", or "positive", this can
be seen as being ordinal.

In a regression analysis, quantitative, semi-quantitative, and ordinal
variables can be modeled directly, e.g.  by including an additive term
$\beta x$ in a regression.  Alternatively, we can use basis functions
or transformed versions of $x$ to construct a more flexible model.  A
purely nominal variable cannot be included directly in a regression,
as it must be "coded".  The usual way of doing this is to select one
level of the variable as the _reference level_, and then create
"dummy" or "indicator" variables for each of the other levels.  For
example, if a nominal variable $x$ can take on values "A", "B", "C",
and we choose level "A" to be the reference level, then we create two
indicators, $z_1 = I(x=B)$ and $z_2 = I(x=C)$.  We cannot also include
$I(x=A)$, in the same regression, since these three indicators sum to
1 and therefore are colinear with the intercept (we could include all
three indicators explicitly and then not include the intercept, but
then if there were another categorical variable in the model, we would
need to omit one of its categories as a reference level).

There are other ways to code nominal variables in a regression, but
the "reference category" approach described above is by far the most
common, and is the default in most packages.  In fact, all standard
coding schemes are equivalent via linear change of variables, so we
are essentially fitting the same model regardless of which coding
scheme is chosen.

It can be a challenge in practice to interpret the coefficients
corresponding to dummy variables.  They must be interpreted in light
of the coding scheme.  If the standard reference scheme is used, then
the coefficients are interpreted as contrasts between one
non-reference category and the reference category.  For example, in
the example given above, the regression coefficient for $z_1$ captures
the difference in mean values for a case with $x=B$ relative to a case
with $x=A$, when all other covariates in the model are equal.

## Interactions

An "additive regression" is one in which the expected value of the
response variable (possibly after a transformation) is expressed
additively in terms of the covariates.  A linear mean structure is
additive, since

$$
E[y] = \beta_0 + \beta_1x_1 + \cdots + \beta_px_p.
$$

A more general additive model is:

$$
E[y] = g_1(x_1) + \cdots + g_p(x_p),
$$

where the $g_j$ are functions.  Models of the form given above can be
estimated using a framework called "GAM" (Generalized Additive
Models).

In any additive model, the change in the mean $E[y]$ associated with
changing one covariate by a fixed amount is not dependent on the
values of the other covariates.  For example, in the GAM, if we
observe $x_1$ to change from $a$ to $b$, then the expected value of
$y$ changes by $g_1(b) - g_1(a)$.  This change is universal in the
sense that its value does not depend on the values of the other
covariates ($x_2, \ldots x_p$).

An *interaction* arises when the difference of means resulting from a
change in one covariate is not invariant to the values of the other
covariates.  There are many ways that an interaction can arise, but in
practice we often model an interaction by taking a product of two
variables.  For example, we may have the mean structure

$$
E[y] = \beta_1x_1 + \beta_2x_2 + \beta_3x_1x_2.
$$

In this model, the parameters $\beta_1$ and $\beta_2$ are the *main
effects* of $x_1$ and $x_2$, respectively.  If we observe $x_1$ to
change from 0 to 1, then $E[y]$ changes by $\beta_1 + \beta_3x_2$
units.  Note that in this case, the change in $E[y]$ corresponding to
a specific change in $x_1$ depends on the value of $x_2$, so is not
universal in the sense described above.

Including products of covariates in a statistical model is the most
common way to model an interaction.  But note that the notion of an
interaction, as defined above, is much more general than what can be
expressed just by including products of covariates in the linear
predictor.

Focusing on interactions of the product type, a regression model with
interactions can be represented by including products of two, three,
or more variables, or by including products of transformed variables.
For example $\log(x_1)\cdot \sqrt{x-2}$ is an interaction between
$x_1$ and $x_2$.  If basis functions or categorical variables are
present, things can get complicated:

* If $x_1$ is categorical and $x_2$ is continuous, then $x_1$ will be
represented in the model through dummy variables $z_1, \ldots, z_q$.
The interaction of $x_1$ and $x_2$ is the set of products
$x_2z_1, x_2z_2, \ldots, x_2z_q$.

* If $x_1$ and $x_2$ are both categorical, and we represent $x_1$ with
dummy variables $w_1, \ldots, w_q$, and we represent $x_2$ with dummy
variables $z_1, \ldots, z_{q^\prime}$, then the interaction is the set
of all $q \cdot q^\prime$ products
$w_1\cdot z_1, w_1\cdot z_2, \ldots, w_2\cdot z_1, w_2\cdot z_2, \ldots, w_q\cdot z_{q^\prime}.$

* If $x_1$ is represented using three basis functions $f_1$, $f_2$,
and $f_3$, then the interaction of $x_1$ with another continuous
variable $x_2$ is represented by the terms
$f_1(x_1)\cdot x_2, f_2(x_1)\cdot x_2, f_3(x_1)\cdot x_2.$

One challenge that arises when working with interactions is that
people struggle to interpret the regression parameters (slopes) of the
fitted models. This problem can be minimized by centering all the
covariates (or at least by centering the covariates that are present
in interactions).

If the covariates are centered, and we work with the mean structure
$E[y] = \beta_1x_1 + \beta_2x_2 + \beta_3x_1x_2$, then $\beta_1$ is
the rate at which $E[y]$ changes as $x_1$ changes, as long as
$x_2 \approx 0$.  Similarly, $\beta_2$ is the rate at which $E[y]$ changes
as $x_2$ changes, as long as $x_1 \approx 0$.  Roughly speaking, when
$x_1$ and $x_2$ are close to their means (which are both zero due to
centering), then $\beta_1$ and $\beta_2$ can be interpreted like main
effects in a model without interactions.  As we move away from the
mean, we need to consider the interaction, so the change in $E[y]$
corresponding to a unit change in $x_1$ is $\beta_1 + \beta_3x_2$, and
the change in $E[y]$ corresponding to a unit change in $x_2$ is
$\beta_2 + \beta_3x_1$.

There is a connection between interactions and derivatives.  The
"regression effect" of $x_j$ can be defined in very general terms as
the derivative $d E[y]/dx_j$.  In an additive model, $d E[y]/dx_j$ is
a constant, i.e. it does not depend on the value of $x_k$ for
$k \ne j$.  If an interaction between $x_j$ and $x_k$ is present, then
$d E[y]/dx_j$ will depend on $x_k$.

There are two main reasons why it is often a good idea to center
covariates that are to be included in interactions:

* If the covariates are centered, then the main effects in a model
with interactions have clear interpretations about the rate of change
of $E[y]$ corresponding to a unit change in one variable, when the
other variables are close to their means.

* When covariates are not centered, variables formed as products,
e.g. $x_1x_2$, have complex colinearity properties with other
variables, especially with $x_1$ and $x_2$. This can lead to very
large standard errors for the main effects, or to settings where
models converge slowly or not at all.  Often these convergence
problems can be easily resolved by centering variables.

If the covariates are not centered, it is common that people
misinterpret the estimates of main effects.  In general, main effects
have no meaningful interpretation if interactions are present and the
covariates are not centered.  For example, suppose that $y$ is blood
pressure, $x_1$ is body mass index (BMI), and $x_2$ equals 1 for
females and 0 for males.  We then fit the working model

$$
E[y] = \beta_0 + \beta_1x_1 + \beta_2x_2 + \beta_3x_1x_2.
$$

In this case, if we do not center the covariates, then the main effect
$\beta_2$ would mathematically represent the expected difference in
blood pressure between a female with BMI equal to zero and a male with
BMI equal to zero. Since it is not possible to have BMI equal to zero,
this interpretation is meaningless.  On the other hand, if we were to
center the covariates (including the binary covariate indicating
female sex), then $\beta_1$ would be equal to the weighted average of
the rate of change in $E[y]$ per unit change in $x_2$ for females and
for males, when $x_1$ is near its mean value, weighted by the
proportions of females and males.

We can work through the above example in more detail.  If $\bar{f}$ is
the proportion of females, then after centering the sex variable, the
coding for $x_2$ becomes $x_2=1-\bar{f}$ for females, and
$x_2=-\bar{f}$ for males.  The regression equation for females can be
rearranged to

$$
E[y] = \beta_0 + (1-\bar{f})\beta_2 + (\beta_1 + \beta_3(1-\bar{f}))x_1
$$

and similarly for males we get

$$
E[y] = \beta_0 - \bar{f}\beta_2 + (\beta_1 - \beta_3\bar{f})x_1
$$

The weighted averages of the BMI slopes for females and males are

$$
\bar{f} (\beta_1 + \beta_3(1-\bar{f})) + (1 - \bar{f})(\beta_1 - \beta_3\bar{f}) = \beta_3.
$$

Another important thing to note is that the interpretation of the
interaction coefficient itself is completely unrelated to how the
variables are centered.  As shown below, regardless of how we center
$x_1$ and $x_2$, $\beta_3$ is always the coefficient of $x_1x_2$.

$$
E[y] = \beta_1(x_1-c_1) + \beta_2(x_2-c_2) +
\beta_3(x_1-c_1)(x_2-c_2)\\ = \beta_3c_1c_2 -\beta_1c_1 - \beta_2c_2
+ (\beta_1 - c_2\beta_3)x_1 + (\beta_2-c_1\beta_3)x_2 + \beta_3x_1x_2.
$$

Another debate that comes up when working with interactions is whether
it is necessary to include all nested "lower order terms" when
including an interaction term in a model.  For example, if $x_1x_2$ is
included in a model, must we also include $x_1$ and $x_2$ as main
effects?  There are different points of view on this.  One argument is
that $x_1$, $x_2$, and $x_1x_2$ are just three covariates, and can be
selected or excluded from a model independently.  However, many
variable selection procedures enforce a _hereditary constraint_ in
which main effects cannot be dropped in a model selection process if
their interaction is included.

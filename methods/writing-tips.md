# Writing tips

This page discusses some observations that we have made when reading
assignments for courses like this in the past. Many of these points apply to
any writing about data, but a few of them are specific to this course. It is
very important to keep in mind that most of the points made here are general
guidelines, not absolute rules. In a specific context, you should use your
judgment and violate any of these principles if you feel that it is
appropriate to do so.

## Content/overall approach

- Writing for this course should narrowly focus on data analysis, which
  usually includes (i) motivating your analysis, (ii) describing the relevant
  structure of the data that you are analyzing, (iii) precisely describing the
  methods that you use, including their relevant theoretical properties, and
  (iv) interpreting the results of your analysis and stating your findings.
  Part (i) provides some room for broader discussion, e.g. about the
  scientific or social context, but for this course, this type of discussion
  should be extremely brief. Nearly every sentence of your writing should
  pertain to your data analysis and the rigorous interpretation of your
  findings.

- Scientific writing follows a constrained style, yet there is plenty of
  opportunity to develop your own voice. This is the reason that we do not
  provide templates or example memos for this course. It will be impossible
  for you to write well in an appropriate scientific/academic style if you
  have not read a lot of scientific writing. We will provide some pointers to
  good sources of scientific writing.

- The goal of these memos is to document a specific research question that you
  have devised, and an analysis that you have conducted using the data that we
  provide. You will then present an interpretation of your findings in
  relation to your research question. Little-to-no space should be devoted to
  broader discussion of the scientific domain or statistical issues not
  directly related to your analysis and findings. Since the length of your
  memo is strictly limited, any such discussion will take away from the
  quality and quantity of your specific findings and will lead to the ideas
  presented in your memo being underdeveloped.

- Write in complete sentences, with the content organized logically into
  paragraphs that are usually 3-5 sentences long. Do not submit an outline or
  list of points. You do not need to include section headers but may do so if
  you wish.

- Your memo should be structured around a narrative that is driven by your
  research aim, hypotheses, and findings. Methods, plots, and raw results
  should follow and support the narrative rather than the other way around.
  Concretely, this means that your memo should not read like "First we did X,
  then we did Y, then we did Z...", or "Figure 1 shows, Figure 2 shows, Figure
  3 shows...".

- Most memos will include a brief introductory paragraph containing an
  explicit statement of the research aim along with essential background
  information, one or two paragraphs discussing the methods used in your
  analysis, one or two paragraphs of results, and one or two paragraphs of
  discussion. In the discussion paragraphs you can be a bit more interpretive
  than the results paragraphs, which should be more "matter of fact".

- Do not discuss "future work" or things that you could have done but did not
  do. It may be appropriate to discuss some limitations of your analysis but
  these should be very specific to your analysis and conclusions, not a
  generic list of limitations that would apply in almost any setting.

- Aim for an academic tone. Avoid the excessive use of self-congratulatory
  adjectives like "comprehensive", "nuanced", "advanced", and "intricate" to
  describe your analysis unless you can specifically state why these terms are
  relevant. Let the reader be the judge of whether these terms are applicable.

- You can write in the first person "I found..." or using the "majestic we"
  (presenting your results using the pronoun "we" even though the analysis was
  done by one person).

- Make specific claims and assertions, then back them up with specific
  evidence and arguments. The evidence and arguments should be specifically
  grounded in the data that you have analyzed and the statistical methods that
  you employ.

- Focus on findings that follow directly from your analysis of the data. For
  this class, you do not need to provide independent sources to corroborate or
  support your findings.

- Don't use "hedge words" unnecessarily, but do use qualifications whenever a
  definitive statement cannot be justified. If you find yourself hedging too
  frequently, your findings may be too weak to report.

- Balance confidence and modesty when discussing your findings and their
  implications. Be realistic about what can be accomplished when addressing a
  complex and challenging question with a single limited data set. Avoid
  making sweeping conclusions or suggesting without qualification that your
  findings support specific policy actions to be made in the "real world" or
  reveal major unexpected scientific insights. Remember that one brief
  analysis of one data set will only very rarely allow definitive or deeply
  novel insights to be gained.

- Avoid "straw man" arguments, for example, spending a lot of time discussing
  the difficulty of answering questions that are tangential to the main
  analysis goal.

- Don't enumerate things that are impossible to do with the data that you
  have. Don't waste space writing about things that you did not do, even if
  you think that you have an interesting reason for not having done them.

- Avoid presenting your work as if you followed a set script or recipe. There
  are no pre-defined scripts for data analysis. The flow of your report should
  follow the logic of the question you are aiming to address, and the evidence
  and interpretation that you use to support your claims.

- Avoid presenting your methods as a "laundry list", e.g. "In this memo we use
  least squares regression, residual diagnostics, and principal components
  analysis". Each method should be used in the memo for a specific reason, and
  this reason is more important than the method itself. The writing should
  reflect this, therefore, instead of a laundry list, you might use something
  like "to assess ... we use ..., and then to understand ... we use ...".

- Discussion of how you "cleaned" the data is only relevant if it is needed
  for the reader to understand your question, approach, and conclusions. There
  is no need to document menial data processing tasks that involve no
  statistical judgment. Avoid use of programming language syntax in your
  writing. This includes variable names -- refer to variables using plain
  language.

- Try to make your writing as self-contained as reasonably possible. For this
  course, you can assume that the reader is familiar with the statistical
  methodology that you use -- you do not need to teach the reader about
  statistics from the ground up. Nevertheless, it is important to reinforce
  key aspects of your approach, by explaining the strengths and weaknesses of
  the methods that you employ, and presenting the statistical reasoning behind
  your conclusions.

- Do not simply summarize the lectures or other course materials. These
  materials are not organized around focused research questions, and contain
  much didactic material that is not appropriate for your memos. You can
  borrow many of the ideas from the course materials to support your
  arguments, but your memo should be structured around a specific scientific
  question, which is not usually the way that the lectures and other course
  materials are organized.

- To make the writing more realistic, do not refer to our class as such, or to
  prior assignments.

# Analysis

- Avoid over-summarization of data. In most cases, if you report a statistic
  such as the sample mean based on more than a few hundred observations, you
  have missed an opportunity to delve deeper into the data, e.g. by reporting
  the means for relevant subgroups.

- Provide quantitative context around your key findings. For example, if your
  primary point is about the mean of some measure, it may also be important to
  report its standard deviation and the sample size. Context can also be
  provided by reporting the mean and standard deviation of the variable in
  relevant subgroups, not only in the dataset as a whole.

- Define any measures or variables that you use unambiguously, including their
  scales or measurement units. For example if you are looking at crime, is it
  the absolute number of crimes or the crime rate? When reporting rates,
  always be explicit about the denominator of the rate (e.g. violent crimes
  per 10,000 population per year).

- If you are presenting a regression analysis, unambiguously state what is the
  dependent variable, and what are the independent variables.

- Avoid presenting lists of descriptive statistics without using them to
  support a larger point, e.g. don't report the mean of every variable in a
  dataset without stating what we learn from knowing these means.

- Avoid over-simplified analyses -- most systems in the real-world have
  nonlinear and interactive behavior, and exhibit various forms of
  heterogeneity. When building models, you should generally aim to fit as
  complicated of a model as you have statistical power to estimate. If
  presenting summary statistics for heterogeneous data, present them at an
  appropriately disaggregated or stratified level.

- Don't overstate or misstate the role of Gaussianity in statistical analyses.
  Often it is a minor consideration. It is almost never important that the
  data are Gaussian. To the extent that Gaussianity matters, it is usually
  because it is easier to calibrate tests and confidence intervals when
  aggregateed quantities like descriptive statistics or parameter estimates
  are approximately Gaussian.

## Organization and structure

- Always include a title. The title should convey something specific about the
  main focus of your writing. It should reflect the main substantive question
  that you are addressing, and may also refer to the data that you will be
  using. You generally will not want to mention the analytic methodology (e.g.
  regression analysis) in the title. The title should not have the format
  "Analyzing..." or "Exploring...". You want to give the reader a reason to
  expect that they will learn something specific by reading your memo. Your
  analyses and explorations of the data will not be of interest to a reader
  unless you can tell them up-front what question you will be posing and
  attempting to resolve.

- The introductory paragraph should convey the main topic and focus of your
  writing. Information about the data set, and the methodology that you plan
  to use should also be discussed, but the primary substantive question is the
  most important information and should be covered first. You may need to
  provide a bit of background information to motivate your question, but keep
  the background discussion in the first paragraph brief.

- It should be possible to understand your main ideas by reading the report
  once from top to bottom. Avoid making statements that are only meaningful
  after reading something that comes later in the report (even if it is in the
  very next sentence).

## Language and style

- Avoid the passive voice

- Avoid presenting your work like a recipe or diary, e.g. "First I did ...,
  then I did ...".

- Avoid excessive use of "would", e.g. "the analysis would provide us with
  some insights..." is better written as "the analysis will provide us with
  some insights..." or "the analysis provides us with some insights..."

- Limit ambiguous internal cross-references, e.g. "above", "below", "before"
  (used to refer to something that appears elsewhere in your report). Only do
  this if it is unambiguous what you are referring to, otherwise be more
  specific (e.g. provide a section number if you are using them, or use
  language such as "in the regression analysis presented above we saw
  that...").

- Favor short sentences and paragraphs. If a sentence has more than two
  clauses, consider reworking it into multiple shorter sentences. If a
  paragraph has more than five sentences, consider splitting it into multiple
  shorter paragraphs.

- Don't start a sentence with "And". Avoid starting sentences with "Such", or
  starting non-interrogative sentences with "Which".

- Avoid generic and low-information statements such as stating that something
  is "reasonable" or "interesting".

- Rarely use "etc." unless it is unambiguous what it could be replaced with.

- Use consistent tense. If you are describing something concrete that you did
  (like conducting an analysis of a data set), it may be appropriate to use
  past tense (e.g. "I dropped all rows with missing data"). If you are
  describing a research finding, it is usually preferable to use present tense
  ("smoking and age are positively associated" not "smoking and age were
  positively associated").

## Terminology

- Don't be afraid to use technical terms as needed, but use them precisely and
  only where the exact definition is important. Avoid "jargon" (technical
  language that is used informally and needlessly).

- Since the word "significant" has a specific technical meaning in statistics,
  it is best to avoid using it in its colloquial sense. In general, it is best
  to use the phrase "statistically significant" when you mean that a formal
  statistical uncertainty assessment has been used to support a finding. Using
  "significant" on its own can be ambiguous since it has both a
  common-language meaning and a technical meaning.

- The word "valid" is best avoided. In addition to being vague, this term
  reinforces a binary interpretation of statistical evidence, i.e. that there
  is a sharp line between analytic methods that are right or wrong ("valid" or
  "invalid"). Most statisticians prefer to think in more continuous terms. An
  analysis should be meaningful (otherwise don't present it), but nearly
  always will have limitations. It is neither fully valid nor fully invalid.

- The term "assumption" is very commonly used when describing statistical
  methods, but is arguably over-used. A good alternative term is "condition".
  Your findings will be more meaningful if the conditions approximately hold,
  and less meaningful if they are strongly violated. If you do choose to refer
  to "assumptions", in applied statistics we usually want to argue that the
  assumptions hold "approximately", not "exactly".

- It is important to distinguish "testable assumptions" (which need not really
  be assumptions at all if you carefully assess them), versus "untestable
  assumptions", which you must take on faith. Of course, most assumptions are
  neither completely testable or completely untestable, but fall somewhere in
  between.

- Keep in mind the distinction between "methods" and "models". A "method" is
  any technique used for analyzing and gaining insight from data. A "model" is
  a specific mathematical or computational formalism that aims to describe how
  a system works. For any given class of models, there is a corresponding
  class of methods that can be used to fit those models to data. For example,
  we have "linear models" and "ordinary least squares" -- the former is a
  class of models, the latter is a class of techniques for fitting models to
  data. It is usually not appropriate to write "OLS models", since it confuses
  these two distinct ideas.

- The term "estimation" is used in statistics to refer to the setting where we
  sample data from a population and use the sample to estimate population
  parameters. Because the sample is random, the estimate is random, and we use
  probability to evaluate the performance of the estimator. The term
  "approximation" is used in mathematics to refer to a calculated value that
  is close to some other value, like using 3.14 as an approximation to Pi. The
  error in an approximation is usually deterministic (not random). Some
  statistical methods involve approximation (in addition to estimation), but
  it is usually much more correct to describe the result of a statistical
  procedure as being an estimate not an approximation.

- You will frequently need to discuss interactions in regression models. If
  the response variable is $y$ and covariates $x$ and $z$ interact, you could
  simply write that "there is an interaction between $x$ and $z$". However
  this assumes that the reader has a lot of experience interpreting
  interactions. It is more intuitive to write that "the role of $x$ differs
  based on the value of $z$", or equivalently that "the role of $z$ differs
  based on the value of $x$". If, for example, $z$ is categorical, you could
  write that "the fitted slope of $E[y | x, z]$ on $x$ differs between the
  groups defined by $z$". You would need to adjust this language if, for
  example, there is a link function.

## Causality

- Causality is a fundamental issue in most statistical analysis. In some cases
  such as an controlled experiment with randomization, it is possible to
  interpret associations in causal terms. More often, we have imperfect
  experiments or data from observational studies, so are unable to be strongly
  confident about the extent to which any findings reflect mechanisms in a
  causal sense.

- Many common-language terms have a connotation of causality and should be
  avoided unless you have done something to demonstrate that your findings
  have a causal interpretation. This includes words such as: "impacts",
  "causes", "affects", "determines", and "changes".

- Certain phrases such as "the effect of x on y" can be re-worded to convey
  that the findings are not meant to be interpreted causally. For example, it
  is possible to discuss "the relationship between x and y", or to state that
  "x predicts y".

- When discussing a regression model, it is better to refer to the
  "coefficient for x" or the "slope parameter for x" rather than the "effect
  of x", since the latter implies causality.

- While you should be very careful about casually implying that your findings
  have a causal interpretation, it is usually not necessary to eliminate all
  hints of causal thinking from your writing. Science is primarily concerned
  with identifying causal relationships. Evidence for causality does not need
  to come solely from the data that you are analyzing. In some cases, there
  may be a basis for interpreting results in a causal manner even when the
  data themselves do not demonstrate causality.

- There is a branch of statistics called "causal inference". Methods from this
  field, such as weighting, stratification, and matching can be useful ways to
  assess whether a data set supports causal interpretations of findings.
  However, there is not a bright line between methods from "causal inference"
  and other statistical methods. A regression model with suitable control
  variables may support a causal interpretation, even though it is not viewed
  as being a technique from causal inference.

- In general, demonstrating causality is a trade-off between rigor and power.
  Causal inference tends to favor rigor even when a great sacrifice of power
  results. Mainstream regression analysis favors a different balance in which
  some effort is made to promote a causal interpretation (e.g. by controlling
  for known confounders), but efforts to achieve causality are balanced with
  efforts to maintain power.

## Graphics

- Only include plots or tables that are essential for conveying your message
  and that you explicitly refer to in your report. Do include any graphs that
  are discussed in your report and that are used to support your main
  conclusions. In most cases 1-4 graphs is appropriate for a memo of two
  pages.

- Number all plots and tables and always use these numbers when referring to
  plots/tables in the body of your report. Do not refer to graphs as being
  "above" or "below" a particular location in the text.

- Refer to plots using proper and specific terminology. For example, is the
  plot a scatterplot, a histogram, a time-series plot, a quantile-quantile
  plot, a graph of a function, etc. The horizontal and vertical axes should be
  labeled using terms that connect the graph to the data being analyzed.

- In a statistical analysis, graphs typically are used to reinforce specific
  statistical claims. For example, histograms reflect "probability density",
  and scatterplots often reflect "conditional mean", "heteroscedasticity", or
  "conditional variance". Usually it is appropriate to make these connections
  explicit when writing about a graph.

- Just as with other forms of statistical analysis, graphs are most effective
  when used comparatively. Thus, it is often effective to overlay or juxtapose
  graphs corresponding to comparisons in your analysis. For example, you can
  usually overlay 2-3 histograms or density estimates, or you can capture
  multiple conditional mean relationships by coloring the points in a
  scatterplot corresponding to different subgroups.

- Avoid pie charts in almost all situations

- Bar graphs are a very elementary and limited form of statistical graphic. It
  is almost always possible to create something more informative. For example,
  side-by-side box plots show the mean (which is often what is shown in a bar
  graph), but also convey information about the dispersion. Or, you can
  identify another attribute of the data underlying each bar, and then show a
  scatterplot of the primary feature against this secondary attribute. Another
  good alternative to bar graphs is a dot plot.

## Grammar

- In general, it is better to avoid contractions (e.g. "it's") in this type of
  writing.

- Either place a blank line between paragraphs, or indent each paragraph (or
  do both).

- If you use parenthetical expressions, place a space before each opening
  parenthesis.

## Statistical content

- Statistical analysis usually aims to address a question about a population,
  based on a sample of data from the population. Be sure that it is clear what
  is the relevant population for your analysis, and how the sample was
  obtained.

- Most empirical research involves comparisons, although occasionally the
  primary insight is an "absolute" value. In general, you should aim to make
  comparisons. For example, we might like to compare the influenza death rate
  in children to adults, rather than estimate the absolute death rate in
  either group.

- In almost any report, you should state the sample size for any analyses that
  you are reporting.

- Avoid fitting overly-simple models to large datasets. Underfitting is just
  as serious a problem as overfitting.

- When referring to a p-value always make clear what hypothesis the p-value
  evaluates.

- A correlation is a specific type of statistical measure describing the
  relationship between two variables. In most cases, if you are not talking
  about one of the standard correlation measures (Pearson, Spearman, Kendall)
  it is best to refer to the statistic as an "association".

- Be explicit when discussing percentage changes in variables that are
  themselves percentages. The term "percentage points" can clarify your
  meaning in this setting. For example, if you say that the unemployment rate
  increased by 3% (say, from 3%), it is not clear if you mean that it
  increased by 3 percentage points (to 6%), or by 3% of 3% (to 3.09%).

- Categorical variables must usually be coded into multiple dummy variables
  before being considered in a regression model. The interpretation of the
  coefficients for these dummy variables is only meaningful relative to the
  reference category (if there is one, or more generally it can only be
  interpreted in light of the "coding scheme"). It is critical to consider and
  clearly report the coding scheme when discussing parameter estimates for
  categorical variables.

- It is fine to do a certain amount of "descriptive analysis", but avoid
  providing descriptive statistics for no clear reason.

- Avoid referring to single data values (e.g. outliers, extremes). Statistical
  data analysis is almost always about abstracting away from the specific data
  that you have observed and saying something about the population. Single
  data values in your data set rarely provide much evidence regarding
  properties of the population.

- One of the main contributions that someone with advanced statistical
  training can make in a data analysis is to reveal what the data say after
  removing or controlling for certain forms of unwanted variation, or to
  uncover how various subgroups relate to (or differ from) each other, or to
  adjust for confounding, selection biases, or other measurement artifacts and
  inadequacies. This can be accomplished in many ways, for example, by
  reporting summary statistics on stratified data, or by using a model-based
  adjustment (e.g. regression).

- Be aware of which of your quantitative findings have units and which are
  dimension free. If your results depend on units make sure that the units are
  clearly reported, and avoid drawing conclusions about a value being "big" or
  "small" if the scale of the value is dependent on (arbitrary) measurement
  units.

- Administrative geographical units (e.g. US states and counties) are
  generally very unbalanced in terms of population. When using such an
  unbalanced partition as a grouping variable, it is rarely of interest to
  report on the absolute amount of any variable. For example, in the case of
  U.S. states, virtually everything (GDP, murders, bankruptcies, airports,
  ...) will track with state population. If you rank the states by any of
  these features, you will usually get a list that starts with California,
  Texas, Florida, etc. and ends with Wyoming and Alaska. In nearly all cases,
  when working with unbalanced groups you should report some type of rate or
  other normalized value. Note that this does not imply that you should never
  model absolute quantities. But whether you are reporting modeled results, or
  raw data, you should generally present your findings on a relative scale.

# fitl: A Self-contained fit of linear models with regression diagnostics

I do a lot of performance analysis and other data analysis.  Simple regression
with a handful of columns is often adequate, but has many hazards addressable
through various mitigations & diagnostics I have not found bundled into one,
self-contained tool.  So, I wrote this program, pronounced almost like "fiddle"
originally in C but have now ported to Nim.

It is both a Nim library and command-line utility which is made easy using the
[cligen](https://github.com/c-blake/cligen/) framework.  `nimble install fitl`
should just work, but if it fails then just git clone cligen & spfun and then
`nim c -d:danger --path:cligen --path:spfun fitl` to build the CLI.

# Example use

XXX TODO Do me.

# Truly self-contained - even the linear algebra

BLAS/LAPACK from the Fortran world have many implementations - Netlib-Fortran
compiled, vendor optimized, CUDA/GPU optimized, etc.  This is nice for large
problems, but overkill for simple linear models with <100s of predictor and
<many millions of data points.  Very much *because of* implementation diversity,
there is much potential for "linker/dependency hell".

Meanwhile, `svdx.nim` adequately does a column-major layout SVD in ~42 lines of
non-comment code (46 if you count `basicLA.dot`).  For a long time, for small (p
<~40, moderate n) models, this actually ran faster than the Intel MKL SVD.  By
now Intel may have special cased such use cases.  The point is not that it is
"fast", but that it is not very slow for its use case - hence "adequate".

# Statistical Background

While a couple pages of text cannot realistically replace a course in stats, I
can try to give an overview of "whys" behind certain `fitl` choices.

Ordinary least squares linear models are usually both wrong and practical.  This
is true of linear statistics like "means" in general.  The linearity yields much
computational simplicity & efficiency.  This is why they are practical and this
draws people in to rationalize the approach.  There is kind of a "cottage
industry" to justify the same idea by various avenues.

These justifications are often not fully mathematically adequate which is why
they are "wrong".  E.g., max likelihood and Gaussian error distributions are one
way, while minimum variance unbiasedness in the space of linear estimators is
another.  These are all nice properties, but not enough.  Max likelihood mostly
has nice asymptotic traits, but most samples are not so close to infinite sample
size.  The central limit theorem motivates Gaussian "assumptions", yet real
world distributions are mercilessly contaminated with outliers/mixtures/heavy
tails.  Meanwhile, min var unbiased is nice, but trestricting to only linear
estimators when non-linear exist is merely a computational argument, not more
statistically fundamental.

# Various Mitigations & Diagnostics

Since all the justifications are convenience/wrong, the default assumption for
linear models should be that they are wrong.  This means you need diagnostics to
tell you how wrong.  You also probably want mitigations to limit the impact of
wrongness on conclusions.  `fitl` has a bare minimum of both and I am open to
adding more (and what I consider a bare minimum is more than most provide/do).

## Calculational:

Many real world data sets have linearly co-dependent predictor columns.  This
can arise simply from some being statistically flat and co-linearity with an
almost always included constant intercept column.  More subtle mechanisms also
exist.  Whatever the causes, this breaks simpler-than-SVD approaches to best fit
coefficient determination.  So, SVD is basically mandatory with some kind of
singular value clipping to create a pseudo-inverse based on the design matrix.

While often only flat weights are available, weighted OLS linear regression is a
simple transformation of the flat weight case.  This is supported by the first
non-option argument to `fitl` being a weight or "sigma" column.  `fitl` does not
estimate "errors in variables" models which (but for effective variance style)
require non-linear computational methods.  (I may someday do a non-linear
variant but differentiation, auto- or otherwise is a big enough boost to tempt
one into a more "dependent" program.)

## Inferential:

Estimated parameter covariance is often critical to assessment of fit quality or
inference on parameter/predictor significance.  While an estimate is almost
available for free from the linear OLS solution, this estimate can be very
misrepresentative.  The distribution of y/response errors must really be very
Gaussian for the formula to not mislead.  The computationally expensive but
statistically easy & reliable way to adapt to this problem is estimating Cov(b)
via bootstrap resampling `fitl --boot=25` does at "only" 25x the CPU cost.

## Estimational:

Outliers, long or heavy tails or mixture distributions or non-linear functional
dependencies are so much more common in real world data that they should likely
just be the default assumption.  Since assuming so quickly (immediately?)
destroys computational simplicity & efficiency, people instead punt.

The field of robust statistics has a long history of simple data filtering
workarounds for outliers such as
[trimming](https://en.wikipedia.org/wiki/Trimmed_estimator) and
[Winsoration](https://en.wikipedia.org/wiki/Winsorizing).  `fitl` currently just
does basic trimming in appropriate units with options like `--trim=4.5` &
`--its=2`.  Beware that many iteration can trim away almost all your data!

## Diagnostical:

The most basic goodness-of-fit statistic is R^2 or the achieved Chi-squared.
Mispecified models often show strong serial autocorrelation in fit residuals
from the simple effect of many regions where the points cluster far from the
model.  Failure of trimming to achieve a Gaussian-like distribution also results
in a non-Gaussian raw distribution of residuals.

So, if you are using linear regression the least you should do is check for the
main departures with `fitl --gof` and `--acf`.  For aid in interpretation, the
statistical significance of autocorrelation function is also reported.  (Rank
autocorrelation is not yet done but would not take much convincing!)

Ideally, if you are not in an automated setting, you should visually inspect
residuals with `fitl --resids`.

## Feature Discovery-(ical):

In recent times, "machine learning" has become popular.  This is basically
applied stats with a "Just predict, by hook or by crook" focus rather than more
traditional/inference-focused/"science-oriented" stats.  Consequently, the ethos
is to throw every possible predictor (and sometimes technique in general) at
problems and see what features are discovered.  In the context of multi-linear
regression, this puts a pressure on preventing overfitting often called
[regularization](https://en.wikipedia.org/wiki/Regularization_%28mathematics%29)
or sometimes [model selection](https://en.wikipedia.org/wiki/Model_selection).

Since SVD is already done for reasons of basic estimation sanity, [ridge
regression](https://en.wikipedia.org/wiki/Ridge_regression) aka [Tikhonov
regularization](https://en.wikipedia.org/wiki/Tikhonov_regularization) is the
most natural (and efficient) way.  While derivable from a quadratic penalty on
the parameter vector in the least squares optimization, this boils down to
"adjusting" singular values.

This is sometimes called "shrinkage" since large b.b are penalized.  Since the
ridge penalty is proportional to a flat b.b value, there is also usually a
desire to steer towards comparably scaled b_i via standardizing Xs & ys (e.g.
mean 0, variance 1) which `fitl` also supports.

The question remains how to pick the ridge parameter (or penalty coefficient).
`fitl` supports both manual (as a fraction of the largest singular value) and
automatic model selection via the `--sv` and `--xv` options.  Any negative value
for `sv` is flipped and used directly.  Meanwhhile `sv==0.0` implies a cross
validation approach by minimizing the score specified in `xv`.

Right now only PRESS/Leave-One-Out cross validation & Generalized Cross
Validation scores are supported.  Cross validation is splitting data into
estimation (or training) and prediction (or test) sets { and probably other
terms as well).  Both GCV & LOO/PRESS scores are efficiently computable from the
SVD solution, but some reasonably advocate chunkier splitting and `fitl` may
grow such an option. { Indeed, from a purely scientific perspective, "in the
wild" scant data is already often wildly "over (re-)used".  Simple K-fold cross
validation might perform better for "true out of sample". }

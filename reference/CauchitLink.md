# S7 Class for the Cauchit Link

Carries the cauchit transformation \\\eta = \tan(\pi(\theta - 1/2))\\ on
\\(0, 1)\\, the Cauchy quantile function, with inverse \\\theta = 1/2 +
\arctan(\eta)/\pi\\.

Its tails are far heavier than the logit's or the probit's, so an
extreme linear predictor moves the probability much less. That makes it
the choice when a few observations would otherwise drive the fit to a
boundary.

## Usage

``` r
CauchitLink(
  link_name = character(0),
  link_bounds = integer(0),
  link_params = NULL
)
```

## Arguments

- link_name:

  A character string naming the link, set by the constructor and shown
  by [`print()`](https://rdrr.io/r/base/print.html).

- link_bounds:

  A length-two numeric vector, the open interval the parameter lives in.
  Set by the constructor; see Value for this link's.

- link_params:

  A list of the link's own parameters, empty where it has none. Set by
  the constructor.

## Value

An S7 object of class `CauchitLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(0, 1)` and it carries no link
parameters.

## Methods

Ten methods are registered on this class:
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
and
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md),
and the four derivative orders in each direction,
[`dlinkfun()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkfun.md)
through
[`d4linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d4linkfun.md)
going out and
[`dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md)
through
[`d4linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d4linkinv.md)
coming back.
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
and
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
delegate to [`stats::qcauchy()`](https://rdrr.io/r/stats/Cauchy.html)
and [`stats::pcauchy()`](https://rdrr.io/r/stats/Cauchy.html), which
stay accurate in both tails. The eight derivatives come from a compiled
kernel, one call per order and direction.

## See also

[`cauchit_link()`](https://statmodels7.github.io/linkfunctions7/reference/cauchit_link.md),
the constructor users call.

# S7 Class for the InverseSq Link

Carries the inverse-square transformation \\\eta = 1/\theta^2\\ on \\(0,
\infty)\\, with inverse \\\theta = 1/\sqrt{\eta}\\.

It is the canonical link of the inverse Gaussian family. Its image is
\\(0, \infty)\\, so a negative linear predictor has no parameter behind
it, and like
[`inverse_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_link.md)
the map is decreasing.

## Usage

``` r
InverseSqLink(
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

An S7 object of class `InverseSqLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(0, Inf)` and it carries no link
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
coming back. The two directions are not the same map, \\1/\theta^2\\
going forward and \\\eta^{-1/2}\\ coming back, so the two sets of
derivatives are written out separately: integer falling factorials one
way, half-integer powers the other.

## See also

[`inverse_sq_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_sq_link.md),
the constructor users call.

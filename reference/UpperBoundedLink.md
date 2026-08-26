# S7 Class for an Upper Bounded Link

Carries the reflected log on \\(-\infty, u)\\: \\\eta = \log(u -
\theta)\\, with inverse \\\theta = u - e^{\eta}\\.

The reflection makes the map decreasing, so the odd-order derivatives
change sign against
[LowerBoundedLink](https://statmodels7.github.io/linkfunctions7/reference/LowerBoundedLink.md)'s
while the even ones do not. It is the mirror image of
[`LowerBoundedLink()`](https://statmodels7.github.io/linkfunctions7/reference/LowerBoundedLink.md),
the log of the distance below `upr`.

## Usage

``` r
UpperBoundedLink(
  link_name = character(0),
  link_bounds = integer(0),
  link_params = NULL,
  upr = integer(0)
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

- upr:

  The upper endpoint.

## Value

An S7 object of class `UpperBoundedLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(-Inf, upr)`, and its
`link_params` holds `upr`.

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
coming back. The reflection makes the map decreasing, so every inverse
derivative is the negative of
[`exp_floored()`](https://statmodels7.github.io/linkfunctions7/reference/exp_floored.md)
of \\\eta\\, and the forward ones are the log's read at \\\mathrm{upr} -
\theta\\, the odd orders carrying the sign the reflection introduces and
the even ones not.

## See also

[`bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md),
the constructor users call.

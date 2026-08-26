# S7 Class for the Sqrt Link

Carries the square-root transformation \\\eta = \sqrt{\theta}\\ on \\(0,
\infty)\\, with inverse \\\theta = \eta^2\\.

Its image is only \\(0, \infty)\\, so a negative linear predictor has no
parameter behind it. That is a property of the link, and
[`check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
reports the failing invertibility check as expected for exactly this
reason.

## Usage

``` r
SqrtLink(
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

An S7 object of class `SqrtLink`, inheriting from
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
coming back. The forward derivatives are half-integer falling
factorials. The inverse map is \\\eta^2\\, so its derivatives terminate:
the second is the constant two and the third and fourth are exactly
zero, both built by
[`const_like()`](https://statmodels7.github.io/linkfunctions7/reference/const_like.md)
so a missing value still propagates.

## See also

[`sqrt_link()`](https://statmodels7.github.io/linkfunctions7/reference/sqrt_link.md),
the constructor users call.

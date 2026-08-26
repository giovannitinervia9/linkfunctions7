# S7 Class for the Softplus Link

Carries the softplus transformation on \\(0, \infty)\\, whose inverse
\\\theta = \log(1 + e^{a\eta})/a\\ is a smooth approximation of
\\\max(0, \eta)\\ that sharpens as the scale \\a\\ grows.

It is the alternative to the log link for a positive parameter: the log
link maps a large negative \\\eta\\ to something indistinguishable from
zero, while the softplus approaches zero linearly and stays numerically
alive there. The scale is stored in `link_params`, so one class serves
every \\a\\.

## Usage

``` r
SoftplusLink(
  link_name = character(0),
  link_bounds = integer(0),
  link_params = NULL,
  a = integer(0)
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

- a:

  The scale parameter, strictly positive.

## Value

An S7 object of class `SoftplusLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
and carrying its three properties `link_name`, `link_bounds` and
`link_params`. Its `link_bounds` are `c(0, Inf)` and its `link_params`
holds `a`.

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
coming back. The softplus is an antiderivative of the logistic, so its
inverse-link derivatives are the logistic ones shifted a place,
\\h^{(k+1)} = a^k \sigma^{(k)}\\, and they reuse the same polynomials
the logit does. The forward set comes from a compiled kernel written in
\\u = -\mathrm{expm1}(-a\theta)\\, which is finite where the
\\\mathrm{expm1}(a\theta)\\ form overflows.

## See also

[`softplus_link()`](https://statmodels7.github.io/linkfunctions7/reference/softplus_link.md),
the constructor users call.

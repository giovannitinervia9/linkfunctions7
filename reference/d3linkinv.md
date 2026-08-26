# 3rd Derivative of an Inverse Link Function

The third derivative of the inverse link \\g^{-1}(\eta)\\ with respect
to the linear predictor. This is the direction a modeling routine
working on the unconstrained scale needs: it is the chain-rule factor
that carries a derivative of the log-likelihood from \\\theta\\ onto
\\\eta\\.

## Usage

``` r
d3linkinv(x, eta)
```

## Arguments

- x:

  An object of class `link`.

- eta:

  A numeric vector of linear predictors. Any finite value is admissible;
  the inverse link clamps its result strictly inside `x@link_bounds`
  before this derivative is taken.

## Value

A numeric vector of the same length as `eta`, missing wherever `eta` is.

## Details

Every link answers this generic. A link whose class registers no method
for it gets the base class's numerical one, which applies a single
central stencil to the highest order that link does supply analytically,
never a chain of lower-order differences.
[`link_fallback_orders()`](https://statmodels7.github.io/linkfunctions7/reference/link_fallback_orders.md)
says which orders of a given link are exact, and
[`check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
leaves a fallback order unchecked, since comparing it against a
difference of itself would agree however wrong the link is.

Call this generic directly in a hot loop.
[`linkderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md)
and
[`linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md)
route by order and so dispatch twice, once on themselves and once here,
which is about a third of the cost of the call.

## See also

[`linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md),
which routes to this generic by order, and
[`d3linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d3linkfun.md)
for the same order in the other direction.

## Examples

``` r
# Every derivative of exp is exp, so the log link's inverse gives the
# same number at every order.
d3linkinv(log_link(), 1) - exp(1)
#> [1] 0

# The logit's inverse derivatives are polynomials in theta. At eta = 0 the
# logistic is symmetric about 1/2, so its even-order derivatives vanish
# there while the first is the Bernoulli variance, 1/4.
d3linkinv(logit_link(), 0)
#> [1] -0.125
```

# 5th Derivative of an Inverse Link Function

The fifth derivative of the inverse link \\g^{-1}(\eta)\\ with respect
to the linear predictor.

## Usage

``` r
d5linkinv(x, eta)
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

This is the order a score-driven filter's outer curvature reaches. The
chain rule that carries a log-likelihood derivative of order \\k\\ from
\\\theta\\ onto \\\eta\\ contracts the family's components against the
partial Bell polynomials in \\h', \ldots, h^{(k)}\\, so an order-5
quantity on the unconstrained scale needs the fifth derivative of the
inverse link and nothing beyond it.

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
[`d5linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d5linkfun.md)
for the same order in the other direction.

## Examples

``` r
# Every derivative of exp is exp, so the log link's inverse gives the
# same number at every order.
d5linkinv(log_link(), 1) - exp(1)
#> [1] 0

# At eta = 0 the logistic is symmetric about 1/2, so its even-order
# derivatives vanish there while the odd ones do not.
d5linkinv(logit_link(), 0)
#> [1] 0.25
```

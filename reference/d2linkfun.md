# 2nd Derivative of a Link Function

The second derivative of the link \\g(\theta)\\ with respect to the
parameter, on the parameter scale. This is the direction a delta-method
standard error is carried in, from a variance on \\\theta\\ to one on
\\\eta\\.

## Usage

``` r
d2linkfun(x, theta)
```

## Arguments

- x:

  An object of class `link`.

- theta:

  A numeric vector of parameter values, inside `x@link_bounds`. A value
  outside gives `NaN` or `NA` according to the link, and nothing is
  thrown.

## Value

A numeric vector of the same length as `theta`, missing wherever `theta`
is.

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

[`linkderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md),
which routes to this generic by order, and
[`d2linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d2linkinv.md)
for the same order in the other direction.

## Examples

``` r
# The log link's forward derivatives are -1 / t^2, so at theta = 2:
d2linkfun(log_link(), 2) - (-1 / 2^2)
#> [1] 0

# Missingness propagates instead of being filled in.
d2linkfun(logit_link(), c(0.5, NA))
#> [1]  0 NA
```

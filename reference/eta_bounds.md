# The Range of Predictors a Link Admits

The image of the link's parameter bounds under
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md),
which is the set of predictors the inverse link is defined on. Where a
link maps onto the whole real line the answer is `c(-Inf, Inf)`; where
it does not, the two finite ends are the boundary of what a caller may
hand to
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md).

## Usage

``` r
eta_bounds(x)
```

## Arguments

- x:

  A
  [`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
  object.

## Value

A numeric vector of length two, sorted, with `-Inf` or `Inf` in either
position where that end is unbounded.

## Details

A link need not map onto the whole real line: the square root reaches
only the positive half, and so do
[`inverse_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_link.md),
[`inverse_sq_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_sq_link.md)
and
[`power_link()`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md)
at a positive exponent. The bounds are returned sorted, since a
decreasing link reverses them, and are infinite in the directions where
they cannot be established.

## Why a caller outside this package asks

The internal use is to keep a finite-difference grid inside the set the
inverse link is defined on, a stencil straying outside returning `NaN`
and so making a numerical derivative missing rather than inaccurate. The
use from outside is a different question with the same answer:
[link_bounds()](https://statmodels7.github.io/linkfunctions7/reference/link.md)
says what a link maps **onto**, and a consumer that carries an
unconstrained vector needs to know what it maps **from**.

Both ends matter, and only together. A family that reads a free vector
in \\\mathbb{R}^d\\ and applies an inverse link to each coordinate needs
the map to be defined and injective there, and a link with finite eta
bounds is neither:
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
of a square root link is even, so `-2` and `2` both give 4 and the round
trip returns the absolute value. Asking
`all(is.infinite(eta_bounds(link)))` is how such a family rejects one at
construction, where the message can name the link.

## See also

[`link_bounds_clamp()`](https://statmodels7.github.io/linkfunctions7/reference/link_bounds_clamp.md),
which keeps
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)'s
result strictly inside the parameter bounds, and
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md),
the map whose image this is.

## Examples

``` r
# A link from the whole real line onto the positive half of the theta axis.
eta_bounds(log_link())
#> [1] -Inf  Inf

# The square root reaches only the positive half of the eta axis, so its
# inverse is even there and the round trip returns the absolute value.
eta_bounds(sqrt_link())
#> [1]   0 Inf
linkfun(sqrt_link(), linkinv(sqrt_link(), -2))
#> [1] 2

# Which is the question a consumer carrying an unconstrained vector asks.
links <- list(log_link(), softplus_link(), logit_link(), sqrt_link(),
              inverse_link(), inverse_sq_link(), power_link(0.5))
vapply(links, function(l) all(is.infinite(eta_bounds(l))), logical(1))
#> [1]  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE
```

# The Body Shared by Every Numerical Fallback

Computes the order-`order` derivative of a link, in either direction, by
differentiating once the highest order the link supplies analytically.

## Usage

``` r
fallback_deriv(x, v, order, inverse)
```

## Arguments

- x:

  An object of class `link`.

- v:

  A numeric vector: \\\theta\\ going forward, \\\eta\\ coming back.

- order:

  The derivative order wanted, 1 to 5.

- inverse:

  Logical; `TRUE` for the inverse-link direction.

## Value

A numeric vector of the same length as `v`.

## Details

The whole design is in the two lines that pick `m` and `gap`: never
differentiate numerically more than the number of orders actually
missing. A link analytic to the second order asks for a first difference
to reach the third, not three; a link supplying nothing but
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
is the only case in which a fourth-order stencil is applied to the
function itself.

Note the recursion is only apparent. The base function is fetched
through
[`linkderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md),
which dispatches to the link's own method for an order it implements, so
the chain always terminates on analytic code, and never on another
fallback.

## Methods

The ten registrations on the base class
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
have this function as their whole body:
[`dlinkfun()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkfun.md)
through
[`d5linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d5linkfun.md)
going out and
[`dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md)
through
[`d5linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d5linkinv.md)
coming back, each passing its order and its direction. A link inherits
them for the orders it does not implement itself, so a link defined with
nothing but
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
and
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
can still answer every derivative generic. S7 requires a method's
formals to match the generic's, and the two directions name their
argument differently, so the eight wrappers are written out rather than
generated.

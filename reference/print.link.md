# Print Method for S7 Link Objects

Prints a link's name, the open interval its parameter lives in, and the
link parameters it carries, if any. Three lines at most, and two for a
link with no parameters of its own.

## Usage

``` r
# S3 method for class 'link'
print(x, ...)
```

## Arguments

- x:

  An object of class `link`.

- ...:

  Additional arguments passed to methods, currently unused.

## Value

`x`, invisibly. Called for the printed output.

## Details

The domain is shown as the open interval it is, so a probability link
reads `(0, 1)` and never `[0, 1]`: no link ever returns an endpoint, and
[`link_bounds_clamp()`](https://statmodels7.github.io/linkfunctions7/reference/link_bounds_clamp.md)
is what keeps that true in double precision.

The parameter line appears only for a link that has parameters, and
names them, so `power(lambda=2)` and `bounded(lwr=0, upr=10)` report the
values they were constructed with.

## Examples

``` r
print(logit_link())
#> S7 Link Object: logit
#>   - Parameter domain (theta): (0, 1)

# links carrying parameters report them too
print(power_link(2))
#> S7 Link Object: power(lambda=2)
#>   - Parameter domain (theta): (0, Inf)
#>   - Link parameters: lambda = 2
print(bounded_link(0, 10))
#> S7 Link Object: bounded(lwr=0, upr=10)
#>   - Parameter domain (theta): (0, 10)
#>   - Link parameters: lwr = 0, upr = 10
```

# The five logistic derivative polynomials are written twice: in R, as
# logistic_deriv(), and in C++, as the static inline logistic_poly() that the
# three links reach through their compiled kernels. Nothing called the R one,
# so nothing held the two together; this is that test, and it is what makes
# logistic_deriv() the R statement of the kernel rather than dead code.
#
# The comparison carries a tolerance rather than asking expect_identical().
# Both forms are Horner, so both contain multiply-adds, and a compiler is free
# to contract one into an FMA and drop the intermediate rounding. Measured on
# this machine the two agree to the bit at every order; that is a property of
# one compiler and not of the arithmetic, and the toolkit has three recorded
# occasions of an identity assertion over compiled floating point going red on
# one platform alone. The tolerance below is tight enough that a reordered or
# mistranscribed polynomial fails it by orders of magnitude.

grid <- seq(0.001, 0.999, length.out = 501)

test_that("the R polynomials match the compiled ones at every order", {
  for (k in 1:5) {
    expect_equal(
      logistic_deriv(grid, k), lk_logistic_poly_cpp(grid, k),
      tolerance = 1e-15
    )
  }
})

test_that("a mistranscribed polynomial would fail the comparison", {
  # the negative control: order 3 with one coefficient wrong by a unit
  wrong <- grid * (1 - grid) * (1 + grid * (-5 + 6 * grid))
  expect_gt(max(abs(wrong - lk_logistic_poly_cpp(grid, 3L))), 1e-3)
})

test_that("the three links reach the polynomials the page says they do", {
  z <- seq(-6, 6, length.out = 401)

  # the logit uses them as they stand
  p <- stats::plogis(z)
  for (k in 1:5) {
    expect_equal(linkinvderiv(logit_link(), z, k), logistic_deriv(p, k),
                 tolerance = 1e-15)
  }

  # a doubly bounded link scales them by the interval width
  w <- 3
  bl <- bounded_link(lwr = 1, upr = 1 + w)
  for (k in 1:5) {
    expect_equal(linkinvderiv(bl, z, k), w * logistic_deriv(p, k),
                 tolerance = 1e-15)
  }

  # the softplus is an antiderivative of the logistic, so it uses them one
  # order down, scaled by a power of its own steepness
  a <- 2
  q <- stats::plogis(a * z)
  for (k in 2:5) {
    expect_equal(linkinvderiv(softplus_link(a), z, k),
                 a^(k - 1) * logistic_deriv(q, k - 1L), tolerance = 1e-15)
  }
})

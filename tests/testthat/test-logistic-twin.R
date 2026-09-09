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
# one compiler and not of the arithmetic, and the toolkit has four recorded
# occasions of an identity assertion over compiled floating point going red on
# one platform alone.
#
# The tolerance is taken against each polynomial's own SCALE and not pointwise,
# because these polynomials pass through zero inside the grid and a relative
# comparison there divides by nearly nothing. expect_agrees_on_scale() in
# helper-links.R carries the measurement and the margin.

grid <- seq(0.001, 0.999, length.out = 501)

test_that("the R polynomials match the compiled ones at every order", {
  for (k in 1:5) {
    expect_agrees_on_scale(
      logistic_deriv(grid, k), lk_logistic_poly_cpp(grid, k),
      label = sprintf("logistic_deriv at order %d", k)
    )
  }
})

test_that("a mistranscribed polynomial would fail the comparison", {
  # the negative control: order 3 with one coefficient wrong by a unit
  wrong <- grid * (1 - grid) * (1 + grid * (-5 + 6 * grid))
  right <- lk_logistic_poly_cpp(grid, 3L)
  expect_gt(max(abs(wrong - right)), 1e-3)

  # and it must fail the comparison the tests above actually run, which is what
  # says the scale tolerance kept its teeth: the gap is 1.19 against that
  # order's own scale, where the arithmetic's own disagreement is 1e-13
  expect_failure(expect_agrees_on_scale(wrong, right))
  expect_gt(max(abs(wrong - right)) / max(abs(right)), 1)
})

test_that("a disagreement the size a runner produces still passes", {
  # The other edge of the same margin. On the ARM64 macOS runner, whose
  # compiler contracts a multiply-add where the one here does not, the worst
  # absolute gap between the two routes over 965 differing points was 9.1e-15.
  # A tolerance that failed on that would make the suite red on one platform
  # for arithmetic that is not a defect, so the check runs in both directions.
  right <- lk_logistic_poly_cpp(grid, 3L)
  expect_agrees_on_scale(right + 9.1e-15, right)
})

test_that("the three links reach the polynomials the page says they do", {
  z <- seq(-6, 6, length.out = 401)

  # the logit uses them as they stand
  p <- stats::plogis(z)
  for (k in 1:5) {
    expect_agrees_on_scale(linkinvderiv(logit_link(), z, k),
                           logistic_deriv(p, k),
                           label = sprintf("logit at order %d", k))
  }

  # a doubly bounded link scales them by the interval width
  w <- 3
  bl <- bounded_link(lwr = 1, upr = 1 + w)
  for (k in 1:5) {
    expect_agrees_on_scale(linkinvderiv(bl, z, k), w * logistic_deriv(p, k),
                           label = sprintf("doubly bounded at order %d", k))
  }

  # the softplus is an antiderivative of the logistic, so it uses them one
  # order down, scaled by a power of its own steepness
  a <- 2
  q <- stats::plogis(a * z)
  for (k in 2:5) {
    expect_agrees_on_scale(linkinvderiv(softplus_link(a), z, k),
                           a^(k - 1) * logistic_deriv(q, k - 1L),
                           label = sprintf("softplus at order %d", k))
  }
})

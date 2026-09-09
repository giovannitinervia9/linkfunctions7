# Every link the package ships, in one place. Lives in a helper rather than in
# a test file because more than one file needs it, and testthat gives each test
# file its own environment: a function defined at the top of one is not visible
# from another.

all_links <- function() {
  list(
    identity   = identity_link(),
    log        = log_link(),
    logit      = logit_link(),
    probit     = probit_link(),
    cloglog    = cloglog_link(),
    loglog     = loglog_link(),
    cauchit    = cauchit_link(),
    rhobit     = rhobit_link(),
    sqrt       = sqrt_link(),
    inverse    = inverse_link(),
    inverse_sq = inverse_sq_link(),
    power_2    = power_link(2),
    power_half = power_link(0.5),
    softplus_1 = softplus_link(1),
    softplus_3 = softplus_link(3),
    lower_b    = bounded_link(lwr = 2),
    upper_b    = bounded_link(upr = 5),
    both_b     = bounded_link(lwr = 2, upr = 5)
  )
}

# Two floating-point routes to the same polynomial, compared against that
# polynomial's own scale rather than pointwise.
#
# A pointwise relative comparison is the wrong instrument here, and the reason
# is the polynomials rather than the arithmetic. The logistic derivative
# polynomials pass through zero inside any grid spanning the unit interval:
# orders 2 and 4 reach exactly zero and orders 3 and 5 come within 3e-4 of it,
# with two sign changes at orders 2 and 3 and four at orders 4 and 5. Near a
# root the denominator of a relative comparison is nearly nothing, so one
# contracted multiply-add reads there as a relative error of 2.3e-12 while the
# same absolute difference elsewhere reads as 1e-16.
#
# Measured on the ARM64 macOS runner, whose compiler contracts a multiply-add
# where the one used here does not, the worst ABSOLUTE disagreement over 965
# differing points is 9.1e-15. Against each order's own scale that is between
# 3.6e-14 and 9.5e-14, and one mistranscribed coefficient costs 1.19 on the
# same scale. The default below sits ten times above the worst measured and
# twelve orders below the defect the comparison exists to catch.
#
# Where the expected vector is identically zero the scale is zero and the
# comparison becomes exact, which is the right reading of that case.
expect_agrees_on_scale <- function(object, expected, tolerance = 1e-12,
                                   label = NULL) {
  testthat::expect_lt(
    max(abs(object - expected)), tolerance * max(abs(expected)),
    label = label
  )
}

#' @title S7 Class for the Sqrt Link
#'
#' @description
#' Carries the square-root transformation \eqn{\eta = \sqrt{\theta}} on
#' \eqn{(0, \infty)}, with inverse \eqn{\theta = \eta^2}.
#'
#' Its image is only \eqn{(0, \infty)}, so a negative linear predictor has no
#' parameter behind it. That is a property of the link, and [check_link()]
#' reports the failing invertibility check as expected for exactly this
#' reason.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval the
#'   parameter lives in. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#'
#' @return An S7 object of class `SqrtLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(0, Inf)` and it carries no link parameters.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. The forward derivatives are half-integer falling factorials. The
#' inverse map is \eqn{\eta^2}, so its derivatives terminate: the second is
#' the constant two and the third and fourth are exactly zero, both built by
#' [const_like()] so a missing value still propagates.
#'
#' @aliases linkfun.SqrtLink
#' @aliases linkinv.SqrtLink
#' @aliases dlinkfun.SqrtLink
#' @aliases d2linkfun.SqrtLink
#' @aliases d3linkfun.SqrtLink
#' @aliases d4linkfun.SqrtLink
#' @aliases d5linkfun.SqrtLink
#' @aliases dlinkinv.SqrtLink
#' @aliases d2linkinv.SqrtLink
#' @aliases d3linkinv.SqrtLink
#' @aliases d4linkinv.SqrtLink
#' @aliases d5linkinv.SqrtLink
#'
#' @seealso [sqrt_link()], the constructor users call.
#' @keywords internal
SqrtLink <- S7::new_class(
  name = "SqrtLink",
  parent = link
)

# --- Methods for SqrtLink ---

# Forward and inverse link functions
S7::method(linkfun, SqrtLink) <- function(x, theta) sqrt(theta)
S7::method(linkinv, SqrtLink) <- function(x, eta) eta^2

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, SqrtLink) <- function(x, theta) {
  1 / (2 * sqrt(theta))
}
S7::method(d2linkfun, SqrtLink) <- function(x, theta) {
  -1 / (4 * (theta^1.5))
}
S7::method(d3linkfun, SqrtLink) <- function(x, theta) {
  3 / (8 * (theta^2.5))
}
S7::method(d4linkfun, SqrtLink) <- function(x, theta) {
  -15 / (16 * (theta^3.5))
}
S7::method(d5linkfun, SqrtLink) <- function(x, theta) {
  105 / (32 * (theta^4.5))
}

# Exact analytical derivatives of the inverse link function (wrt eta)
# 3rd and 4th derivatives uniquely vanish to exactly 0 for this quadratic form.
S7::method(dlinkinv, SqrtLink) <- function(x, eta) 2 * eta
S7::method(d2linkinv, SqrtLink) <- function(x, eta) const_like(eta, 2)
S7::method(d3linkinv, SqrtLink) <- function(x, eta) const_like(eta, 0)
S7::method(d4linkinv, SqrtLink) <- function(x, eta) const_like(eta, 0)
S7::method(d5linkinv, SqrtLink) <- function(x, eta) const_like(eta, 0)

#' @title The Square Root Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The square-root link \eqn{\eta = \sqrt{\theta}} on \eqn{(0, \infty)};
#' its image is \eqn{(0, \infty)}.
#' @details
#' The Square Root link is mathematically defined as \eqn{\eta = \sqrt{\theta}}.
#' Consequently, the inverse link is derived as \eqn{\theta = \eta^2}.
#'
#' Unlike the Log link, this transformation allows \eqn{\theta} to reach 0 exactly.
#' While the inverse function (\eqn{\eta^2}) is mathematically valid for negative 
#' values of \eqn{\eta}, in the specific context of this link function, the linear 
#' predictor \eqn{\eta} is typically constrained to be non-negative. This restriction 
#' preserves a strictly one-to-one mapping with \eqn{\theta}.
#'
#' The strict mathematical domain for \eqn{\theta} is `c(0, Inf)`.
#'
#' @return An S7 object of class `SqrtLink` (inheriting from `link`) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @examples
#' lk <- sqrt_link()
#' lk
#'
#' theta <- c(0.25, 1, 4)
#' eta <- linkfun(lk, theta)
#' eta
#' linkinv(lk, eta)
#'
#' # the inverse is a quadratic, so the third and fourth derivatives vanish
#' d2linkinv(lk, c(1, 2))
#' d3linkinv(lk, c(1, 2))
#'
#' # the same link as the power family at lambda = 1/2
#' linkfun(power_link(0.5), 4)
#'
#' @seealso [link()], [power_link()], [log_link()]
#' @export
sqrt_link <- function() {
  SqrtLink(
    link_name = "sqrt",
    link_bounds = c(0, Inf),
    
    # The standard square root link requires no additional parameters
    link_params = NULL
  )
}

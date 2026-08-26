#' @title S7 Class for the Rhobit Link
#'
#' @description
#' Carries the rhobit transformation
#' \eqn{\eta = \mathrm{atanh}(\theta) = \tfrac{1}{2}\log((1+\theta)/(1-\theta))}
#' on \eqn{(-1, 1)}, with inverse \eqn{\theta = \tanh(\eta)}.
#'
#' This is Fisher's z, the natural chart for a correlation: it carries the open
#' interval onto the whole line, so an optimizer moving freely in \eqn{\eta}
#' never proposes a correlation outside its range.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval the
#'   parameter lives in. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#'
#' @return An S7 object of class `RhobitLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(-1, 1)` and it carries no link parameters.
#'
#' @section Methods:
#' Ten methods are registered on this class: [linkfun()] and [linkinv()],
#' and the four derivative orders in each direction, [dlinkfun()] through
#' [d4linkfun()] going out and [dlinkinv()] through [d4linkinv()] coming
#' back. `linkfun()` and `linkinv()` are `atanh()` and `tanh()`. The eight
#' derivatives come from a compiled kernel, one call per order and
#' direction.
#'
#' @aliases linkfun.RhobitLink
#' @aliases linkinv.RhobitLink
#' @aliases dlinkfun.RhobitLink
#' @aliases d2linkfun.RhobitLink
#' @aliases d3linkfun.RhobitLink
#' @aliases d4linkfun.RhobitLink
#' @aliases dlinkinv.RhobitLink
#' @aliases d2linkinv.RhobitLink
#' @aliases d3linkinv.RhobitLink
#' @aliases d4linkinv.RhobitLink
#'
#' @seealso [rhobit_link()], the constructor users call.
#' @keywords internal
RhobitLink <- S7::new_class(
  name = "RhobitLink",
  parent = link
)

# --- Methods for RhobitLink ---

# Forward and inverse link functions
S7::method(linkfun, RhobitLink) <- function(x, theta) atanh(theta)
S7::method(linkinv, RhobitLink) <- function(x, eta) tanh(eta)

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, RhobitLink) <- function(x, theta) lk_rhobit_fwd_cpp(theta, 1L)
S7::method(d2linkfun, RhobitLink) <- function(x, theta) lk_rhobit_fwd_cpp(theta, 2L)
S7::method(d3linkfun, RhobitLink) <- function(x, theta) lk_rhobit_fwd_cpp(theta, 3L)
S7::method(d4linkfun, RhobitLink) <- function(x, theta) lk_rhobit_fwd_cpp(theta, 4L)

# Exact analytical derivatives of the inverse link function (wrt eta)
# All of these are evaluated as polynomials in t = tanh(eta)
S7::method(dlinkinv, RhobitLink) <- function(x, eta) lk_rhobit_inv_cpp(eta, 1L)
S7::method(d2linkinv, RhobitLink) <- function(x, eta) lk_rhobit_inv_cpp(eta, 2L)
S7::method(d3linkinv, RhobitLink) <- function(x, eta) lk_rhobit_inv_cpp(eta, 3L)
S7::method(d4linkinv, RhobitLink) <- function(x, eta) lk_rhobit_inv_cpp(eta, 4L)

#' @title The Rhobit (Fisher's z) Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The rhobit link
#' \eqn{\eta = \mathrm{atanh}(\theta) = \log((1+\theta)/(1-\theta))/2}
#' on \eqn{(-1, 1)}, Fisher's z; the natural link for a correlation.
#' @details
#' The Rhobit link is defined mathematically using the inverse hyperbolic tangent function:
#' \eqn{\eta = \text{arctanh}(\theta) = \frac{1}{2} \log\left(\frac{1 + \theta}{1 - \theta}\right)}.
#'
#' The inverse link is the hyperbolic tangent function:
#' \eqn{\theta = \tanh(\eta) = \frac{\exp(2\eta) - 1}{\exp(2\eta) + 1}}.
#'
#' The valid mathematical domain of \eqn{\theta} is exactly `c(-1, 1)`.
#'
#' @return An S7 object of class `RhobitLink` (inheriting from `link`) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @examples
#' lk <- rhobit_link()
#' lk
#'
#' # built for parameters constrained to (-1, 1), such as a correlation
#' rho <- c(-0.9, 0, 0.9)
#' eta <- linkfun(lk, rho)    # Fisher's z
#' eta
#' linkinv(lk, eta)
#'
#' # the inverse is tanh, so its first derivative is the squared sech
#' dlinkinv(lk, 0)
#'
#' check_link(lk)
#'
#' @seealso [link()], [logit_link()]
#' @export
rhobit_link <- function() {
  RhobitLink(
    link_name = "rhobit",
    link_bounds = c(-1, 1),
    link_params = NULL
  )
}

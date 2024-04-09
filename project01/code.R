#' Orlagh Keane, S2084384
#' Add your own function definitions on this file.

#' Log-Exponential density
#'
#' Compute the density or log-density for a Log-Exponential (LogExp)
#' distribution
#'
#' @param x vector of quantiles
#' @param rate vector of rates
#' @param log logical; if TRUE, the log-density is returned

dlogexp <- function(x, rate = 1, log = FALSE) {
  result <- log(rate) + x - rate * exp(x)
  if (!log) {
    exp(result)
  }
  result
}

#' Log-Sum-Exp
#'
#' Convenience function for computing log(sum(exp(x))) in a
#' numerically stable manner
#'
#' @param x numerical vector

log_sum_exp <- function(x) {
  max_x <- max(x, na.rm = TRUE)
  max_x + log(sum(exp(x - max_x)))
}


#' wquantile 
#'
#' Calculates empirical sample quantiles with optional weights, for given probabilities. 
#' Like in quantile(), the smallest observation corresponds to a probability of 0 and the largest to a probability of 1. 
#' Interpolation between discrete values is done when type=7, as in quantile(). 
#' Use type=1 to only generate quantile values from the raw input samples.
#'
#' @param x numeric vector whose sample quantiles are wanted
#' NA and NaN values are not allowed in numeric vectors unless na.rm is TRUE
#' @param probs numeric vector of probabilities with values in [0,1]
#' @param na.rm logical; if true, any NA and NaN's are removed from x before the quantiles are computed
#' @param type numeric, 1 for no interpolation, or 7, for interpolated quantiles. Default is 7
#' @param weights	 numeric vector of non-negative weights, the same length as x, or NULL. 
#' The weights are normalised to sum to 1. 
#' If NULL, then wquantile(x) behaves the same as quantile(x), with equal weight for each sample value

wquantile <- function (x, probs = seq(0, 1, 0.25), na.rm = FALSE, type = 7, 
                       weights = NULL, ...) 
{
  if (is.null(weights) || (length(weights) == 1)) {
    weights <- rep(1, length(x))
  }
  stopifnot(all(weights >= 0))
  stopifnot(length(weights) == length(x))
  if (length(x) == 1) {
    return(rep(x, length(probs)))
  }
  n <- length(x)
  q <- numeric(length(probs))
  reorder <- order(x)
  weights <- weights[reorder]
  x <- x[reorder]
  wecdf <- pmin(1, cumsum(weights)/sum(weights))
  if (type == 1) {
  }
  else {
    weights2 <- (weights[-n] + weights[-1])/2
    wecdf2 <- pmin(1, cumsum(weights2)/sum(weights2))
  }
  for (pr_idx in seq_along(probs)) {
    pr <- probs[pr_idx]
    if (pr <= 0) {
      q[pr_idx] <- x[1]
    }
    else if (pr >= 1) {
      q[pr_idx] <- x[n]
    }
    else {
      if (type == 1) {
        j <- 1 + pmax(0, pmin(n - 1, sum(wecdf <= pr)))
        q[pr_idx] <- x[j]
      }
      else {
        j <- 1 + pmax(0, pmin(n - 2, sum(wecdf2 <= pr)))
        g <- (pr - c(0, wecdf2)[j])/(wecdf2[j] - c(0, 
                                                   wecdf2)[j])
        q[pr_idx] <- (1 - g) * x[j] + g * x[j + 1]
      }
    }
  }
  q
}

#' Compute empirical weighted cumulative distribution
#'
#' Version of `ggplot2::stat_ecdf` that adds a `weights` property for each
#' observation, to produce an empirical weighted cumulative distribution function.
#' The empirical cumulative distribution function (ECDF) provides an alternative
#' visualisation of distribution. Compared to other visualisations that rely on
#' density (like [geom_histogram()]), the ECDF doesn't require any
#' tuning parameters and handles both continuous and discrete variables.
#' The downside is that it requires more training to accurately interpret,
#' and the underlying visual tasks are somewhat more challenging.
#'
# @inheritParams layer
# @inheritParams geom_point
#' @param na.rm If `FALSE` (the default), removes missing values with
#'    a warning.  If `TRUE` silently removes missing values.
#' @param n if NULL, do not interpolate. If not NULL, this is the number
#'   of points to interpolate with.
#' @param pad If `TRUE`, pad the ecdf with additional points (-Inf, 0)
#'   and (Inf, 1)
#' @section Computed variables:
#' \describe{
#'   \item{x}{x in data}
#'   \item{y}{cumulative density corresponding x}
#' }
#' @seealso wquantile
#' @export
#' @examples
#' library(ggplot2)
#'
#' n <- 100
#' df <- data.frame(
#'   x = c(rnorm(n, 0, 10), rnorm(n, 0, 10)),
#'   g = gl(2, n),
#'   w = c(rep(1/n, n), sort(runif(n))^sqrt(n))
#' )
#' ggplot(df, aes(x, weights = w)) + stat_ewcdf(geom = "step")
#'
#' # Don't go to positive/negative infinity
#' ggplot(df, aes(x, weights = w)) + stat_ewcdf(geom = "step", pad = FALSE)
#'
#' # Multiple ECDFs
#' ggplot(df, aes(x, colour = g, weights = w)) + stat_ewcdf()
#' ggplot(df, aes(x, colour = g, weights = w)) +
#'   stat_ewcdf() +
#'   facet_wrap(vars(g), ncol = 1)

stat_ewcdf <- function(mapping = NULL, data = NULL,
                       geom = "step", position = "identity",
                       ...,
                       n = NULL,
                       pad = TRUE,
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatEwcdf,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      n = n,
      pad = pad,
      na.rm = na.rm,
      ...
    )
  )
}


#' @title StatEwcdf ggproto object
#' @name StatEwcdf
#' @rdname StatEwcdf
#' @aliases StatEwcdf
#' @format NULL
#' @usage NULL
#' @export
#' @importFrom ggplot2 aes after_stat has_flipped_aes Stat
NULL

StatEwcdf <- ggplot2::ggproto(
  "StatEwcdf", ggplot2::Stat,
  required_aes = c("x|y", "weights"),
  dropped_aes = c("weights"),     
  
  default_aes = ggplot2::aes(y = ggplot2::after_stat(y)),
  
  setup_params = function(data, params) {
    params$flipped_aes <-
      ggplot2::has_flipped_aes(data,
                               params,
                               main_is_orthogonal = FALSE,
                               main_is_continuous = TRUE)
    
    has_x <- !(is.null(data$x) && is.null(params$x))
    has_y <- !(is.null(data$y) && is.null(params$y))
    if (!has_x && !has_y) {
      rlang::abort("stat_ewcdf() requires an x or y aesthetic.")
    }
    has_weights <- !(is.null(data$weights) && is.null(params$weights))
    #    if (!has_weights) {
    #      rlang::abort("stat_ewcdf() requires a weights aesthetic.")
    #    }
    
    params
  },
  
  compute_group = function(data, scales, n = NULL, pad = TRUE, flipped_aes = FALSE) {
    data <- flip_data(data, flipped_aes)
    # If n is NULL, use raw values; otherwise interpolate
    if (is.null(n)) {
      x <- unique(data$x)
    } else {
      x <- seq(min(data$x), max(data$x), length.out = n)
    }
    
    if (pad) {
      x <- c(-Inf, x, Inf)
    }
    if (is.null(data$weights)) {
      data_ecdf <- ecdf(data$x)(x)
    } else {
      data_ecdf <-
        spatstat.geom::ewcdf(
          data$x,
          weights = data$weights / sum(abs(data$weights)) 
        )(x)
    }
    
    df_ecdf <- vctrs::new_data_frame(list(x = x, y = data_ecdf), n = length(x))
    df_ecdf$flipped_aes <- flipped_aes
    ggplot2::flip_data(df_ecdf, flipped_aes)
  }
)



# MY ANSWERS

# QUESTION 2

#' Negative Log Likelihood
#'
#' Calculate and returns the negated log likelihood for the specified model
#'
#' @param beta model parameters
#' @param data data frame containing required variables
#' @param model A or B
#' 
#' @return value for the neg log like based on model
neg_log_like <- function(beta, data, model) {
  xi <- data$CAD_Weight
  yi <- data$Actual_Weight
  if(model == "A") {
    sd <- sqrt(exp(beta[3] + beta[4] * xi))
  } else if(model == "B") {
    sd <- sqrt(exp(beta[3]) + exp(beta[4]) * xi^2)
  }
  mean <- beta[1] + beta[2] * xi
  return(-sum(dnorm(yi, mean , sd, log = TRUE)))
}


#' Data Optimised Estimation
#'
#' Uses optim and neg_log_like to estimate the two models and return the best 
#' set of parameters and the estimate of the Hessian Matrix at the solution
#'
#' @param data data frame containing required variables
#' @param model A or B
#' 
#' @return  the two optimal model params and its Hessian
filament1_estimate <- function(data, model) {
  # Define initial values for beta
  if(model == "A") {
    initials <- c(-0.1, 1.07, -2, 0.05)
  } else if(model == "B") {
    initials <- c(-0.15, 1.07, -13.5, -6.5)
  }
  
  # Run optimization
  fit <- optim(par = initials, 
               fn = neg_log_like, 
               data = data,
               model = model,
               method = "BFGS", 
               hessian = TRUE)
  
  # Return best set of parameters found and estimate of the Hessian
  return(list(parameters = fit$par, hessian = fit$hessian))
}








# Question 3



#' Log-density for prior distribution
#'
#' Uses dnorm and dlogexp to evaluate the log of the joint prior density for the 
#' four theta parameters
#'
#' @param theta parameter vector
#' @param params vector of γ parameters
#' 
#' @return sum of the log prior densities
log_prior_density <- function(theta, params) { 
  result <-  dnorm(theta[1], mean = 0, sd = sqrt(params[1]), log = TRUE) +
    dnorm(theta[2], mean = 1, sd = sqrt(params[2]), log = TRUE) +
    dlogexp(theta[3], rate = params[3], log = TRUE) +
    dlogexp(theta[4], rate = params[4], log = TRUE)
  return(result)
}

#' Observation log-likelihood
#'
#' Uses dnorm to evauluate the log of the post density
#'
#' @param theta parameter vector
#' @param x models x - CAD weight
#' @param y models y - actual weight
#' 
#' @return sum of the logs of the prior densities
log_like <- function(theta, x, y) { 
  priors <- dnorm(y,
   mean = theta[1] + theta[2] * x,
   sd = sqrt(exp(theta[3]) + exp(theta[4]) * x^2),
   log = TRUE
  )
  return(sum(priors))
}


#' Log-density for the posterior distribution
#'
#' Uses log_prior_density and log_like to evaluate the observation 
#' log-likelihood for the model
#'
#' @param theta parameter vector
#' @param x models x - CAD weight
#' @param y models y - actual weight
#' @param params vector of γ parameters
#' 
#' @return log density for post distribution of theta
log_posterior_density <- function(theta, x, y, params) { 
  density <- log_prior_density(theta, params) +
    log_like(theta, x, y)
  return(density)
}

#' Posterior mode
#'
#' Uses log_posterior_density to fond mode mu of the log-posterior-density
#' evaluates the Hessian at the mode, as well as the negated Hessian S
#'
#' @param theta_start parameter vector
#' @param x models x - CAD weight
#' @param y models y - actual weight
#' @param params vector of γ parameters
#' 
#' @return list with elements mode, hessian, S
posterior_mode <- function(theta_start, x, y, params) { 
  opt <- optim(theta_start, 
               log_posterior_density,
                x = x, 
                y = y, 
                params = params,
                control = list(fnscale = -1),
                hessian = TRUE
                )
  elements_list <- (list(
    mode = opt$par,
    hessian = opt$hessian,
    S = solve(-opt$hessian)
  )) 
  return(elements_list)
}

#' Importance sampler
#'
#' Uses log_posterior_density to fond mode mu of the log-posterior-density
#' evaluates the Hessian at the mode, aswell as the negated Hessian S
#'
#' @param N sample size
#' @param mu 
#' @param S negated Hessian matrix
#' @param x models x - CAD weight
#' @param y models y - actual weight
#' @param params vector of γ parameters
#' 
#' @return list with elements mode, hessian, S
do_importance <- function(N, mu, S, x, y, params = c(1, 1, 1, 1)) { 
  imp_sample <- mvtnorm::rmvnorm(N, mean = mu, sigma = S)
  log_weights <- rep(0, N)
  # calc the log weights for each
  for (i in 1:N) {
    log_weights[i] <-
      log_posterior_density(imp_sample[i, ], x = x, y = y, params = params) -
      mvtnorm::dmvnorm(imp_sample[i, ], mean = mu, sigma = S, log = TRUE)
  }
  # then normalise log weight and convert thetas 
  log_weights <- log_weights - log_sum_exp(log_weights)
  imp_sample[, 3:4] <- exp(imp_sample[, 3:4])
  # change col names
  colnames(imp_sample) <- c("beta1", "beta2", "beta3", "beta4")
  # Combine samples with log-weights into a data frame
  result <- cbind(as.data.frame(imp_sample), data.frame(log_weights = log_weights))
  return(result)
}




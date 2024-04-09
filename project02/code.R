#' Orlagh Keane, s2084384
#' Add your own function definitions on this file.

#' neg_log_lik
#
#' @description Evaluate the negated log-likelihood for model A and B
#' @param beta A vector with the beta parameters
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' Must have columns `CAD_Weight` and `Actual_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model

neg_log_lik <- function(beta, data, model){
  
  mu <- beta[1] + beta[2]*data[["CAD_Weight"]]
  
  # distinguish between the two models to find the particular standard deviation for the betas
  if(model == "A") {
    sigma <- sqrt(exp(beta[3] + beta[4]*data[["CAD_Weight"]]))
  }else{
    sigma <- sqrt(exp(beta[3])+exp(beta[4]) * (data[["CAD_Weight"]]^2))
  }
  - sum(dnorm(data[["Actual_Weight"]],
              mean = mu,
              sd=sigma,
              log = TRUE))
  
}

#' filament_estimate
#
#' @description Estimate filament models with different variance structure
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' Must have columns `CAD_Weight` and `Actual_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @return An estimation object suitable for use with [filament1_predict()]

filament1_estimate <- function(data, model) {
  model <- match.arg(model, c("A", "B"))
  if (model == "A") {
    beta_start <- c(-0.1, 1.07, -2, 0.05)
  } else {
    beta_start <- c(-0.15, 1.07, -13.5, -6.5)
  }
  opt <- optim(beta_start,
               neg_log_lik,
               data = data,
               model = model,
               hessian = TRUE,
               method = "Nelder-Mead",
               control = list(maxit = 5000)
  )
  fit <- list(
    model = model,
    par = opt$par,
    hessian = opt$hessian
  )
  class(fit) <- c("filament1_estimate", "list")
  fit
}

#' filament1_aux_EV
#' 
#' @description Evaluate the expectation and variance for model A and B
#' @param beta A vector with the beta parameters
#' @param data A `data.frame` containing the required predictors, including `CAD_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @param Sigma_beta : If not NULL, an estimate of the covariance matrix for
#                 the uncertainty of estimated betas
#' @return A list with four elements:
#     E : E(y|beta,x)
#     V : Var(y|beta,x)
#     VE : Var(E(y|beta,x)|x) or NULL
#     EV : E(Var(y|beta,x)|x) or NULL

filament1_aux_EV <- function(beta, data, model = c("A", "B"),
                             Sigma_beta = NULL) {
  
  model <- match.arg(model)
  if (model == "A") {
    
    ZE.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZV.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZE = cbind(ZE.0, ZV.0 * 0) 
    ZV = cbind(ZE.0 * 0, ZV.0)
    
    VE <- EV <- NULL
    if (!is.null(Sigma_beta)) {
      # E(Var(y|beta,x)|x)
      EV <- exp(ZV %*% beta + rowSums(ZV * (ZV %*% Sigma_beta)) / 2)
      # Var(E(y|beta,x)|x)
      VE <- rowSums(ZE * (ZE %*% Sigma_beta))
    }
    out <- list(
      E = ZE %*% beta,
      V = exp(ZV %*% beta),
      VE = VE,
      EV = EV
    )
  } else {
    
    ZE.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZV.0 <- model.matrix( ~ 1 + I(CAD_Weight^2), data = data)
    ZE = cbind(ZE.0, ZV.0 * 0) 
    ZV = cbind(ZE.0 * 0, ZV.0)
    
    VE <- EV <- NULL
    if (!is.null(Sigma_beta)) {
      # E(Var(y|beta,x)|x)
      # (pmin: Ignore large Sigma_beta values)
      EV <- ZV %*% exp(beta + pmin(0.5^2, diag(Sigma_beta)) / 2)
      # Var(E(y|beta,x)|x)
      VE <- rowSums(ZE * (ZE %*% Sigma_beta))
    }
    out <- list(
      E = ZE %*% beta,
      V = ZV %*% exp(beta),
      VE = VE,
      EV = EV
    )
  }
  out
}

#' filament1_predict
#' 
#' @description Uses estimated model to predict new observations
#' @param est filament1_estimate()
#' @param newdata A `data.frame` to predict
#' @param model Must be A or B
#' @return A data.frame of `pred_mean`, `pred_sd`, `lwr`, and `upr`, from the prediction for each row of `newdata`.
#' 
filament1_predict <- function(est, newdata, model=c("A", "B")) {
  theta <- est$par
  Sigma_theta <- solve(est$hessian)
  aux_EV <- filament1_aux_EV(beta = theta, data = newdata, model = model, Sigma_beta = Sigma_theta)
  pred_mean <- aux_EV$E
  pred_var <- aux_EV$V
  pred_sd <- (aux_EV$EV + aux_EV$VE)^0.5
  q <- qt(1 - 0.05 / 2, df = Inf)
  lwr <- pred_mean - q * pred_sd
  upr <- pred_mean + q * pred_sd
  preds <- data.frame(
    pred_mean = pred_mean,
    pred_sd = pred_sd,
    lwr = lwr,
    upr = upr
  )
  return(preds)
}

#' compute_SE
#'
#' @description Compute the Squared Error (SE) score
#' @param actual A numeric vector of actual weights
#' @param predicted_mean A numeric vector of predicted means
#' @return The SE score
#'
compute_SE <- function(actual, predicted_mean) {
  se <- (actual - predicted_mean)^2
  return(se)
}

#' compute_DS
#'
#' @description Compute the Dawid-Sebastiani (DS) score
#' @param actual A numeric vector of actual weights
#' @param predicted_mean A numeric vector of predicted means
#' @param sd A numeric vector of standard deviations
#' @return The DS score
#'
#'
compute_DS <- function(actual, predicted_mean, sd) {
  ds <- ((actual - predicted_mean)^2) / sd^2 + 2 * log(sd)
  return(ds)
}


#' leave1out
#' 
#' @description Performs leave-one-out cross-validation for the selected model for each observation
#' @param data A data.frame containing the required predictors
#' @param model Must be A or B
#' @return A data.frame that extends the original data frame with four additional columns mean, sd, se and ds of leave-one-out prediction means
#' 
leave1out <- function(data, model=c("A", "B")) {
  data <- data %>% mutate(mean = NA_real_, sd = NA_real_)
  for (i in seq_len(nrow(data))) {
    est <- filament1_estimate(data[-i,], model)
    pred <- filament1_predict(est, newdata = data[i,], model)
    pred_mean <- pred$pred_mean
    pred_sd <- pred$pred_sd
    data[i, "mean"] <- pred_mean
    data[i, "sd"] <- pred_sd
  }
  data <- data %>%
    mutate(
      se = compute_SE(Actual_Weight, mean),
      ds = compute_DS(Actual_Weight, mean, sd)
    )
  data
}

#' arch_loglike
#' 
#' @description If a data.frame with columns N and phi is provided, the log-likelihood for each row-pair (N, φ) should be returned.
#' @param data a data.frame with columns N and phi
#' @param y observations
#' @return log-likelihood for each row-pair (N, φ)
#' 
arch_loglike <- function(data, y) {
  if (is.data.frame(data) && "N" %in% names(data) && "phi" %in% names(data)) {
    log_like <- rep(-Inf, nrow(data))
    for (i in seq_len(nrow(data))) {
      N <- data$N[i]
      phi <- data$phi[i]
      y1 <- y[1]
      y2 <- y[2]
      log_like[i] <- -lgamma(y1 + 1) - lgamma(y2 + 1) - lgamma(N - y1 + 1) - lgamma(N - y2 + 1) +
          2 * lgamma(N + 1) + (y1 + y2) * log(phi) + (2 * N - y1 - y2) * log(1 - phi)
    }
    return(log_like)
  }
}

#' estimate
#' 
#' @description Implements Monte Carlo integration method to approximate py(y), E(N|y), and E(phi|y).
#' @param y numeric vector of observations
#' @param xi probability parameter for geom dist
#' @param a shape1 parameter for beta dist
#' @param b shape2 parameter for beta dist
#' @param K number of Monte Carlo samples
#' @return data.frame of p_y, E_N, and E_phi
#' 
estimate <- function(y, xi, a, b, K) {
  N_samples <- rgeom(K, prob = xi)
  phi_samples <- rbeta(K, shape1 = a, shape2 = b)
  
  log_like <- arch_loglike(data.frame(N = N_samples, phi = phi_samples), y = y)
  
  p_y <- mean(exp(log_like))
  E_N <- mean(N_samples * exp(log_like) / p_y)
  E_phi <- mean(phi_samples * exp(log_like) / p_y)
  
  return(data.frame(p_y = p_y, E_N = E_N, E_phi = E_phi))
}

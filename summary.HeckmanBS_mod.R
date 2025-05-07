summary.HeckmanBS_mod <- function(object, ...) {
  est <- object$coefficients
  se  <- object$prop_sigmaBS

  # Trata possíveis NA em erros padrão
  tval <- est / se
  pval <- 2 * pnorm(-abs(tval))
  pval[!is.finite(pval)] <- NA
  tval[!is.finite(tval)] <- NA

  summary_table <- data.frame(
    Estimate     = est,
    `Std. Error` = se,
    `t value`    = tval,
    `Pr(>|t|)`   = pval
  )

  # Adiciona nomes (se existirem em start)
  if (!is.null(names(est)) && length(unique(names(est))) == length(est)) {
    rownames(summary_table) <- names(est)
  } else {
    rownames(summary_table) <- paste0("par", seq_along(est))
  }
  

  # Separação por blocos
  nXS <- object$NXS
  nXO <- object$NXO
  idx_gamma <- 1:nXS
  idx_beta  <- (nXS + 1):(nXS + nXO)
  idx_phi_rho <- (nXS + nXO + 1):length(est)

  coef_sel <- summary_table[idx_gamma, , drop = FALSE]
  coef_out <- summary_table[idx_beta,  , drop = FALSE]
  coef_err <- summary_table[idx_phi_rho, , drop = FALSE]

  cat("--------------------------------------------------------------\n")
  cat(" Birnbaum-Saunders Heckman Model (Modificado)\n")
  cat("--------------------------------------------------------------\n")
  cat("Log-Likelihood:", round(object$loglik, 3), "\n")
  cat("AIC:", round(object$aic, 0), "BIC:", round(object$bic, 0), "\n")
  cat("Number of observations: (", object$N0, "censored and", object$N1, "observed )\n")
  cat(length(est), "parameters ( df =", object$df, ")\n")
  cat("--------------------------------------------------------------\n")
  cat("Probit selection equation:\n")
  print(coef_sel)
  cat("--------------------------------------------------------------\n")
  cat("Outcome equation:\n")
  print(coef_out)
  cat("--------------------------------------------------------------\n")
  cat("Error terms:\n")
  print(coef_err)
  cat("--------------------------------------------------------------\n")

  invisible(summary_table)
}

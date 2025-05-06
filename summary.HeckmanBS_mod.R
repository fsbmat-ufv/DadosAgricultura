summary.HeckmanBS_mod <- function(object, ...) {
  est <- object$coefficients
  se  <- object$prop_sigmaBS
  
  # Trata possíveis NA em erros padrão
  tval <- est / se
  pval <- 2 * pnorm(-abs(tval))
  pval[!is.finite(pval)] <- NA
  tval[!is.finite(tval)] <- NA
  
  summary_table <- data.frame(
    Estimate   = est,
    `Std. Error` = se,
    `t value`    = tval,
    `Pr(>|t|)`   = pval
  )
  
  # Separação por blocos
  nXS <- object$NXS
  nXO <- object$NXO
  coef_sel <- summary_table[1:nXS, , drop = FALSE]
  coef_out <- summary_table[(nXS + 1):(nXS + nXO), , drop = FALSE]
  coef_err <- summary_table[(nXS + nXO + 1):nrow(summary_table), , drop = FALSE]
  
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

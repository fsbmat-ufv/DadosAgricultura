summary.HeckmanBS_mod <- function(object, ...) {
  fisher_infoBS <- object$fisher_info
  prop_sigmaBS <- object$se
  coeffsBS <- object$coefficients
  counts <- object$counts
  value <- -object$loglik
  loglik <- object$loglik
  NObs <- object$nObs
  nParam <- length(coeffsBS)
  df <- object$df
  NXS <- object$NXS
  NXO <- object$NXO
  N0 <- object$N0
  N1 <- object$N1
  aic <- object$aic
  bic <- object$bic
  
  # Ajuste seguro para coefTable
  # tb <- miscTools::coefTable(coeffsBS, prop_sigmaBS, df = df)  # não será usado diretamente
  
  tb1 <- tryCatch(
    miscTools::coefTable(coeffsBS[1:NXS], prop_sigmaBS[1:NXS], df = df),
    error = function(e) matrix(NA, nrow = NXS, ncol = 4,
                               dimnames = list(names(coeffsBS[1:NXS]),
                                               c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
  )
  tb2 <- tryCatch(
    miscTools::coefTable(coeffsBS[(NXS + 1):(NXS + NXO)], prop_sigmaBS[(NXS + 1):(NXS + NXO)], df = df),
    error = function(e) matrix(NA, nrow = NXO, ncol = 4,
                               dimnames = list(names(coeffsBS[(NXS + 1):(NXS + NXO)]),
                                               c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
  )
  tb3 <- tryCatch(
    miscTools::coefTable(coeffsBS[(NXS + NXO + 1):(NXS + NXO + 2)],
                         prop_sigmaBS[(NXS + NXO + 1):(NXS + NXO + 2)], df = df),
    error = function(e) matrix(NA, nrow = 2, ncol = 4,
                               dimnames = list(names(coeffsBS[(NXS + NXO + 1):(NXS + NXO + 2)]),
                                               c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
    )

cat("\n")
cat("--------------------------------------------------------------\n")
cat("   Birnbaum-Saunders Heckman Model (Package: ssmodels)        \n")
cat("--------------------------------------------------------------\n")
cat("--------------------------------------------------------------\n")
cat("Maximum Likelihood estimation \n")
cat("optim function with method BFGS-iterations numbers:", counts, "\n")
cat("Log-Likelihood:", value, "\n")
cat("AIC:", aic, "BIC:", bic, "\n")
cat("Number of observations:", NObs, "(", N0, "censored and", N1, "observed )\n")
cat(nParam, "free parameters", "(", "df=", df, ")\n")
cat("--------------------------------------------------------------\n")
cat("Probit selection equation:\n")
printCoefmat(tb1, signif.stars = TRUE, signif.legend = FALSE, digits = 4)
cat("--------------------------------------------------------------\n")
cat("Outcome equation:\n")
printCoefmat(tb2, signif.stars = TRUE, signif.legend = FALSE, digits = 4)
cat("--------------------------------------------------------------\n")
cat("Error terms:\n")
printCoefmat(tb3, signif.stars = TRUE, signif.legend = TRUE, digits = 4)
cat("--------------------------------------------------------------\n")
}


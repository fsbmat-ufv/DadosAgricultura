HeckmanBS_mod <- function (selection, outcome, data = sys.frame(sys.parent()), 
                           start = NULL) 
{
  # Verificações de tipo
  if (!inherits(selection, "formula")) stop("'selection' deve ser uma fórmula.")
  if (!inherits(outcome, "formula")) stop("'outcome' deve ser uma fórmula.")
  
  # ConmfS# Construção robusta do model.frame da equação de seleção
  mfS <- model.frame(
    formula = selection,
    data = data,
    drop.unused.levels = TRUE,
    na.action = na.pass
  )
  mtS <- terms(mfS)
  XS <- model.matrix(mtS, mfS)
  NXS <- ncol(XS)
  YS <- model.response(mfS)
  YSLevels <- levels(as.factor(YS))
  
  # Construção robusta do model.frame da equação de resultado
  mfO <- model.frame(
    formula = outcome,
    data = data,
    drop.unused.levels = TRUE,
    na.action = na.pass
  )
  mtO <- terms(mfO)
  XO <- model.matrix(mtO, mfO)
  NXO <- ncol(XO)
  YO <- model.response(mfO)
  
  # Valores iniciais, se não fornecidos
  #if (is.null(start)) 
  #  start <- step2(YS, XS, log(YO), XO)
  if (is.null(start)) {
    message("Start não fornecido — utilizando valores iniciais padrão.")
    start <- c(rep(0, ncol(XS) + ncol(XO)), 1, 0)
  }
  
  loglik_BS <- function(par) {
    # Proteção contra parâmetros inválidos
    if (any(!is.finite(par))) return(NA)
    
    n <- length(YO)
    NXS <- dim(model.matrix(~XS))[2] - 1
    NXO <- dim(model.matrix(~XO))[2] - 1
    igamma <- 1:NXS
    ibeta <- seq(tail(igamma, 1) + 1, length = NXO)
    iphi1 <- tail(ibeta, 1) + 1
    irho <- tail(iphi1, 1) + 1
    gamma <- par[igamma]
    beta <- par[ibeta]
    phi1 <- par[iphi1]
    if (!is.finite(phi1) || phi1 < 1e-3 || phi1 > 100) return(NA)
    rho_star <- par[irho]
    rho <- 2 / (1 + exp(-rho_star)) - 1 # mapeia ℝ → (−1, 1)
    phi2 <- 1
    XS0 <- XS[YS == 0, , drop = FALSE]
    XS1 <- XS[YS == 1, , drop = FALSE]
    YO[is.na(YO)] <- 0
    YO1 <- YO[YS == 1]
    XO1 <- XO[YS == 1, , drop = FALSE]
    N0 <- sum(YS == 0)
    N1 <- sum(YS == 1)
    XS0.g <- exp(as.numeric((XS0) %*% gamma))
    XS1.g <- exp(as.numeric((XS1) %*% gamma))
    XO1.b <- exp(as.numeric((XO1) %*% beta))
    term0 <- ((YO1 * (phi1 + 1)/(phi1 * XO1.b))^(1/2) - ((phi1 * 
                                                            XO1.b)/(YO1 * (phi1 + 1)))^(1/2))
    term1 <- exp((-phi1/4) * (term0)^2)#exp(t1)
    term2 <- (((phi1 + 1)/(phi1 * XO1.b * YO1))^(1/2) + ((phi1 * 
                                                            XO1.b)/((phi1 + 1) * (YO1^3)))^(1/2))
    term3 <- (1/(2 * sqrt(2 * pi))) * ((phi1/2)^(1/2))
    term4 <- ((phi2 + 1)/(2 * XS1.g * (1 - rho^2)))^(1/2)
    term5 <- ((phi2 * XS1.g)/(phi2 + 1)) - 1
    term6 <- rho * (phi1/(2 * (1 - rho^2)))^(1/2)
    integrand <- term4 * term5 + term6 * term0
    term7 <- pnorm(integrand, log.p = TRUE)
    term8 <- ((phi2/2)^(1/2)) * (((phi2 + 1)/(phi2 * XS0.g))^(1/2) - 
                                   ((phi2 * XS0.g)/(phi2 + 1))^(1/2))
    FT2 <- pnorm(term8, log.p = TRUE)
    term1[term1 < 1e-300] <- 1e-300
    logterm1 <- log(term1)
    logterm2 <- log(term2)
    
    if (any(!is.finite(logterm1)) || any(!is.finite(logterm2))) return(NA)
    
    ll <- sum(logterm1 + logterm2 + log(term3) + term7) + sum(FT2)

    return(sum(ll))
  }
  
  gradlik_BS <- function(par) {
    library(matrixStats)
    
    n <- length(YO)
    NXS <- dim(model.matrix(~XS))[2] - 1
    NXO <- dim(model.matrix(~XO))[2] - 1
    igamma <- 1:NXS
    ibeta <- seq(tail(igamma, 1) + 1, length = NXO)
    iphi1 <- tail(ibeta, 1) + 1
    irho <- tail(iphi1, 1) + 1
    
    gamma <- par[igamma]
    beta <- par[ibeta]
    phi1 <- par[iphi1]
    rho_star <- par[irho]
    rho <- 2 / (1 + exp(-rho_star)) - 1
    phi2 <- 1
    
    XS0 <- XS[YS == 0, , drop = FALSE]
    XS1 <- XS[YS == 1, , drop = FALSE]
    YO[is.na(YO)] <- 0
    YO1 <- YO[YS == 1]
    XO1 <- XO[YS == 1, , drop = FALSE]
    
    XS0.g <- exp(as.numeric(XS0 %*% gamma))
    XS1.g <- exp(as.numeric(XS1 %*% gamma))
    XO1.b <- exp(as.numeric(XO1 %*% beta))
    
    mu1 <- exp(as.numeric((XO) %*% beta))
    mu2 <- exp(as.numeric((XS) %*% gamma))
    
    term0 <- sqrt(YO1 * (phi1 + 1) / (phi1 * XO1.b)) - 
      sqrt((phi1 * XO1.b) / (YO1 * (phi1 + 1)))
    term1 <- exp((-phi1 / 4) * (term0)^2)
    term2 <- sqrt((phi1 + 1) / (phi1 * XO1.b * YO1)) + 
      sqrt((phi1 * XO1.b) / ((phi1 + 1) * YO1^3))
    term3 <- (1 / (2 * sqrt(2 * pi))) * sqrt(phi1 / 2)
    
    safe_denom4 <- pmax(2 * XS1.g * (1 - rho^2), 1e-12)
    term4 <- sqrt((phi2 + 1) / safe_denom4)
    term5 <- (phi2 * XS1.g) / (phi2 + 1) - 1
    term6 <- rho * sqrt(phi1 / (2 * pmax(1 - rho^2, 1e-12)))
    
    integrand <- term4 * term5 + term6 * term0
    logPhi_int <- pmax(pnorm(integrand, log.p = TRUE), -700)
    lambda_I <- exp(dnorm(integrand, log = TRUE) - logPhi_int)
    
    term8 <- sqrt(phi2 / 2) * 
      (sqrt((phi2 + 1) / (phi2 * XS0.g)) - 
         sqrt((phi2 * XS0.g) / (phi2 + 1)))
    logPhi_t8 <- pmax(pnorm(term8, log.p = TRUE), -700)
    lambda_T8 <- exp(dnorm(term8, log = TRUE) - logPhi_t8)
    
    term9  <- (-1/2) * (sqrt(YO1 * (phi1 + 1) / (phi1 * XO1.b)) + 
                          sqrt((phi1 * XO1.b) / (YO1 * (phi1 + 1)))) * XO1
    term10 <- (1/2) * (sqrt((XO1.b * phi1) / (YO1^3 * (phi1 + 1))) - 
                         sqrt((phi1 + 1) / (YO1 * phi1 * XO1.b))) * XO1
    term11 <- (-1/2) * sqrt(phi2 / 2) * 
      (sqrt((phi2 + 1) / (phi2 * XS0.g)) + 
         sqrt((phi2 * XS0.g) / (phi2 + 1))) * XS0
    
    term12 <- (-1/2) * (sqrt(YO1 / (phi1^3 * XO1.b * (phi1 + 1))) + 
                          sqrt(XO1.b / (YO1 * (phi1 + 1)^3 * phi1)))
    term13 <- (1/2) * (sqrt(XO1.b / (phi1 * (phi1 + 1)^3 * YO1^3)) - 
                         sqrt(1 / (phi1^3 * (phi1 + 1) * XO1.b * YO1)))
    
    safe_denom_df4g <- sqrt(pmax(XS1.g * (1 - rho^2), 1e-12))
    df4g <- lambda_I * (1 / 4) * ((XS1.g + 2) / safe_denom_df4g) * XS1
    
    df5g <- lambda_T8 * term11
    df1b <- (-phi1 / 2) * term0 * term9
    df2b <- term10 / term2
    df4b <- lambda_I * term9 * term6
    
    df1phi1 <- - (term0^2) / 4 - (phi1 / 2) * (term0 * term12)
    df2phi1 <- term13 / term2
    df3phi1 <- 1 / (8 * term3 * sqrt(pi * phi1))
    
    safe_denom14 <- pmax(2 * phi1 * (1 - rho^2), 1e-8)
    term14 <- rho / (2 * sqrt(safe_denom14))
    df4phi1 <- lambda_I * (term0 * term14 + term6 * term12)
    
    safe_rho <- pmax(abs(rho), 1e-6)
    safe_rho2 <- pmax(1 - rho^2, 1e-6)
    
    if (abs(rho) < 1e-8) {
      df4rho_raw <- lambda_I * term0 * sqrt(phi1 / 2)
    } else {
      df4rho_raw <- lambda_I * (term5 * term4 * (rho / safe_rho2) + 
                                  term0 * term6 / (safe_rho * safe_rho2))
    }
    
    df4rho_raw <- pmin(pmax(df4rho_raw, -1e5), 1e5)
    d_rho <- (4 * exp(-rho_star)) / (1 + exp(-rho_star))^2
    gr_rho_star <- df4rho_raw * d_rho
    
    nParam <- length(par)
    gradient <- matrix(0, nrow = n, ncol = nParam)
    gradient[YS == 0, igamma] <- df5g
    gradient[YS == 1, igamma] <- df4g
    gradient[YS == 1, ibeta]  <- df1b + df2b + df4b
    gradient[YS == 1, iphi1]  <- df1phi1 + df2phi1 + df3phi1 + df4phi1
    gradient[YS == 1, irho]   <- gr_rho_star
    
    grad_final <- colSums(gradient)
    if (any(!is.finite(grad_final))) return(rep(NA, length(par)))
    return(grad_final)
  }
  
  
  
  
  
  
  
  theta_BS <- optim(start, loglik_BS, gradlik_BS, method = "BFGS", hessian = TRUE,
                    control = list(fnscale = -1, maxit = 1000))
  names(theta_BS$par) <- c(colnames(XS), colnames(XO), "sigma", 
                           "rho")
  
  a <- start
  a1 <- theta_BS$par
  a2 <- theta_BS$value
  a3 <- theta_BS$counts[2]
  a4 <- theta_BS$hessian
  
  # Tentativa de inversão da Hessiana
  a5 <- tryCatch(solve(-a4), error = function(e) matrix(NA, nrow = length(a1), ncol = length(a1)))
  a6 <- suppressWarnings(sqrt(diag(a5)))
  a6[!is.finite(a6)] <- NA
  
  # Ajuste de rho e seu erro padrão via delta method
  rho_star <- a1[irho]
  rho_hat <- 2 / (1 + exp(-rho_star)) - 1
  g_prime <- 4 * exp(-rho_star) / (1 + exp(-rho_star))^2
  se_rho_star <- a6[irho]
  se_rho <- abs(g_prime) * se_rho_star
  
  # Substituindo o valor de rho e seu erro padrão na saída
  a1[irho] <- rho_hat
  a6[irho] <- se_rho
  
  # Informações adicionais
  a7 <- levels(as.factor(YS))
  a8 <- length(YS)
  a9 <- length(start)
  a10 <- sum(YS == 0)
  a11 <- sum(YS == 1)
  a12 <- ncol(XS)
  a13 <- ncol(XO)
  a14 <- a8 - a9
  a15 <- -2 * a2 + 2 * a9
  a16 <- -2 * a2 + a9 * log(a8)
  
  # Resultado final
  result <- list(
    coefficients     = a1,
    value            = a2,
    loglik           = -a2,
    counts           = a3,
    hessian          = a4,
    fisher_infoBS    = a5,
    prop_sigmaBS     = a6,
    level            = a7,
    nObs             = a8,
    nParam           = a9,
    N0               = a10,
    N1               = a11,
    NXS              = a12,
    NXO              = a13,
    df               = a14,
    aic              = a15,
    bic              = a16,
    initial.value    = start
  )
  class(result) <- "HeckmanBS_mod"
  result
  
}

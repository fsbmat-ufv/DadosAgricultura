library(numDeriv)

# Vetor de parâmetros estimados
theta <- theta_BS$par

# Log-verossimilhança no ponto original
loglik_ana <- loglik_BS(theta)

# Perturbação pequena
eps <- 1e-6

# Inicializa vetor para armazenar os valores perturbados
loglik_eps <- numeric(length(theta))
approx_deriv <- numeric(length(theta))

# Calcula a derivada numérica para cada parâmetro
for (j in 1:length(theta)) {
  delta <- rep(0, length(theta))
  delta[j] <- eps
  loglik_eps[j] <- loglik_BS(theta + delta)
  approx_deriv[j] <- (loglik_eps[j] - loglik_ana) / eps
}

# Resultado em data.frame
result <- data.frame(
  ParamIndex         = 1:length(theta),
  LogLik_Analitica   = rep(loglik_ana, length(theta)),
  LogLik_Perturbada  = loglik_eps,
  Derivada_Numérica  = approx_deriv
)

# Exibe
print(result, digits = 6)

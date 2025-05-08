library(numDeriv)

# Avaliação numérica do gradiente
grad_num <- grad(func = loglik_BS, x = theta_BS$par)

# Avaliação analítica do gradiente
grad_ana <- gradlik_BS(theta_BS$par)

# Comparação lado a lado
cbind(Gradiente_Analitico = round(grad_ana, 6),
      Gradiente_Numerico  = round(grad_num, 6),
      Diferenca           = round(grad_ana - grad_num, 6))


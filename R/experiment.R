# # TEST FILE
#
# # parallel montecarlo
# nsim <- 20 #1000
# n <- 2000
# TT <- 5
# K <- 2
# d <- 1
# gamma <- 0.5
# dhat <- 2
# latent <- matrix(cbind(1, -2), nrow = d, ncol = K)   # d = 1
# #
# # ## Balanced case
# pi <- rep(1/K, K)
# block_size <- round(pi * n)
# blocks <- c()
# for (k in 1:length(block_size)) {
#   blocks <- c(blocks, rep(k, block_size[k]))
# }
# #
# latent_positions <- generate_latent_posSBM(latent, d, block_size)
# data_list <- generate_data(latent_positions, gamma, TT, 100)
# A <- data_list$A
# A0 <- data_list$A0
#
# # M = matrix(c(1:12), ncol=3)
# #
# # incrementByN <- function(x) {
# #   return(x^2)
# # }
# #
# # (y = apply(M, c(1,2),incrementByN))
#
# #
# gamma_hat <- CMLE_est(A)
# P_hat <- estimate_GRDPG(A, d, sim = 100)
# # Omega_init <- latent_positions %*% t(latent_positions)
# # xtol = 1e-5
# #
# #
# alpha_hat <- alpha_est(P_hat, A, gamma_hat, d, est_method = "avg")$alpha_hat
# # alpha_hat_stat <- alpha_est_stat(data_list$A0, gamma_hat, d, lb = -10, ub = 10)
# #
#

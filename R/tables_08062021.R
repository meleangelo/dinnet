# # Tables
#
# results_folder <- "/Users/Angelo/Dropbox/network_dynamics/"
#
# # monte carlo K=2, n=2000, T=4
# K <- 2
# n = 2000
# TT = 4
# # latent pos = (1,-2)
# # gamma = (0.1,0.2,0.3,0.4,0.5)
# gamma <- seq(0.1,0.5, by=0.1)
# # nsim = 100
# nsim = 100
#
#
# # table with summary results
# table_results <- data.frame(matrix(NA, nrow = 5, ncol = 7))
# table_results[,1] <- gamma
#
# for (j in 1:5){
#   namefile <- paste(results_folder, "mc_nsim",nsim,"_n",n,"_T", TT,"_K",K, "_gamma", gamma[j],".Rdata", sep = "")
#   load(namefile)
#   rmse_gamma <- sqrt( mean( (table_gamma[,1] - gamma[j])^2 ) )
#   se_gamma <- sd(table_gamma[,1])
#   table_results[j,2] <- rmse_gamma
#   table_results[j,3] <- se_gamma
#
#   rmse_alpha1 <- sqrt( mean( (-table_alpha_g[,1] - 1)^2 ) )
#   se_alpha1 <- sd(table_alpha_g[,1])
#   table_results[j,4] <- rmse_alpha1
#   table_results[j,5] <- se_alpha1
#
#   rmse_alpha2 <- sqrt( mean( (-table_alpha_g[,2] - (-2))^2 ) )
#   se_alpha2 <- sd(table_alpha_g[,2])
#   table_results[j,6] <- rmse_alpha2
#   table_results[j,7] <- se_alpha2
#
# }
#
# names(table_results) <- c("true gamma", "RMSE gamma", "se gamma",
#                           "RMSE alpha_1", "se alpha_1",
#                           "RMSE alpha_2", "se alpha_2")
#
# library(xtable)
# xtable(table_results, digits = 4)
#
#
#
#
#
#
# # monte carlo K=4, n=2000, T=4
# K = 4
# n = 2000
# TT = 4
# # latent pos = (1, 0.5, -0.7, -1.2)
#
# # gamma = (0.1,0.2,0.3,0.4,0.5)
# gammavals <- seq(0.1,0.5, by=0.1)
# # nsim = 100
# nsim = 100
#
# # Table with summary results
# table_results <- data.frame(matrix(NA, nrow = 5, ncol = 11))
#
#
# table_results[,1] <- gammavals
#
# for (j in 1:5){
#   #gammavals <- seq(0.1,0.5,by = 0.1)
#   results_folder <- "/Users/Angelo/Dropbox/network_dynamics/"
#   namefile <- paste(results_folder, "mc_nsim",nsim,"_n",n,"_T", TT,"_K",K, "_gamma", gammavals[j],".Rdata", sep = "")
#   load(namefile)
#   rmse_gamma <- sqrt( mean( (table_gamma[,1] - gammavals[j])^2 ) )
#   se_gamma <- sd(table_gamma[,1])
#   table_results[j,2] <- rmse_gamma
#   table_results[j,3] <- se_gamma
#
#   rmse_alpha1 <- sqrt( mean( (-table_alpha_g[,1] - 1)^2 ) )
#   se_alpha1 <- sd(table_alpha_g[,1])
#   table_results[j,4] <- rmse_alpha1
#   table_results[j,5] <- se_alpha1
#
#   rmse_alpha2 <- sqrt( mean( (-table_alpha_g[,2] - (.5))^2 ) )
#   se_alpha2 <- sd(table_alpha_g[,2])
#   table_results[j,6] <- rmse_alpha2
#   table_results[j,7] <- se_alpha2
#
#   rmse_alpha3 <- sqrt( mean( (-table_alpha_g[,3] - (-0.7))^2 ) )
#   se_alpha3 <- sd(table_alpha_g[,3])
#   table_results[j,8] <- rmse_alpha3
#   table_results[j,9] <- se_alpha3
#
#   rmse_alpha4 <- sqrt( mean( (-table_alpha_g[,4] - (-1.2))^2 ) )
#   se_alpha4 <- sd(table_alpha_g[,4])
#   table_results[j,10] <- rmse_alpha4
#   table_results[j,11] <- se_alpha4
#
# }
#
# names(table_results) <- c("true gamma", "RMSE gamma", "se gamma",
#                           "RMSE alpha_1", "se alpha_1",
#                           "RMSE alpha_2", "se alpha_2",
#                           "RMSE alpha_3", "se alpha_3",
#                           "RMSE alpha_4", "se alpha_4")
#
# library(xtable)
# xtable(table_results, digits = 4)
#
#
#
#
#
#
# # monte carlo K=4, n=5000, T=4
# K = 4
# n = 5000
# TT = 4
# # latent pos = (1, 0.5, -0.7, -1.2)
#
# # gamma = (0.1,0.2,0.3,0.4,0.5)
# gammavals <- seq(0.1,0.5, by=0.1)
# # nsim = 100
# nsim = 100
#
# # Table with summary results
# table_results <- data.frame(matrix(NA, nrow = 5, ncol = 11))
#
#
# table_results[,1] <- gammavals
#
# for (j in 1:5){
#   #gammavals <- seq(0.1,0.5,by = 0.1)
#   results_folder <- "/Users/Angelo/Dropbox/network_dynamics/"
#   namefile <- paste(results_folder, "mc_nsim",nsim,"_n",n,"_T", TT,"_K",K, "_gamma", gammavals[j],".Rdata", sep = "")
#   load(namefile)
#   rmse_gamma <- sqrt( mean( (table_gamma[,1] - gammavals[j])^2 ) )
#   se_gamma <- sd(table_gamma[,1])
#   table_results[j,2] <- rmse_gamma
#   table_results[j,3] <- se_gamma
#
#   rmse_alpha1 <- sqrt( mean( (-table_alpha_g[,1] - 1)^2 ) )
#   se_alpha1 <- sd(table_alpha_g[,1])
#   table_results[j,4] <- rmse_alpha1
#   table_results[j,5] <- se_alpha1
#
#   rmse_alpha2 <- sqrt( mean( (-table_alpha_g[,2] - (.5))^2 ) )
#   se_alpha2 <- sd(table_alpha_g[,2])
#   table_results[j,6] <- rmse_alpha2
#   table_results[j,7] <- se_alpha2
#
#   rmse_alpha3 <- sqrt( mean( (-table_alpha_g[,3] - (-0.7))^2 ) )
#   se_alpha3 <- sd(table_alpha_g[,3])
#   table_results[j,8] <- rmse_alpha3
#   table_results[j,9] <- se_alpha3
#
#   rmse_alpha4 <- sqrt( mean( (-table_alpha_g[,4] - (-1.2))^2 ) )
#   se_alpha4 <- sd(table_alpha_g[,4])
#   table_results[j,10] <- rmse_alpha4
#   table_results[j,11] <- se_alpha4
#
# }
#
# names(table_results) <- c("true gamma", "RMSE gamma", "se gamma",
#                           "RMSE alpha_1", "se alpha_1",
#                           "RMSE alpha_2", "se alpha_2",
#                           "RMSE alpha_3", "se alpha_3",
#                           "RMSE alpha_4", "se alpha_4")
#
# library(xtable)
# xtable(table_results, digits = 4)

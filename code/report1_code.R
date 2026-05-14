library(xtable)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(RColorBrewer)

mse_var_estimators <- function(n, theta, nsim = 5000) {
  
  true_var <- theta * (1 - theta)
  
  sims_mse <- replicate(nsim, {
    x <- rbinom(n, size = 1, prob = theta)
    xbar <- mean(x)
    s2_0 <- sum((x - xbar)^2) / n
    s2   <- sum((x - xbar)^2) / (n - 1) 
    c(s2_0 = s2_0, s2 = s2)
  })
  
  est_means <- rowMeans(sims_mse)
  biases <- est_means - true_var
  variances <- apply(sims_mse, 1, var)
  mses <- variances + biases^2
  
  list(
    n = n,
    theta = theta,
    true_var = true_var,
    mean = est_means,
    bias = biases,
    var = variances,
    mse = mses
  )
}

params <- expand.grid(
  n = c(10, 20, 50),
  theta = seq(0.01, 0.99, by = 0.01) 
)

results <- apply(params, 1, function(row) {
  n <- as.numeric(row[1])
  theta <- as.numeric(row[2])
  res <- mse_var_estimators(n, theta, nsim = 5000) 
  data.frame(
    n = n,
    theta = theta,
    true_var = res$true_var,
    S2_0_mse = res$mse["s2_0"],
    S2_mse = res$mse["s2"]
  )
})

results <- do.call(rbind, results)


results_long <- pivot_longer(
  results,
  cols = c(S2_0_mse, S2_mse), 
  names_to = "Estimator",     
  values_to = "MSE"          )

  
p1 <- ggplot(results_long %>% filter(n == 10), aes(x = theta, y = MSE, color = Estimator)) +
  geom_line(linewidth = 1) +
  labs(title = "MSE comparison for n = 10",
       x = "Probability (theta)", 
       y = "Mean Squared Error (MSE)",
       color = "Estimator") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")


p2 <- ggplot(results_long %>% filter(n == 20), aes(x = theta, y = MSE, color = Estimator)) +
  geom_line(linewidth = 1) +
  labs(title = "MSE comparison for n = 20",
       x = "Probability (theta)", 
       y = "Mean Squared Error (MSE)",
       color = "Estimator") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")


p3 <- ggplot(results_long %>% filter(n == 50), aes(x = theta, y = MSE, color = Estimator)) +
  geom_line(linewidth = 1) +
  labs(title = "MSE comparison for n = 50",
       x = "Probability (theta)", 
       y = "Mean Squared Error (MSE)",
       color = "Estimator") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p1)
print(p2)
print(p3)



risk_entropy_da_b <- function(n, theta, a, b) {
  
  x <- 0:n
  d_x <- (x + a) / (n + a + b)
  
  L <- theta * log(theta / d_x) + (1 - theta) * log((1 - theta) / (1 - d_x))
  
  p_x <- dbinom(x, size = n, prob = theta)
  
  R <- sum(p_x * L)
  return(R)
}

risk_curve <- function(n, thetas, a, b) {
  sapply(thetas, function(th) risk_entropy_da_b(n = n, theta = th, a = a, b = b))
}

n_values <- c(5, 20, 100) 
thetas <- seq(0.01, 0.99, by = 0.01)
param_sets <- list(
  list(a = 1, b = 1, label = "a=1, b=1"),
  list(a = NA, b = NA, label = "a=n/2, b=n/2 (shrink to 1/2)")
)

df_list <- list()
for (n in n_values) {
  for (ps in param_sets) {
    if (is.na(ps$a)) {
      a <- n / 2
      b <- n / 2
      label <- paste0("a=", n/2, ", b=", n/2) 
    } else {
      a <- ps$a; b <- ps$b; label <- ps$label
    }
    Rvals <- risk_curve(n = n, thetas = thetas, a = a, b = b)
    df_list[[length(df_list) + 1]] <- data.frame(
      n = n,
      theta = thetas,
      risk = Rvals,
      estimator = label
    )
  }
}
df <- bind_rows(df_list)


plot_list <- list()
for (n_val in n_values) {
  df_subset <- filter(df, n == n_val)
  
  p <- ggplot(df_subset, aes(x = theta, y = risk, color = estimator)) +
    geom_line(linewidth = 1) +
    labs(
      title = paste("Risk fuction R(theta, d_{a,b}) for n =", n_val),  
      subtitle = "Estimators comparison d_{a,b}(X) = (X+a)/(n+a+b)",
      x = expression(theta),
      y = "R(theta, d)",
      color = "Estimator"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
  
  plot_list[[paste0("n_", n_val)]] <- p
}


for (p in plot_list) {
  print(p)
}




r_pi <- function(n, a, b, alpha = 1, beta = 1, rel.tol = 1e-6) {
  integrand <- function(theta) {
    sapply(theta, function(th) risk_entropy_da_b(n = n, theta = th, a = a, b = b)) *
      dbeta(theta, alpha, beta)
  }
  
  res <- try(integrate(integrand, lower = 0, upper = 1, rel.tol = rel.tol, 
                       subdivisions = 200L), silent = TRUE)
  res$value
}

grid_search_ab <- function(n, alpha = 1, beta = 1, 
                           a_seq = seq(0.01, n, length.out = 40), 
                           b_seq = seq(0.01, n, length.out = 40)) {
  best_val <- Inf
  best_ab <- c(NA, NA)
  for (a in a_seq) {
    for (b in b_seq) {
      val <- r_pi(n, a, b, alpha, beta)
      if (val < best_val) {
        best_val <- val
        best_ab <- c(a, b)
      }
    }
  }
  list(a = best_ab[1], b = best_ab[2], risk = best_val)
}

refine_opt_ab <- function(n, a0, b0, alpha = 1, beta = 1) {
  fun_log <- function(par) {
    a <- exp(par[1]); b <- exp(par[2])
    r_pi(n, a, b, alpha, beta)
  }
  start <- log(c(a0, b0))
  opt <- optim(start, fn = fun_log, method = "Nelder-Mead",
               control = list(maxit = 500))
  a_opt <- exp(opt$par[1]); b_opt <- exp(opt$par[2])
  list(a = a_opt, b = b_opt, risk = opt$value, optim_result = opt)
}

find_optimal_ab <- function(n, alpha = 1, beta = 1,
                            a_seq = seq(0.01, n, length.out = 40),
                            b_seq = seq(0.01, n, length.out = 40)) {
  grid <- grid_search_ab(n, alpha, beta, a_seq, b_seq)
  refined <- refine_opt_ab(n, grid$a, grid$b, alpha, beta)
  list(grid = grid, refined = refined)
}

n_values <- c(5, 20, 50)

results <- lapply(n_values, function(n) {
  res <- find_optimal_ab(n, alpha = 1, beta = 1,
                         a_seq = seq(0.01, n, length.out = 40),
                         b_seq = seq(0.01, n, length.out = 40))
  risk_1_1 <- r_pi(n, 1, 1, 1, 1)
  list(n = n, grid = res$grid, refined = res$refined, risk_a1b1 = risk_1_1)
})

  



results_df <- do.call(rbind, lapply(results, function(res) {
  data.frame(
    n = res$n,
    a_opt = res$refined$a,
    b_opt = res$refined$b,
    risk_opt = res$refined$risk,
    risk_a1b1 = res$risk_a1b1
  )
}))

print(results_df)

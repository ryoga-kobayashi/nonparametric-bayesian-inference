rm(list = ls(all = TRUE))

# install.packages("mvtnorm")
library(mvtnorm)
# install.packages("MCMCpack")
library(MCMCpack)

data(iris)
x = as.matrix(iris[, 1:4])


delta = function(z, k){
    ifelse(test = z == k, yes = 1, no = 0)
}

rmulti = function(n = 1, p){
    K = length(p)
    random_variables = rmultinom(n = n, size = 1, prob = p)
    out = rep(0, n)
    for(i in 1:n){
        for(k in 1:K){
            if(random_variables[k, i] == 1) out[i] = k
        }
    }
    return(out)
}

# use the gibbs_sampling_p function
calculate_n_k = function(z, K){
    n_k = rep(0, K)
    for(k in 1:K){
        n_k[k] = sum(delta(z, k))
    }
    return(n_k)
}

calculate_x_bar_k = function(x, z, K){
    d = ncol(x)
    n_k = rep(0, K)
    n_k = calculate_n_k(z = z, K = K)
    x_bar_k = matrix(0, nrow = K, ncol = d)
    for(j in 1:d){
        for(k in 1:K){
            x_bar_k[k, j] = (1/n_k[k])*sum(x[z == k, j])
        }
    }
    return(x_bar_k = x_bar_k)
}

gibbs_sampling_z = function(x, mu, tau, p, K = K){
    n = nrow(x)
    d = ncol(x)
    rv = matrix(0, nrow = n, ncol = K)
    for(i in 1:n){
        for(k in 1:K){
            rv[i, k] = dmvnorm(x[i, ], mean = mu[k, ], sigma = (tau)^{-1}*diag(d))*p[k]
        }
    }
    p_z = matrix(0, nrow = n, ncol = K)
    for(i in 1:n){
        for(k in 1:K){
            p_z[i, k] = rv[i, k]/sum(rv[i, ])
        }
    }
    z = rep(0, n)
    for(i in 1:n){
        z[i] = rmulti(n = 1, p = p_z[i, ])
    }
    return(z)
}

gibbs_sampling_mu = function(x, z, tau, p, K){
    n = nrow(x)
    d = ncol(x)
    n_k = calculate_n_k(z = z, K = K)
    x_bar_k = calculate_x_bar_k(x = x, z = z, K = K)
    mu_k = matrix(0, nrow = K, ncol = d)
    for(k in 1:K){
        mean_k = (n_k[k]/(n_k[k] + rho_0))*x_bar_k[k, ] + (rho_0/(n_k[k] + rho_0))*mu_0
        sigma_k = (tau*(n_k[k] + rho_0))^(-1)*diag(d)
        mu_k[k, ] = rmvnorm(n = 1, mean = mean_k, sigma = sigma_k)
    }
    return(mu_k)
}

gibbs_sampling_tau = function(x, z, mu_0, rho_0, a_0, b_0, K){
    n = nrow(x)
    d = ncol(x)
    a_n = a_0 + (n*d)/2
    out = rep(0, K)
    n_k = calculate_n_k(z = z, K = K)
    x_bar_k = calculate_x_bar_k(x = x, z = z, K = K)
    for(k in 1:K){
        out1 = 0
        for(i in 1:n){
            out1 = out1 + t(x[i, ] - x_bar_k[k, ])%*%(x[i, ] - x_bar_k[k, ])*delta(z = z, k = k)[i]
        }
        out2 = n_k[k]*rho_0/(2*(rho_0 + n_k[k]))*t(x_bar_k[k, ] - mu_0)%*%(x_bar_k[k, ] - mu_0)
        out[k] = (1/2)*out1 + out2
    }
    b_n = b_0 + sum(out)
    tau = rgamma(n = 1, shape = a_n, scale = 1/b_n)
    return(tau)
}

gibbs_sampling_p = function(z, alpha, K){
    n_k = calculate_n_k(z = z, K = K)
    alpha_k = alpha + n_k
    p = rdirichlet(n = 1, alpha = alpha_k)
    return(p)
}

z_map = function(z, K){
    n = nrow(z)
    z_count = matrix(0, nrow = n, ncol = K)
    for(i in 1:n){
        for(k in 1:K){
            z_count[i, k] = length(z[i, z[i, ]==k])
        }
    }
    z_map = rep(0, n)
    for(i in 1:n){
        for(k in 1:K){
            if(z_count[i, k] == max(z_count[i, ])) z_map[i] = k
        }
    }
    return(z_map)
}


gaussian_mixture_model = function(x, K, mu_0 = rep(0, ncol(x)), rho_0 = 0, a_0 = 1, b_0 = 1, p_0 = rep(1/K, K), tau_0 = 1, alpha_0 = rep(1, K)
, iter_max = 10000, burn_in = iter_max*0.1){
    # K = 3
    # mu_0 = rep(0, 4); rho_0 = 0; a_0 = 1; b_0 = 1
    # p_0 = rep(1/K, K); tau_0 = 1
    # alpha_0 = rep(1, 3)
    # iter_max = 100

    n = nrow(x)
    d = ncol(x)

    z = matrix(0, nrow = n, ncol = iter_max)
    mu = array(0, dim = c(K, d, iter_max))
    tau = rep(0, iter_max)
    p = matrix(0, nrow = K, ncol = iter_max)

    x_random = sample(1:n, size = K, replace = FALSE)
    for(k in 1:K){
        mu[k, , 1] = x[x_random[k], ]
        cat(mu[k, , 1], "\n")
    }
    tau[1] = tau_0
    p[, 1] = p_0
    z[, 1] = gibbs_sampling_z(x = x, mu = mu[, , 1], tau = tau[1], p = p[, 1], K = K)

    for(s in 2:iter_max){
        mu[, , s] = gibbs_sampling_mu(x = x, z = z[, s-1], tau = tau[s-1], p = p[, s-1], K = K)
        tau[s] = gibbs_sampling_tau(x = x, z = z[, s-1], mu_0 = mu_0, rho_0 = rho_0, a_0 = a_0, b_0 = b_0, K = K)
        p[, s] = gibbs_sampling_p(z = z[, s-1], alpha = alpha_0, K = K)
        z[, s] = gibbs_sampling_z(x = x, mu = mu[, , s], tau = tau[s], p = p[, s], K = K)
        cat("number of iteration is:", s, "\n")
        cat("class is:", z[, s], "\n")
    }

    z_map = z_map(z = z[, burn_in:iter_max], K = K)

    out = list(z = z[, burn_in:iter_max], z_map = z_map, mu = mu, x = x)
    return(out)
    message("the process has finished")
}

data(iris)
x = as.matrix(iris[, 1:4])
fit_gmm = gaussian_mixture_model(x = x, K = 3)
z_true = iris[, 5]
table(fit_gmm$z_map, z_true)

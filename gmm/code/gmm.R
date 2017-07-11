rm(list = ls(all = TRUE))

# install.packages("mvtnorm")
library(mvtnorm)

rmulti = function(p, n = 1){
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

sample_mean_k = function(x, z, K){
    d = ncol(x)
    n_k = rep(0, K)
    for(k in 1:K){
        n_k[k] = length(z[z == k])
    }
    x_bar_k = matrix(0, nrow = K, ncol = d)
    for(j in 1:d){
        for(k in 1:K){
            x_bar_k[k, j] = (1/n_k[k])*sum(x[z == k, j])
        }
    }
    out = list(n_k = n_k, x_bar_k = x_bar_k)
    return(out)
}

gibbs_sampling_z = function(x, mu, K){
    n = nrow(x)
    d = ncol(x)
    z = rep(0, n)
    for(i in 1:n){
        p = rep(0, length = K)
        for(k in 1:K){
            p[k] = dmvnorm(x[i, ], mean = mu[k, ], sigma = diag(d))
        }
        z[i] = rmulti(p = p/sum(p))
    }
    return(z)
}

gibbs_sampling_mu = function(x, z, K){
    d = ncol(x)

    sample_mean_k = sample_mean_k(x = x, z = z, K = K)
    n_k = sample_mean_k$n_k
    x_bar_k = sample_mean_k$x_bar_k
    mu = matrix(0, nrow = K, ncol = d)
    for(k in 1:K){
        mu[k, ] = rmvnorm(1, mean = (n_k[k]/(n_k[k]+1))*x_bar_k[k, ], sigma = (1/(n_k[k]+1))*diag(d))
    }
    mu[is.nan(mu[, ])] = 0
    return(mu)
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

gaussian_mixture_model = function(x, K, iter_max = 1000, burn_in = 0){
    n = nrow(x)
    d = ncol(x)

    mu = array(0, dim = c(K, d, iter_max))
    z  = array(0, dim = c(n, iter_max))
    cat("number of iteration is:", 1, "\n")

    x_random = sample(1:n, size = K, replace = FALSE)
    for(k in 1:K){
        mu[k, , 1] = x[x_random[k], ]
        cat(mu[k, , 1], "\n")
    }
    z[, 1] = gibbs_sampling_z(x = x, mu = mu[, , 1], K = K)
    cat("class is:", z[, 1], "\n")

    # gibbs sampling
    for(s in 2:iter_max){
        mu[, , s] = gibbs_sampling_mu(x = x, z = z[, s-1], K = K)
        z[, s] = gibbs_sampling_z(x = x, mu = mu[, , s], K = K)
        cat("number of iteration is:", s, "\n")
        cat("class is:", z[, s], "\n")
    }

    z_map = z_map(z = z, K = K)

    out = list(z = z[, burn_in:iter_max], z_map = z_map, mu = mu, x = x)
    return(out)
    message("the process has finished")
}

data(iris)
x = as.matrix(iris[, 1:4])
fit_mgm = gaussian_mixture_model(x = x, K = 3)
z_true = iris[, 5]
table(fit_mgm$z_map, z_true)

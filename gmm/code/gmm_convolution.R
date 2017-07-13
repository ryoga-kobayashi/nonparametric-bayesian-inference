rm(list = ls(all = TRUE))

# install.packages("mvtnorm")
library(mvtnorm)

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

calculate_x_bar = function(x){
    n = nrow(x)
    d = ncol(x)
    x_bar = rep(0, d)
    for(j in 1:d){
        x_bar[j] = (1/n)*sum(x[, j])
    }
    return(x_bar)
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

gaussian_mixture_model = function(x, K, iter_max = 1000, mu_0 = rep(0, ncol(x)), rho_0 = 1, a_0 = 1, b_0 = 1, p_0 = rep(1/K, K), tau_0 = 1, alpha_0 = rep(1, K), alpha = 1, burn_in = iter_max*0.1){

    x = x;
    iter_max = 1000; mu_0 = rep(2, ncol(x)); rho_0 = 1; a_0 = 1; b_0 = 1;
    # p_0 = rep(1/K, K);
    tau_0 = 1;
    # alpha_0 = rep(1, K);
    alpha = 1;
    burn_in = iter_max*0.1
    n = nrow(x)
    d = ncol(x)

    z = matrix(0, nrow = n, ncol = iter_max)
    z[1, ] = 1
    z[, 1]
    n_k = calculate_n_k(z[1:i, 1])
    (n - 1)/(n - 1 + alpha)
    alpha/(n - 1 + alpha)

    m_k.i = matrix(0, nrow = K, ncol = d)
    for(k in 1:K){
        m_k.i[k, ] = n_k.i[k]/(n_k.i[k] + rho_0)*x_bar_k.i[k, ] + rho_0/(n_k.i[k] + rho_0)*mu_0
    }


    cat("number of iteration is:", 1, "\n")
    cat("class is:", z[, 1], "\n")
    for(s in 2:iter_max){
        for(i in 1:n){
            n_k = calculate_n_k(z[, s-1], K = K)
            n_k.i = calculate_n_k(z[-i, s-1], K = K)
            x_bar_k = calculate_x_bar_k(x = x, z = z[, s-1], K = K)
            x_bar_k.i = calculate_x_bar_k(x = x[-i, ], z = z[-i, s-1], K = K)
            m_k.i = matrix(0, nrow = K, ncol = d)
            for(k in 1:K){
                m_k.i[k, ] = n_k.i[k]/(n_k.i[k] + rho_0)*x_bar_k.i[k, ] + rho_0/(n_k.i[k] + rho_0)*mu_0
            }
            a_n.i = 2*a_0 + (n-1)*d
            b_n.i_1st_term = 2*b_0
            x_bar = calculate_x_bar(x)
            out = rep(0, n)
            for(ii in 1:n){
                out[ii] = t(x[ii, ] - x_bar)%*%(x[ii, ] - x_bar)
            }
            b_n.i_2nd_term = sum(out[-i])
            x_bar.i = calculate_x_bar(x[-i, ])
            b_n.i_3rd_term = ((n-1)*rho_0/(rho_0 + n - 1))*t(x_bar.i - mu_0)%*%(x_bar.i - mu_0)
            b_n.i = b_n.i_1st_term + b_n.i_2nd_term + b_n.i_3rd_term
            b_n.i = c(b_n.i)
            sigma = (1 + (1/(n - 1 + rho_0)))*b_n.i*diag(d)
            out = rep(0, K)
            for(k in 1:K){
                out[k] = dmvt(x = x[i, ], delta = m_k.i[k, ], df = a_n.i, sigma = solve(sigma), type = "shifted", log = FALSE)*(n_k.i[k] + alpha)/sum(n_k.i + alpha)
            }a
            if(sum(out) == 0){
                p = rep(1/K, K)
            } else {
                p = rep(0, K)
                for(k in 1:K){
                    p[k] = out[k]/sum(out)
                }
            }
            z[i, s] = rmulti(n = 1, p = p)
        }
        cat("number of iteration is:", s, "\n")
        cat("class is:", z[, s], "\n")
    }
    z_map = z_map(z = z[, ], K = 3)
    out = list(z = z[, burn_in:iter_max], z_map = z_map, x = x)
}

data(iris)
x = as.matrix(iris[, 1:4])
fit_gmm = gaussian_mixture_model(x = x, K = 3, mu_0 = colMeans(x), rho_0 = 2)
z_true = iris[, 5]
table(z_map, z_true)

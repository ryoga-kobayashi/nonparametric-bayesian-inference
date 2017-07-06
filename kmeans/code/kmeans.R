rm(list = ls(all = TRUE))

optimize_mu = function(z, K, d){
    n_k = rep(0, K)
    mu = matrix(0, nrow = d, ncol = K)
    for(k in 1:K){
        n_k[k] = length(z[z==k])
    }
    for(j in 1:d){
        for(k in 1:K){
            mu[j, k] = (1/n_k[k])*sum(x[z==k, j])
        }
    }
    return(mu)
}

optimize_z = function(x, z, K, mu){
    n = nrow(x)
    for(i in 1:n){
        opt = rep(0, K)
        for(k in 1:K){
            opt[k] = t(x[i, ] - mu[, k])%*%(x[i, ] - mu[, k])
        }
        z[i] = (1:K)[opt==min(opt)]
    }
    return(z)
}

loss_function = function(x, z, K, mu){
    n = nrow(x)
    out = 0
    for(i in 1:n){
        out = out + t(x[i, ] - mu[, z[i]])%*%(x[i, ] - mu[, z[i]])
    }
    return((1/n)*out)
}

kmeans = function(x, K, iter_max = 100, epsilon = 0.0000001){
    n = nrow(x)
    d = ncol(x)
    l = rep(0, iter_max+1)
    z = sample(x = 1:K, size = n, replace = TRUE)
    l[1] = 0
    for(i in 1:iter_max){
        mu = optimize_mu(z = z, K = K, d = d)
        z = optimize_z(x = x, z = z, K = K, mu = mu)
        l[i+1] = loss_function(x = x, z = z, K = K, mu = mu)
        if(abs(l[i+1] - l[i]) < epsilon){
            iter_max = i
            message("converged")
            break
        }
        cat("iteration number is", i, "\n")
    }
    out = list(x = x, z = z, l = l[2:iter_max], iter_num = iter_max)
    return(out)
}

data(iris)
x = as.matrix(iris[, 1:4])
z_true = iris[, 5]
fit_kmeans = kmeans(x = x, K = 3)
z_fit = fit_kmeans$z
table(z_true, z_fit)




generate_datasets <- function(N, cov_matrix, n_repeats, sig_fraction, sig_means, null_means) {
  
  
  list_of_D <- list()  
  
  for(i in 1:n_repeats){
    
    
    rn <- runif(1)
    
    sig_hit <- rn < sig_fraction
    M <- if (sig_hit) {
      rmvnorm(N, sig_means, cov_matrix)
    } else {
      rmvnorm(N, null_means, cov_matrix)
    }
    
    D_cov <- data.frame(
      dv = as.vector(t(M)),
      t = factor(rep(1:4, N)),
      person = factor(rep(1:N, each = 4)), 
      sh  = rep(sig_hit, length(4*N))
    )
    
    list_of_D[[i]] <- D_cov
  }
  return(list_of_D)
}



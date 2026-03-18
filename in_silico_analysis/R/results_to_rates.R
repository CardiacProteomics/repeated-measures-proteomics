

get_rates <- function(df_in, label){


  cols <- colnames(df_in)
  pv_labels <- cols[grepl("pv_", cols)]
  
  
  sh <- df_in$sig_hit
  null <- !(df_in$sig_hit)
  
  
  df_out <- data.frame( n=length(sh), n_sig = sum(sh), n_nsig=sum(null) )
  
  for(label in pv_labels){
    
    tpr <- mean((df_in[[label]]<0.05)[sh])
    fpr <- mean((df_in[[label]]<0.05)[null])
    
    df_out[[paste0("fpr_", label)]] <- fpr
    df_out[[paste0("tpr_", label)]] <- tpr
    
  } 
  
  return(df_out)
  
}




get_power <- function( cor_matrix, sds, sig_means, n_subjects ){
  l_power <- c()
  
  for(i in 1:3)  
    for(j in (i+1):4){
      
      s1 <- sds[i]
      s2 <- sds[j]
      
      cc <- cor_matrix[i, j]
      
      sd_pooled <- (s1^2 + s2^2 - 2*cc*s1*s2)^0.5
      mean_diff <- sig_means[j]-sig_means[i]
      #effect_size <- mean_diff/sd_pooled
      
      
      res <- power.t.test(n = n_subjects, delta = mean_diff, sd = sd_pooled, sig.level = 0.05, type = "paired", alternative = "two.sided", strict = T)
      l_power <- c(l_power, res$power)
      
      
    }
  
  
  l_power <- c(l_power, NA)
  
  return(l_power)
}


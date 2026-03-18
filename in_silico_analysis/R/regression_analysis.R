



run_tests <- function(list_of_D){
  
  methods <- c("UN", "CS")
  
  contrasts <- c(t12 = "t2=0", 
                 t13 = "t3=0", 
                 t14 = "t4=0", 
                 t23 = "t3-t2=0", 
                 t24 = "t4-t2=0", 
                 t34 = "t4-t3=0")
  
  results <- vector("list", length(list_of_D))
  
  n_sets <- length(list_of_D)
  
  for(i in 1:n_sets){
  
    df_i <- list_of_D[[i]]
    
    row <- data.frame(sh = unique(df_i$sh))
    
    
    for(method in methods){
  
      fit <- lmm(dv ~ t, repetition = ~t|person, structure = method, data = df_i)
      an <- anova(fit)
      
      row[[paste0("pv_", tolower(method))]] <- an$multivariate$p.value
      
      for(name in names(contrasts)){
        row[[paste0("pv_",name,"_", tolower(method))]] <- anova(fit, effects = contrasts[[name]])$multivariate$p.value
      }

    } # methods
    results[[i]] <- row
    } 
  
  df_results <- do.call(rbind, lapply(results, as.data.frame))
  
  return(df_results)
  
}



run_tests_withPBS <- function(df_list, N, cov_matrix, n_repeats, sig_fraction, sig_means, null_means, F_repeats, withCS) {
  
  
  # Set up parallel backend
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
  contrasts <- c(t12 = "t2=0", 
                 t13 = "t3=0", 
                 t14 = "t4=0", 
                 t23 = "t3-t2=0", 
                 t24 = "t4-t2=0", 
                 t34 = "t4-t3=0")
  
  if(withCS){
    methods <- c("UN", "CS")
  }else{
    methods <- c("UN")
  }
  
  list_of_dfs <- list()
  
  
  results_list <- foreach(i = 1:n_repeats, .packages = c("mvtnorm", "LMMstar")) %dopar% {
    tryCatch({
      df_i <- df_list[[i]]
      sig_hit <- unique(df_i$sh)
      
      df_out <- data.frame(sig_hit)
      
      for(method in methods){
        
        fit <- lmm(dv ~ t, repetition = ~t | person, structure = method, data = df_i)
        an <- anova(fit)
        
        # Bootstrap null distribution under UN
        cp_i <- lmm(dv ~ t, repetition = ~t | person, structure = method, data = df_i, df = FALSE, method = "REML")
        cp2_i <- lmm(dv ~ 1, repetition = ~t | person, structure = method, data = df_i, df = FALSE, method = "REML")
        cov_mat_est <- sigma(cp_i)
        mean_est <- rep(as.numeric(cp2_i$param[1]), 4)
        
        F_val <- rep(NA, F_repeats)
        
        for (rep in 1:F_repeats) {
          try({
            M_rep <- rmvnorm(N, mean_est, cov_mat_est)
            df_rep <- data.frame(
              dv = as.vector(t(M_rep)),
              t = factor(rep(1:4, N)),
              person = factor(rep(1:N, each = 4))
            )
            
            cp_rep <- lmm(dv ~ t, repetition = ~t | person, structure = method, data = df_rep, df = F, method = "REML")
            q_rep <- anova(cp_rep)
            F_val[rep] <- q_rep$multivariate$statistic
            
          })
        }
        F_val_real <- na.omit(F_val)
        F_org <- anova(cp_i)$multivariate$statistic
        pv_boot <- sum(F_val_real >= F_org) / length(F_val_real)
        
        list_of_pvalues <- list()
        # Timepoint-specific p-values
        for( ic in 1:length(contrasts)){
          label_ic <- names(contrasts)[ic]
          contrast_ic <- contrasts[ic]
          
          list_of_pvalues[[paste0("pv_", label_ic, "_", method)]] <- anova(fit, effects = c(contrast_ic))$multivariate$p.value
          
        }
        
        df_out <- cbind(df_out, as.data.frame(list_of_pvalues))
        df_out[[paste0("pv_F_ordinary_", method)]] <- an$multivariate$p.value
        df_out[[paste0("pv_F_bootstrap_", method)]] <- pv_boot
        
        
      } #method 
      return(df_out)
      
      
    }, error = function(e) {
      message("Error in iteration ", i, ": ", e$message)
      return(NULL)
    })
  }
  
  df_out <- do.call(rbind, results_list)
  
  stopCluster(cl)
  return(df_out)
}











run_tests_small_sample <- function(df_list, N, cov_matrix, n_repeats, sig_fraction, sig_means, null_means,
                                   b_missing, missing_percentage, F_repeats) {
  
  
  # Set up parallel backend
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
  results_list <- foreach(i = 1:n_repeats, .packages = c("mvtnorm", "LMMstar", "mmrm"), .combine = rbind) %dopar% {
    tryCatch({
      df_i <- df_list[[i]]
      sig_hit <- unique(df_i$sh)
      
      cp_un <- lmm(dv ~ t, repetition = ~t | person, structure = "UN", data = df_i)
      cp_cs <- lmm(dv ~ t, repetition = ~t | person, structure = "CS", data = df_i)
      
      q_un <- anova(cp_un)
      q_cs <- anova(cp_cs)
      
      # Bootstrap null distribution under UN
      cp_un_i <- lmm(dv ~ 1, repetition = ~t | person, structure = "UN", data = df_i, df = FALSE, method = "REML")
      cp_un2_i <- lmm(dv ~ 1, repetition = ~t | person, structure = "UN", data = df_i, df = FALSE, method = "REML")
      
      cov_mat_est <- sigma(cp_un_i)
      mean_est <- rep(as.numeric(cp_un2_i$param[1]), 4)
      
      #F_repeats <- 1000
      F_val <- numeric(F_repeats)
      
      for (rep in 1:F_repeats) {
        try({
          M_rep <- rmvnorm(N, mean_est, cov_mat_est)
          df_rep <- data.frame(
            dv = as.vector(t(M_rep)),
            t = factor(rep(1:4, N)),
            person = factor(rep(1:N, each = 4))
          )
          cp_un_rep <- lmm(dv ~ t, repetition = ~t | person, structure = "UN", data = df_rep,
                           df = F, method = "REML")
          q_rep <- anova(cp_un_rep)
          F_val[rep] <- q_rep$multivariate$statistic
        })
      }
      
      F_org <- q_un$multivariate$statistic
      pv_boot_star <- sum(F_val >= F_org, na.rm = TRUE) / sum(!is.na(F_val))
      
      # Timepoint-specific p-values
      pv_t12_un <- anova(cp_un, effects = c("t2=0"))$multivariate$p.value
      pv_t13_un <- anova(cp_un, effects = c("t3=0"))$multivariate$p.value
      pv_t14_un <- anova(cp_un, effects = c("t4=0"))$multivariate$p.value
      pv_t23_un <- anova(cp_un, effects = c("t3-t2=0"))$multivariate$p.value
      pv_t24_un <- anova(cp_un, effects = c("t4-t2=0"))$multivariate$p.value
      pv_t34_un <- anova(cp_un, effects = c("t4-t3=0"))$multivariate$p.value
      
      
      pv_t12_cs<- anova(cp_cs, effects = c("t2=0"))$multivariate$p.value
      pv_t13_cs<- anova(cp_cs, effects = c("t3=0"))$multivariate$p.value
      pv_t14_cs<- anova(cp_cs, effects = c("t4=0"))$multivariate$p.value
      pv_t23_cs<- anova(cp_cs, effects = c("t3-t2=0"))$multivariate$p.value
      pv_t24_cs<- anova(cp_cs, effects = c("t4-t2=0"))$multivariate$p.value
      pv_t34_cs<- anova(cp_cs, effects = c("t4-t3=0"))$multivariate$p.value
      
      
      
      df_wide <- reshape(df_i, 
                         timevar = "t", 
                         idvar = "person", 
                         direction = "wide")
      
      
      tt12 <- t.test(df_wide$dv.1, df_wide$dv.2, paired = TRUE, var.equal = FALSE, alternative="two.sided")
      tt13 <- t.test(df_wide$dv.1, df_wide$dv.3, paired = TRUE, var.equal = FALSE, alternative="two.sided")
      tt14 <- t.test(df_wide$dv.1, df_wide$dv.4, paired = TRUE, var.equal = FALSE, alternative="two.sided")
      tt23 <- t.test(df_wide$dv.2, df_wide$dv.3, paired = TRUE, var.equal = FALSE, alternative="two.sided")
      tt24 <- t.test(df_wide$dv.2, df_wide$dv.4, paired = TRUE, var.equal = FALSE, alternative="two.sided")
      tt34 <- t.test(df_wide$dv.3, df_wide$dv.4, paired = TRUE, var.equal = FALSE, alternative="two.sided")
      
      pairtt12 <-  tt12$p.value
      pairtt13 <-  tt13$p.value
      pairtt14 <-  tt14$p.value
      pairtt23 <-  tt23$p.value
      pairtt24 <-  tt24$p.value
      pairtt34 <-  tt34$p.value
      
      
      data.frame(
        pv_un = q_un$multivariate$p.value,
        pv_cs = q_cs$multivariate$p.value,
        sh_out = sig_hit,
        pv_t12_un = pv_t12_un,
        pv_t13_un = pv_t13_un,
        pv_t14_un = pv_t14_un,
        pv_t23_un = pv_t23_un,
        pv_t24_un = pv_t24_un,
        pv_t34_un = pv_t34_un,
        pv_t12_cs = pv_t12_cs,
        pv_t13_cs = pv_t13_cs,
        pv_t14_cs = pv_t14_cs,
        pv_t23_cs = pv_t23_cs,
        pv_t24_cs = pv_t24_cs,
        pv_t34_cs = pv_t34_cs,
        pairtt12, 
        pairtt13, 
        pairtt14, 
        pairtt23,
        pairtt24,
        pairtt34,
        pv_boot_star = pv_boot_star, 
        row.names = i
      )
      
      
    }, error = function(e) {
      message("Error in iteration ", i, ": ", e$message)
      na_row <- data.frame(
        pv_un = NA, pv_cs = NA, sh_out = NA,
        pv_t12_un = NA, pv_t13_un = NA, pv_t14_un = NA, pv_t23_un = NA,pv_t24_un = NA,pv_t34_un = NA,
        pv_t12_cs = NA, pv_t13_cs = NA, pv_t14_cs = NA, pv_t23_cs = NA,pv_t24_cs = NA,pv_t34_cs = NA,
        pv_boot_star = NA,
        row.names = i
      )
      return(na_row)
    })
  }
  
  stopCluster(cl)
  return(results_list)
}




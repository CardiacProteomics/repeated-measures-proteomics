

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


results_to_rates_small_sample <- function(df_res){
  
  n_sig <- sum(df_res$sh)
  n_nsig <- sum(!df_res$sh)
  n_tot <- length(df_res$sh)
  
  df_res$sh <- df_res$sh_out
  
  
  tpr_un_boot <- sum(as.numeric(df_res$pv_boot[!is.na(df_res$pv_boot)])<=0.05 & df_res$sh_out[!is.na(df_res$pv_boot)])/sum(df_res$sh_out[!is.na(df_res$pv_boot)])
  tpr_un_F <- sum(df_res$pv_un<=0.05 & df_res$sh)/sum(df_res$sh)
  tpr_un_t2 <- sum(df_res$pv_t2_un<=0.05 & df_res$sh)/sum(df_res$sh)
  tpr_un_t3 <- sum(df_res$pv_t3_un<=0.05 & df_res$sh)/sum(df_res$sh)
  tpr_un_t4 <- sum(df_res$pv_t4_un<=0.05 & df_res$sh)/sum(df_res$sh)
  tpr_un_t5 <- sum(df_res$pv_t5_un<=0.05 & df_res$sh)/sum(df_res$sh)
  #tpr_un_boot_star <- sum(as.numeric(df_res$pv_boot[!is.na(df_res$pv_boot)])<=0 & df_res$sh_out) 
  #tpr_un_boot_mmrm <- sum(as.numeric(df_res$pv_boot_mmrm[!is.na(df_res$pv_boot_mmrm)])==0 & df_res$sh_out) 
  
  tpr_cs_F <- sum(df_res$pv_cs<=0.05 & df_res$sh)/sum(df_res$sh)
  tpr_cs_t2 <- sum(df_res$pv_t2_cs<=0.05 & df_res$sh)/sum(df_res$sh)
  tpr_cs_t3 <- sum(df_res$pv_t3_cs<=0.05 & df_res$sh)/sum(df_res$sh)
  tpr_cs_t4 <- sum(df_res$pv_t4_cs<=0.05 & df_res$sh)/sum(df_res$sh)
  tpr_cs_t5 <- sum(df_res$pv_t5_cs<=0.05 & df_res$sh)/sum(df_res$sh)
  
  
  fpr_un_F <- sum(as.numeric(df_res$pv_un[!is.na(df_res$pv_un)])<=0.05 & !df_res$sh_out[!is.na(df_res$pv_un)])/sum(!df_res$sh)
  fpr_un_t2 <- sum(df_res$pv_t2_un<=0.05 & !df_res$sh)/sum(!df_res$sh)
  fpr_un_t3 <- sum(df_res$pv_t3_un<=0.05 & !df_res$sh)/sum(!df_res$sh)
  fpr_un_t4 <- sum(df_res$pv_t4_un<=0.05 & !df_res$sh)/sum(!df_res$sh)
  fpr_un_t5 <- sum(df_res$pv_t5_un<=0.05 & !df_res$sh)/sum(!df_res$sh)
  fpr_un_boot <- sum(as.numeric(df_res$pv_boot_star[!is.na(df_res$pv_boot_star)])<=0.05 & !df_res$sh_out[!is.na(df_res$pv_boot_star)])/sum(!df_res$sh_out[!is.na(df_res$pv_boot_star)])
  
  
  #fpr_un_boot_mmrm <- sum(as.numeric(df_res$pv_boot_mmrm[!is.na(df_res$pv_boot_mmrm)])<=0.05 & !df_res$sh_out[!is.na(df_res$pv_boot_mmrm)])/sum(!df_res$sh)
  
  
  fpr_tt_t2 <- sum(df_res$pairtt2<=0.05 & !df_res$sh)/sum(!df_res$sh)
  fpr_tt_t3 <- sum(df_res$pairtt3<=0.05 & !df_res$sh)/sum(!df_res$sh)
  fpr_tt_t4 <- sum(df_res$pairtt4<=0.05 & !df_res$sh)/sum(!df_res$sh)
  fpr_tt_t5 <- sum(df_res$pairtt5<=0.05 & !df_res$sh)/sum(!df_res$sh)
  
  tpr_tt_t2 <- sum(df_res$pairtt2<=0.05 & df_res$sh)/sum(df_res$sh)
  tpr_tt_t3 <- sum(df_res$pairtt3<=0.05 & df_res$sh)/sum(df_res$sh)
  tpr_tt_t4 <- sum(df_res$pairtt4<=0.05 & df_res$sh)/sum(df_res$sh)
  tpr_tt_t5 <- sum(df_res$pairtt5<=0.05 & df_res$sh)/sum(df_res$sh)
  
  fpr_cs_F <- sum(df_res$pv_cs<=0.05 & !df_res$sh)/sum(!df_res$sh)
  fpr_cs_t2 <- sum(df_res$pv_t2_cs<=0.05 & !df_res$sh)/sum(!df_res$sh)
  fpr_cs_t3 <- sum(df_res$pv_t3_cs<=0.05 & !df_res$sh)/sum(!df_res$sh)
  fpr_cs_t4 <- sum(df_res$pv_t4_cs<=0.05 & !df_res$sh)/sum(!df_res$sh)
  fpr_cs_t5 <- sum(df_res$pv_t5_cs<=0.05 & !df_res$sh)/sum(!df_res$sh)
  
  #fpr_dream <- sum(df_res$pv_dream<=0.05 & !df_res$sh)/sum(!df_res$sh)
  #tpr_dream <- sum(df_res$pv_dream<=0.05 & df_res$sh)/sum(df_res$sh)
  #fpr_limma <- sum(df_res$pv_limma<=0.05 & !df_res$sh)/sum(!df_res$sh)
  #tpr_limma <- sum(df_res$pv_limma<=0.05 & df_res$sh)/sum(df_res$sh)
  
  data.frame(
    fpr_un = c(fpr_un_t2, fpr_un_t3, fpr_un_t4, fpr_un_t5, fpr_un_F, fpr_un_boot),
    fpr_cs = c(fpr_cs_t2, fpr_cs_t3, fpr_cs_t4, fpr_cs_t5, fpr_cs_F, "Not performed"),
    fpr_tt = c(fpr_tt_t2, fpr_tt_t3, fpr_tt_t4, fpr_tt_t5, "Not performed", "Not performed"), 
    tpr_un = c(tpr_un_t2, tpr_un_t3, tpr_un_t4, tpr_un_t5, tpr_un_F, tpr_un_boot),
    tpr_cs = c(tpr_cs_t2, tpr_cs_t3, tpr_cs_t4, tpr_cs_t5, tpr_cs_F, "Not performed"), 
    tpr_tt = c(tpr_tt_t2, tpr_tt_t3, tpr_tt_t4, tpr_tt_t5, "Not performed", "Not performed"), 
    extra_info = c(n_tot, n_sig, n_nsig, "Not performed", "Not performed", "Not performed")
  )
}




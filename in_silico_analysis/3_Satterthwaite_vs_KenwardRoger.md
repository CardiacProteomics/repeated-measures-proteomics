``` r
library(ggplot2)
library(ggpubr)
library(MASS)
library(LMMstar)
library(mvtnorm)
library(lavaSearch2)
library(nlme)
library(foreach)
library(doParallel)
library(compiler)
library(foreach)
library(doParallel)
library(mmrm)
library(here)
library(MBESS)
library(mvtnorm)
library(MASS)
rm(list = ls())

here::i_am("in_silico_analysis/3_Satterthwaite_vs_KenwardRoger.Rmd")
source_path <- here("R")

source(paste0(source_path, "/generate_datasets.R"))
```

``` r
cov_matrix <- toeplitz(c(1,0.8,0.5, 0.2))
sig_means <- c(1, 1.5, 1.5, 1.5) # 0.35
null_means <- c(1 ,1 ,1, 1)

sig_fraction <- 0.5

n_repeats <- 10

for(j in 1:10){

sh <- c()
pv_KR <- c()
pv_S <- c()
nn <- c()

for(N in c(7,8,9,10,11,12,13,14,15)){
  
list_of_dfs <- generate_datasets(N, cov_matrix, n_repeats, sig_fraction, sig_means, null_means)

for(df_cov in list_of_dfs){
try({


  fit_KR <- mmrm(
    formula = dv ~ t + us(t | person),
    data = df_cov,
    optimizer = "BFGS", 
    method = "Kenward-Roger"
  )
  
  fit_S <- mmrm(
    formula = dv ~ t + us(t | person),
    data = df_cov,
    optimizer = "BFGS", 
    method = "Satterthwaite"
  )
  
  
  q_KR <- car::Anova(fit_KR, type="II")
  q_S <- car::Anova(fit_S, type="II")
  
  pv_KR <- c(pv_KR, q_KR$`Pr(>F)`)
  pv_S <- c(pv_S, q_S$`Pr(>F)`)
  
  nn <- c(nn, N)
  
  sh <- c(sh, unique(df_cov$sh))

})
}
}

df_res <- data.frame(nn, sh, pv_KR, pv_S)

rn <- floor( runif(1, 1, 1e9) )

data_path <- here(  "Numerical_data", "S_vs_KR")

saveRDS(df_res, paste0(data_path, "/df_SvsKR_", rn ,".rds"))

}
```

``` r
data_path <- here(  "Numerical_data", "S_vs_KR")
folder <- data_path
files <- list.files(path = folder)

res <- lapply( files, function(file){
  
      d_in <- readRDS(paste0(data_path,"/", file))
      
    
}) # lapply 

df_res <- do.call(rbind, res)


df_plot <- df_res %>%
      dplyr::group_by(nn) %>%
      dplyr::summarise(
        
        tpr_KR = sum(pv_KR<=0.05 & sh)/sum(sh), 
        tpr_S = sum(pv_S<=0.05 & sh)/sum(sh), 
        fpr_KR = sum(pv_KR<=0.05 & !sh)/sum(!sh), 
        fpr_S = sum(pv_S<=0.05 & !sh)/sum(!sh), 
        .groups = "drop"
      ) %>%
  dplyr::rename(n_person = nn)


g1 <- ggplot(df_plot) + 
  geom_point(aes(x = n_person, y = fpr_KR, color= 'Kenward-Roger' )) + 
  geom_line(aes(x = n_person, y = fpr_KR, color= 'Kenward-Roger' )) + 
  geom_point(aes(x = n_person, y = fpr_S, color= 'Satterthwaite' )) + 
  geom_line(aes(x = n_person, y = fpr_S, color= 'Satterthwaite' )) + 
  scale_color_manual( 
  name = "", 
  values = c( 
           "Kenward-Roger" = "steelblue", 
           "Satterthwaite" = "red4"))  + 
  theme_bw() + 
  theme(
    legend.text = element_text(size = 14)
  ) + 
  labs(x = "Number of subjects", y = "False positive rate", title = "False positive rate")+ 
  geom_hline(yintercept = 0.05, linetype = "dashed") 

g2 <- ggplot(df_plot) + 
  geom_point(aes(x = n_person, y = tpr_KR, color= 'Kenward-Roger' )) + 
  geom_line(aes(x = n_person, y = tpr_KR, color= 'Kenward-Roger' )) + 
  geom_point(aes(x = n_person, y = tpr_S, color= 'Satterthwaite' )) + 
  geom_line(aes(x = n_person, y = tpr_S, color= 'Satterthwaite' )) + 
  scale_color_manual( 
  name = "", 
  values = c( 
           "Kenward-Roger" = "steelblue", 
           "Satterthwaite" = "red4"))  + 
  theme_bw() + 
  theme(
    legend.text = element_text(size = 14)
  ) + 
  labs(x = "Number of subjects", y = "True positive rate", title = "Power")


g <- ggarrange(g1,g2, common.legend = TRUE, legend = "bottom")



data_path <- here( "Figures")
ggsave(plot = g, filename= paste0(data_path, "/sFig_KR_vs_S.pdf"), width=8, height=4)
```

    ## systemfonts and textshaping have been compiled with different versions of Freetype. Because of this, textshaping will not use the font cache provided by systemfonts

``` r
g
```

![](3_Satterthwaite_vs_KenwardRoger_files/figure-gfm/plot%20analysis-1.png)<!-- -->

# Analysis of tumor proteomes with both a time and space component

This dataset contains samples from 34 patients with breast cancer. For
each patient there is sample of the tumor prior to treatment, a sample
post treatment and a sample from healthy tissue close to the tumor.
Since there are three measurements from the same patient, the
measurements are not independent and should be analyzed with
repeated-measures statistics. First the publicly data is imported and
prepared for analysis.

    library(openxlsx)
    library(ggplot2)
    library(stringr)
    library(factoextra)
    library(cluster)
    library(pheatmap)
    library(ggpubr)
    library(ggpointdensity)
    library(LMMstar)
    library(here)
    #source("XXXX/functional_analysis_functions.R") ## importing homemade overrepresentation analysis function

    ## please access the data from: 

    # Import data
    M <- read.xlsx("D:/Project_longitudinal_methods/Geiger_cancer_cells/msb209443-sup-0003-datasetev1.xlsx", sheet="EV1", startRow = 3)
    D <- read.xlsx("D:/Project_longitudinal_methods/Geiger_cancer_cells/msb209443-sup-0004-datasetev2.xlsx", sheet="EV2", startRow = 2)

    # Extraction of data columns
    df_dat <- D[, grepl("Ratio", colnames(D))]             # retrieve proteomics data from dataframe
    rownames(df_dat) <- D$Majority.protein.IDs             # add protein group as row name
    df_dat[] <- lapply(df_dat, function(x) as.numeric(x))  # make the data numeric
    df_dat[df_dat=="NaN"] <- NA                            #Replace string by actual missing value

    # Maps from patient id to clinical characteristics
    patient2tumor_size  <- setNames(M$Tumor.size, M$Patient.code )         
    patient2lymph_node  <- setNames(M$Lymph.node.involved, M$Patient.code)
    patient2age         <- setNames(M$Age, M$Patient.code)

    # Construct dataframe with meta data (df_meta) 
    # Make sure that the data is aligned with df_dat

    patient <- c()
    tissue <- c()
    col_keep <- c()
    tumor_size <- c()
    lymph_node <- c()
    age <- c()

    for(col_i in colnames(df_dat)){
      
      tmp <- strsplit(col_i, "normalized.")[[1]][2]
      num <- str_extract(tmp, "\\d+")
      chars <- str_extract(tmp, "[A-Za-z]+")  # Extract alphabetic part
      
      
      if(nchar(chars)==1 & !(num %in% c(34,35,36,37,38)) ){
      
        patient <- c(patient, num)
        tissue <- c(tissue, chars)
        col_keep <- c(col_keep, col_i)
        tumor_size <- c(tumor_size, unname(patient2tumor_size[num]))
        lymph_node <- c(lymph_node, unname(patient2lymph_node[num]))
        age <- c(age, unname(patient2age[num]))

      }
    }


    df_meta <- data.frame(patient, tissue, tumor_size, lymph_node, age)
    df_dat <- df_dat[, col_keep]

    df_fit <- df_meta


    rows <- rownames(df_dat)


    res <- lapply(rows, function(row){
      
      y <- as.numeric(df_dat[row, ])
      
      if(mean(!is.na(y))<0.7) return(NULL)  # filter out proteins with more than 30% missing values
      
        df_fit$y <- y # add the protein intensity to the dataframe for fitting 
        fit <- lmm(y~tissue + lymph_node, repetition = ~tissue|patient, structure = "UN", data=df_fit) # run the covariance pattern analysis/regression 
        an <- anova(fit) # run F-tests
        su <- summary(fit)
        
        pv_all <- an$multivariate$p.value[1] # p-value for change across time and tissue
        pv_ln <-  an$multivariate$p.value[2] # p-value for lymph_node
        
        ## p-values for all pairwise comparisons of time and tissue 
        pv_tt_AvsB <- anova(fit, effects = c("tissueB = 0"))$multivariate$p.value
        pv_tt_AvsC <- anova(fit, effects = c("tissueC = 0"))$multivariate$p.value
        pv_tt_BvsC <- anova(fit, effects = c("tissueB - tissueC = 0"))$multivariate$p.value
        
        ## parameters for all pairwise comparisons of time and tissue 
        par_AvsB <- su$mean['tissueB', 'estimate']
        par_AvsC <- su$mean['tissueC', 'estimate']
        par_BvsC <- su$mean['tissueC', 'estimate'] -  su$mean['tissueB', 'estimate']
        
        data.frame(
          row = row, 
          pv_all = pv_all, 
          pv_ln = pv_ln, 
          pv_tt_AvsB = pv_tt_AvsB, 
          pv_tt_AvsC = pv_tt_AvsC, 
          pv_tt_BvsC = pv_tt_BvsC, 
          par_AvsB = par_AvsB, 
          par_AvsC = par_AvsC, 
          par_BvsC = par_BvsC
        )
        

    }) #lapply

    df_tt <- do.call(rbind, res)


    df_tt$min_par <- apply(df_tt[, grepl("par", colnames(df_tt))],1, min)
    df_tt$max_par <- apply(df_tt[, grepl("par", colnames(df_tt))],1, max)

    df_tt$fdr_AB <- p.adjust(df_tt$pv_tt_AvsB, method = "BH")
    df_tt$fdr_AC <- p.adjust(df_tt$pv_tt_AvsC, method = "BH")
    df_tt$fdr_BC <- p.adjust(df_tt$pv_tt_BvsC, method = "BH")

    df_tt$sig_AB <- df_tt$fdr_AB<=0.05
    df_tt$sig_AC <- df_tt$fdr_AC<=0.05
    df_tt$sig_BC <- df_tt$fdr_BC<=0.05


    data_path <- here(  "Results_real_world", "tumor")
    saveRDS(df_tt, file = paste0(data_path,"/df_tt_geiger.rds"))

    data_path <- here(  "Results_real_world", "tumor")

    df_tt <- readRDS(df_tt, file = paste0(data_path, "/df_tt_geiger.rds"))



    ymax <- max(-log10(df_tt$min))

    g1 <- ggplot(df_tt[!df_tt$sig_AB, ]) + geom_pointdensity(aes(x=par_AvsB, y = -log10(pv_tt_AvsB)),alpha=0.5, shape = 16, size = 1.4) + 
      geom_point(data = df_tt[df_tt$sig_AB, ], aes(x=par_AvsB, y = -log10(pv_tt_AvsB)), color='black', fill = "brown2", alpha=1, shape = 21, size = 1.4) + 
      theme_bw() + ylim(0, ymax) + 
      labs(y = "p-value [-log10]", x="log2(pre_treat) - log2(post_treat)", 
           title = paste0("Down=", 
                     sum(df_tt$sig_AB & df_tt$par_AvsB<0), 
                     " Up=", 
                     sum(df_tt$sig_AB & df_tt$par_AvsB>0)
                     )) + theme(legend.position = "none") 

    g2 <- ggplot(df_tt[!df_tt$sig_AC, ]) + geom_pointdensity(aes(x=par_AvsC, y = -log10(pv_tt_AvsC)), alpha=0.5, shape = 16, size = 1.4) + 
      geom_point(data = df_tt[df_tt$sig_AC, ], aes(x=par_AvsC, y = -log10(pv_tt_AvsC)), color='black', fill = "brown2", alpha=1, shape = 21, size = 1.4) + 
      theme_bw()+ ylim(0, ymax) + 
      labs(y ="p-value [-log10]",
                                       x = "log2(pre_treat) - log2(normal_tissue)", 
                                       title=paste0("Down=", 
                     sum(df_tt$sig_AC & df_tt$par_AvsC<0), 
                     " Up=", 
                     sum(df_tt$sig_AC & df_tt$par_AvsC>0)
                     ))+ theme(legend.position = "none") 

    g3 <- ggplot(df_tt[!df_tt$sig_BC, ]) + geom_pointdensity(aes(x=par_BvsC, y = -log10(pv_tt_BvsC)), alpha=0.5, shape = 16, size = 1.4) + 
      geom_point(data = df_tt[df_tt$sig_BC, ], aes(x=par_BvsC, y = -log10(pv_tt_BvsC)), color='black', fill = "brown2", alpha=1, shape = 21, size = 1.4) + 
      theme_bw()+ ylim(0, ymax) + 
      labs(y = "p-value [-log10]", 
           x = "log2(post_treat) - log2(normal_tissue)", 
           title = paste0("Down=", 
                     sum(df_tt$sig_BC & df_tt$par_BvsC<0), 
                     " Up=", 
                     sum(df_tt$sig_BC & df_tt$par_BvsC>0)
                     ))+ theme(legend.position = "none") 

    g <- ggarrange(g1,g2,g3, ncol=3)

    data_path <- here(  "Figures")

    ggsave(paste0(data_path, "/Geiger_VPs.pdf"), width=8, height=3)

    data_path <- here(  "Results_real_world", "tumor")

    df_tt <- readRDS(df_tt, file = paste0(data_path, "/df_tt_geiger.rds"))

    ind7 <- df_tt$sig_AB | df_tt$sig_AC | df_tt$sig_BC

    ind <- rownames(df_dat) %in% df_tt[ind7, ]$row


    df_dat_sig <- df_dat[ind, ]

    mean_A <- c()
    mean_B <- c()
    mean_C <- c()

    rn <- c()
    for(row in rownames(df_dat_sig)){
      
      rn <- c(rn, row)
      mean_A <- c(mean_A, mean(as.numeric(df_dat_sig[row, df_meta$tissue == "A"]), na.rm=T))
      mean_B <- c(mean_B, mean(as.numeric(df_dat_sig[row, df_meta$tissue == "B"]), na.rm=T))
      mean_C <- c(mean_C, mean(as.numeric(df_dat_sig[row, df_meta$tissue == "C"]), na.rm=T))

    }

    df_mean <- data.frame(mean_A, mean_B, mean_C)

    df_sc <- t(scale(t(df_mean)))

    rownames(df_sc) <- rownames(df_dat_sig)
    data_path <- here(  "Figures")

    p <- pheatmap(df_sc, cutree_rows = 6, 
                  filename = paste0(data_path, "/pheatmap_geiger.pdf"), 
                  show_rownames = FALSE, show_colnames = FALSE)

    f <- fviz_nbclust(1-df_sc, hcut, method = "wss")

    ord <- p$tree_row$order
    clusters <- cutree(p$tree_row, k=6)
    arrangement <- unique(clusters[ord])


    term_clusters <- data.frame(
      #GO_Term = names(clusters),
      #geneset_size, 
      Cluster = clusters, 
      order = ord#, 
      #cluster_mean
    )


    scaled_df <- as.data.frame(lapply(df_dat_sig, function(x) {
      if(is.numeric(x)) {
        (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
      } else {
        x
      }
    }))

    rownames(scaled_df) <- rownames(df_dat_sig)



    count <- 0
    dd <- list()
    lop <- list()
    for(cluster_i in unique(clusters)){
      
      count <- count + 1
      
      y <- c()
      ti <- c()
      patient <- c()
      names <- rownames(term_clusters[term_clusters$Cluster==cluster_i, ])
      nn <- c()
      for(name_i in names){
        
        #ord <- order(df_fit$patient)
        
        #y <- c(y, as.numeric(scaled_df[name_i, ])[ord])
        #ti <- c(ti, df_fit$tissue[ord])
        #patient <- c(patient, df_fit$patient[ord])
        #nn <- c(nn, rep(name_i, length(ord)))
        
        y <- c(y, as.numeric(df_sc[name_i, ]))
        ti <- c(ti, c("A", "B", "C"))
        nn <- c(nn, rep(name_i,3))
        

        
      }
      
      df_plot <- data.frame(y,ti, nn)
      

      

      
      df_plot$nn <- factor(df_plot$nn, levels = unique(df_plot$nn))
      df_plot$ti <- factor(df_plot$ti, levels = c("A", "B", "C"))
      

     

      
       dd[[cluster_i]] <- df_plot
      
      #df_plot$patient <- as.factor(df_plot$patient)
      
      #df_plot <- df_plot %>%
      #arrange(patient, ti)
       
       df_me <- df_plot %>%
      group_by(ti) %>%
      dplyr::summarise(mean_y = mean(y, na.rm = TRUE)) %>%
      mutate(ti = factor(ti, levels = levels(df_plot$ti)))
       
       df_me <- as.data.frame(df_me)
       
       df_me$ti <- factor(df_me$ti, levels = c("A", "B", "C"))
       
    if(count == 0){
      g <- ggplot(df_plot) + 
        geom_line(aes(ti, y, group = nn), alpha = 0.1) + 
        geom_line(data = df_me, aes(ti, mean_y, group = 1), color = "red", size = 1.2) +  # Add mean trajectory
        labs(x = "", y = "standardized mean protein intensity") + 
        scale_x_discrete(labels = c("A" = "Tumor-pre", "B" = "Tumor-post", "C" = "Normal")) + 
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12 ), 
              panel.grid.major.y = element_blank(),
              panel.grid.minor = element_blank())
      
    } else {
      g <- ggplot(df_plot) + 
        geom_line(aes(ti, y, group = nn), alpha = 0.1) + 
        geom_line(data = df_me, aes(ti, mean_y, group = 1), color = "red", size = 1.2) +  # Add mean trajectory
        labs(x = "", y = NULL) + 
        scale_x_discrete(labels = c("A" = "Tumor-pre", "B" = "Tumor-post", "C" = "Normal")) + 
        theme_minimal() + 
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(), 
              axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
             panel.grid.major.y = element_blank(),
              panel.grid.minor = element_blank())
    }
      
      
      lop[[cluster_i]] <- g
     
    }


    g <- ggarrange(lop[[1]], lop[[2]], lop[[3]], lop[[4]], lop[[5]], lop[[6]], ncol = 6)


    data_path <- here(  "Figures")

    ggsave(plot = g, paste0(data_path, "/tumor_trajectories.pdf"), width=8, height=3)

    ### tmp

    universe <- D$Gene.names

    uid2gene <- setNames(D$Gene.names, D$Majority.protein.IDs)



    list_of_ora_art <- list()
    for(i in unique(clusters)){

      
      ind <- clusters == i 
      names_sub <- names(clusters)[ind]
      
      genes_sub <- c()
      for(name in names_sub){
        genes_sub <- c(genes_sub, uid2gene[[name]])
      }
      
      
      out1 <- genes_ORA_new(genes_sub, universe, "reactome_human", 5,200)

      
      
      list_of_ora_art[[i]] <- out1[out1$pv_bh<=0.05, ]$pw
      
      ora_plot <- out1[out1$pv_bh<=0.05, ]

      ora_plot <- ora_plot[order(ora_plot$gene_ratio), ]
      #ora_plot$pw_nice <- substr(ora_plot$pw_nice, 1, 50)
      #ora_plot$pw_nice <- factor(ora_plot$pw_nice, levels = ora_plot$pw_nice)
      
      
      
      
      #g <- ggplot(ora_plot, aes(x = gene_ratio, y = pw, size = gene_count, color = pv_bh)) +
      #  geom_point() +
      #  scale_color_gradient(low = "red", high = "blue") +
      #  theme_bw() +
      #  labs(title = paste0("ORA ART cluster:", i, "_size:", sum(ind)), x = "Gene Ratio", y = "Description", color = "Adjusted P-value", size = "Count")
      
      #ggsave(filename = paste0("E:/Project_longitudinal_methods/Figures/Geiger ORA/cluster_art_",i, ".png"), plot = g)
      

      
    }

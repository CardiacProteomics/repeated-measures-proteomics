# Repeated measures proteomics 
This project considers statistical methods for analyzing repeated-measures proteomics data. We present two exemplary analyses of publicly available repeated-measures proteomics datasets: one dataset with repeated measurements over time and one dataset that mixes time and space. The two analyses
serve as tutorials for repeated-measures analysis and are available in the folder "real_world_analysis". 

This repository is related to the manuscript "Statistical Approaches for Repeated-Measures Proteomics Data”. In addition to the two tutorials, it contains an analysis of advantages and pitfalls of various statistical repeated-measures methods. This is achieved by analysing in-silico generated idealized repeated-measures data. We introduce covariance pattern models for proteomics and demonstrate their advantages over existing methods. The main issue is to control the false positive rate at the intended level (typically 5%), which is often not the case with existing resources. The analyses can be found in the folder "in_silico_analysis".    
  
Lead developer: Mikkel Skjoldan Svenningsen

## Analysis of two real-world repeated-measures proteomics datasets. 
### Tumor dataset
 - Analysis of the publicly available dataset from "Proteomic patterns associated with response to breast cancer neoadjuvant treatment".
 - Data import and structuring
 - Repeated measures analysis using covariance pattern models with unstructured covariance
 - Volcano plots
 - Trajectory clustering
 - Functional analysis of trajectory clusters
### Weightloss dataset
  


## Analysis of in-silico results presented in the manuscript

### 0_test_performance.Rmd:
 - In silico generation of datasets resembling idealized repeated-measures proteomics datasets
 - Repeated-measures analysis of change over time/site (covariance pattern models)
    - Two covariance structures: Unstructured and compound symmetry    
 - Comparison of t-tests fit with ordinary F-test
 - Produces figure 2 from the manuscript

### 1_small_samples_F_test.Rmd: 
  - In silico generation of datasets resembling idealized repeated-measures proteomics datasets
     - Smaller sample sizes compared to previous analysis (0_test_performance) 
  - Repeated-measures analysis of change over time/site (covariance pattern models)
     - Two covariance structures: Unstructured and compound symmetry
  - Comparison of ordinary F-tests with F-tests based on bootstrapping
  - Produces figure 3 from the manuscript

### 2_small_samples_ttest.Rmd: 
  - In silico generation of datasets resembling idealized repeated-measures proteomics datasets
     - Sample sizes: 10-25  
  - Repeated-measures analysis of change over time/site (covariance pattern models)
     - One covariance structures: Unstructured 
  - Performance of T-test for small sample sizes
  - Produces figure 4 from the manuscript

## Analysis of in-silico results presented in the Supplementary Material

### 3_Satterthwaite_vs_KenwardRoger.Rmd:
  - In silico generation of repeated-measures data
  - Analysis using the mmrm package
  - Comparison of Satterthwaite and Kenward-Roger degrees of freedom estimate for small sample sizes

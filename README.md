# repeated-measures-proteomics
Tutorial for analysis of repeated measures proteomics. Related to the manuscript: "Statistical Approaches for Repeated-Measures Proteomics Data”.  

Lead developer: Mikkel Skjoldan Svenningsen

# 0_test_performance.Rmd:
 - In silico generation of datasets resembling idealized repeated-measures proteomics datasets
 - Repeated-measures analysis of change over time/site (covariance pattern models)
    - Two covariance structures: Unstructured and compound symmetry    
 - Comparison of t-tests fit with ordinary F-test
 - Produces figure 2 from the manuscript

# 1_small_samples_F_test.Rmd: 
  - In silico generation of datasets resembling idealized repeated-measures proteomics datasets
     - Smaller sample sizes compared to previous analysis (0_test_performance) 
  - Repeated-measures analysis of change over time/site (covariance pattern models)
     - Two covariance structures: Unstructured and compound symmetry
  - Comparison of ordinary F-tests with F-tests based on bootstrapping
  - Produces figure 3 from the manuscript




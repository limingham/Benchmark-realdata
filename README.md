
<img width="4798" height="1703" alt="Figure2" src="https://github.com/user-attachments/assets/bb9b3c33-5960-4f8a-980d-e264e3f480d2" />


# Evaluating deconvolution methods using real bulk RNA-expression data for robust prognostic insights across cancer types

### Background:
Deconvolution of bulk RNA-expression data unlocks the cellular complexity of cancer, yet traditional pseudobulk benchmarks may not always be reliable in real-world settings where absolute cell proportions are unknown.

### Results:
Here, we introduce a novel real-data framework, leveraging 18 real bulk RNA-expression cohorts (5,891 samples) across nine cancer types to evaluate five deconvolution methods based on differentially proportioned (DP) and prognosis-related (PR) cell types. Across three innovative benchmark scenarios—consistency with scRNA-seq, reproducibility across cohorts, and reproducibility of prognostic relevance—ReCIDE and BayesPrism stand out as two robust deconvolution methods. Application of a pan-cancer analysis based on the deconvolution of TCGA cohorts identifies matrix cancer-associated fibroblasts (mCAF) as a prognostic marker with consistent effects across multiple cancers. Building on this finding, we find a prognostic indicator combining classical monocytes and mCAF cell proportions to be significant in five TCGA cohorts, which we further validate in five independent GEO cohorts.

### Conclusions:
This study broadens deconvolution benchmarking, offering actionable tools for precision oncology and guiding method selection for translational research.

                
<img width="5224" height="5811" alt="Figure2 - BCDE" src="https://github.com/user-attachments/assets/10f9c5fb-51da-4283-94de-4df63dcc918a" />
            
### Publicly accessible code and results
All deconvolution results and benchmark codes have now been released!  All deconvolution results are located in the 'benchmark_deconvolution_results' folder, and the code and evaluation results corresponding to each benchmark scenario can be found in the respective figure folder. 

### Data Used in the Study
Our research utilized extensive single-cell RNA-seq and bulk RNA-expression data, with their download links provided in the Public_data folder.

### Contact us
limh25@m.fudan.edu.cn



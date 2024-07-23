# Regime-switches-and-Bayesian-wavelet-estimation
This repository contains the code to reproduce the numerical experiments presented in the paper "Identifying regime switches through Bayesian wavelet estimation: application to environmental and genetic data." This study uses a dynamic two-component mixture model, where the mixture weights vary with the data index, enhancing adaptability to diverse datasets. An efficient Gibbs sampling algorithm jointly estimates the mixture component parameters and dynamic weights. The method employs a wavelet-based data augmentation technique, leveraging wavelet bases' properties in curve estimation. Validation through Monte Carlo simulations and real-world applications, including flood detection in a Brazilian river and chromosomal anomaly identification in DNA samples, demonstrates the robustness and versatility of the proposed approach.

# Key Findings

- **Dynamic Mixture Weights:** Traditional models assume constant mixture weights, but our approach allows these weights to vary dynamically, making the model more adaptable to different datasets.
- **Efficient Gibbs Sampling Algorithm:** An efficient Gibbs sampling algorithm was introduced to jointly estimate the mixture component parameters and dynamic weights.
- **Wavelet-based Data Augmentation:** We proposed a wavelet-based version of data augmentation to assess the dynamic behavior of the mixture weights, exploiting the advantageous properties of wavelet bases in curve estimation.
- **Monte Carlo Simulations:** The method was validated through extensive Monte Carlo simulations, demonstrating its robustness and accuracy in a controlled setting.
- **Real-world Applications:** The method was applied to two real-world datasets:
  - **Flood Detection in Southern Brazil:** Time series data of river water levels were analyzed to distinguish between inundation and non-inundation periods.
  - **Chromosomal Anomalies in DNA Samples:** Genetic data (aCGH data) were analyzed to identify regions of genomic imbalances in test DNA samples.
- **Robust Results:** The results aligned with existing findings, showcasing the robustness and effectiveness of our method in different applications.

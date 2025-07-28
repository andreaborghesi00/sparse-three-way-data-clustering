# Sparse Model-Based Clustering of Three-Way Environmental Data

This project explores model-based clustering of three-way data using sparse, penalized matrix-variate Gaussian mixture models. The goal is to extract interpretable pollution patterns across time, sensors, and pollutants from the Madrid air quality dataset.

> This project was developed as part of a graduate-level statistics seminar. All modeling was implemented in **R**, using specialized libraries for multivariate analysis and sparse estimation.

---

## Objective

Given a 3D tensor of air quality measurements:
- **Dimensions**: `[Time × Stations × Pollutants]`
- Perform **unsupervised clustering** of the time dimension, i.e., identify recurring patterns in pollutant concentrations across the city.
- Encourage **sparse estimation of covariance matrices** to improve interpretability and reduce overfitting.

---

## Methods

### Data Preparation
- Raw data collected from multiple air quality stations in Madrid.
- Preprocessing included:
  - Filtering incomplete sensors
  - Reshaping into a time-series tensor
  - Scaling and centering

### Clustering Model
- **Matrix-variate Gaussian Mixture Model**:
  - Captures structured covariance across spatial and variable dimensions.
  - Supports parameter sharing between mixture components.
- **Sparsity-inducing penalties**:
  - Graphical lasso applied to covariance matrices.
  - Penalized log-likelihood used for model selection.

### Model Selection
- Explored a wide range of model configurations (240 total).
- Evaluated via BIC and convergence diagnostics.
- Handled failures due to ill-conditioning, non-convergence, and degeneracy.

---

## Results

- Optimal model identified **4 temporal clusters** reflecting distinct air pollution regimes.
- Sparsity in covariance structures enhanced interpretability.
- Visualizations reveal strong patterns between NOx, O3, and traffic-linked pollutants.
- Weak correlation between cluster assignments and external meteorological variables (e.g., temperature, humidity), suggesting pollution patterns are not trivially explained by weather alone.

---

## Tech Stack

- **Language**: R
- **Key Libraries**: `mixvar`, `glasso`, `ggplot2`, `reshape2`

---

### References
- Friedman, J., Hastie, T., & Tibshirani, R. (2008). *Sparse inverse covariance estimation with the graphical lasso*.



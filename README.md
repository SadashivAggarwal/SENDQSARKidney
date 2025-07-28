# ![SENDQSARKidney Logo](https://img.icons8.com/fluency/48/database.png) SENDQSARKidney

## An R Package for QSAR Modeling with SEND Database

---

![QSAR Modeling](https://img.icons8.com/fluency/96/000000/artificial-intelligence.png)

### Overview

**SENDQSARKidney** is an R package designed to accelerate **Quantitative Structure-Activity Relationship (QSAR)** modeling using data from the Standard for Exchange of Nonclinical Data (**SEND**) database.  
It streamlines:

- Data acquisition
- Pre-processing
- Organ-specific toxicity score calculation
- Molecular descriptor calculation
- Machine learning model development and evaluation

Comprehensive function documentation is available in the **Articles** section on the package’s GitHub site.

---

## Features

- 🚀 **Automated Data Processing**: Fast, efficient workflows for pre-processing SEND data
- 📊 **Comprehensive Analysis**: Calculate z-scores for:
  - Body weight
  - Kidney-to-body weight ratio
  - Laboratory test results
- 🤖 **Machine Learning Integration**: 
  - Classification model training
  - Hyperparameter tuning
  - Model performance assessment
- 📈 **Visualization Tools**:
  - Histograms, bar plots, AUC curves, and more for insightful data exploration

---

## Workflow

1. **Input Database Path**  
   Supply the path to your database or `.xpt` files containing nonclinical SEND-formatted study data.  
   ![Database Icon](https://img.icons8.com/color/48/database.png)
2. **Data Pre-processing**  
   Use modular functions (`f1` to `f8`) to clean, harmonize, and prepare your data for machine learning applications.  
   ![Processing Icon](https://img.icons8.com/office/40/process.png)

---

## Modular Functions

### 1. Data Extraction & Score Calculation

- **f1: `get_compile_data`**  
  Extracts structured data from the specified database path.
- **f2: `get_bw_score`**  
  Calculates body weight z-scores for each animal (calls `f1`).
- **f3: `Kidneytobw_zscore`**  
  Computes kidney-to-body weight z-scores (calls `f1`).
- **f4: `Kidney_lb_score`**  
  Calculates z-scores for laboratory test results (calls `f1`).
- **f5: `Kidney_mi_score`**  
  Computes z-scores for microscopic findings (calls `f1`).

### 2. Aggregation & Machine Learning Preparation

- **f6: `process_kidney_score`**  
  Aggregates LB, MI, and kidney-to-BW ratio z-scores into a unified data frame (uses `f1`–`f5`).
- **f7: `get_col_harmonized_scores_df`**  
  Harmonizes and standardizes column names for downstream consistency (uses `f6`).
- **f8: `get_ml_data_and_tuned_hyperparameters`**  
  Prepares machine learning-ready data and tunes model hyperparameters (uses `f7`).

### 3. Model Building & Evaluation

- **f9: `get_rf_model_with_cv`**  
  Builds a random forest model with cross-validation (uses `f8`).
- **f10: `get_zone_exclusioned_rf_model_with_cv`**  
  Enhances classification accuracy by excluding uncertain predictions (uses `f8`).

---

## Helper Functions

Additional modular helpers support each analysis step, ensuring flexible, readable, and robust workflows.

---

## In Development

- New modules for additional organ endpoints
- Expanded ML workflows
- Advanced visualization and reporting tools

---

## Dependencies

- **R** (≥ 4.1.0)
- Key packages:  
  `dplyr`, `tidyr`, `data.table`, `ggplot2`, `caret`, `xgboost`, `readxl`, `DBI`, `ROCR`  
  *(See the DESCRIPTION file for the complete list.)*

---

## Installation

```r
# Install from GitHub
remotes::install_github("yourusername/SENDQSARKidney")

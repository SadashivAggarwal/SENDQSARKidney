# ![SENDQSARKidney Logo](https://img.icons8.com/fluency/48/database.png) SENDQSARKidney

## SENDQSAR: An R Package for QSAR Modeling with SEND Database

---

![QSAR Modeling](https://img.icons8.com/fluency/96/000000/artificial-intelligence.png)

### About

- **SENDQSARKidney** is an R package to facilitate **Quantitative Structure-Activity Relationship (QSAR)** modeling using data from the Standard for Exchange of Nonclinical Data (**SEND**) database.
- It streamlines:
  - Data acquisition
  - Pre-processing
  - Organ-wise toxicity score calculation
  - Descriptor calculation
  - Machine learning-based model evaluation
- All functions are described in detail in the **“Articles”** section on the package’s GitHub website.

---

## Features

- 🚀 **Automated Data Processing:** Simple and efficient pre-processing workflows.
- 📊 **Comprehensive Analysis:** Z-score calculations for:
  - Body weight
  - Liver-to-body weight ratio
  - Laboratory tests
- 🤖 **Machine Learning Integration:** 
  - Supports classification models, hyperparameter tuning, and performance evaluation.
- 📈 **Visualization Tools:**
  - Histograms, bar plots, AUC curves, and more for insightful data interpretation.

---

## Workflow

1. **Input Database Path:**  
   Provide the path to the database or `.xpt` files with nonclinical study data in SEND format.  
   ![Database Icon](https://img.icons8.com/color/48/database.png)
2. **Data Pre-processing:**  
   Use modular functions (`f1` to `f8`) to clean, harmonize, and prepare your data for machine learning applications.  
   ![Processing Icon](https://img.icons8.com/office/40/process.png)

---

## Modular Functions Overview

#### Liver Toxicity Score Calculation for Individual STUDYID

- **f1: `get_compile_data`**  
  - Fetches structured data from the specified database path.
- **f2: `get_bw_score`**  
  - Calculates body weight z-scores for each animal (uses `f1`).
- **f3: `Kidneytobw_zscore`**  
  - Computes liver-to-body weight z-scores (uses `f1`).
- **f4: `Kidney_lb_score`**  
  - Calculates z-scores for laboratory test results (uses `f1`).
- **f5: `Kidney_mi_score`**  
  - Computes z-scores for microscopic findings (uses `f1`).

#### Liver Toxicity Score Calculation and Aggregation for Multiple STUDYID

- **f6: `process_kidney_score`**
  - Aggregates z-scores for LB, MI, and liver-to-BW ratio into a single data frame.
  - Internally calls `f1`–`f5`.

---

## Helper Functions

- Modular helper functions to support each analysis step, ensuring code clarity and flexibility.

---

## Functions in Development

- Additional modules are planned for:
  - Organ-specific endpoints
  - Expanded ML workflows
  - Enhanced visualization and reporting

---

## Dependencies

- **R (≥ 4.1.0)**
- **Key packages:** `dplyr`, `tidyr`, `data.table`, `ggplot2`, `caret`, `xgboost`, `readxl`, `DBI`, `ROCR`
- See `DESCRIPTION` for a full list

---

## Installation

```R
# Install from GitHub
remotes::install_github("yourusername/SENDQSARKidney")

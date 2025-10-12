# 🧬 Bacteriocin Prediction Through Cross-Validation-Based and Hypergraph-Based Feature Evaluation Approaches

This repository provides the full workflow, datasets, and scripts used in the study:  
**“Bacteriocin Prediction Through Cross-Validation-Based and Hypergraph-Based Feature Evaluation Approaches.”**  
The project introduces an integrated framework combining **feature extraction**, **Pearson correlation filtering**, and **feature evaluation** using **cross-validation-based (CVFS)** and **hypergraph-based (HFE)** approaches to enhance bacteriocin prediction performance.

---

## 🧩 Workflow Overview

### 1. **Data Collection**
Two sequence sets were curated:
- **Bacteriocin sequences** (`positiveBac.faa`)
- **Non-bacteriocin sequences** (`negativeBac.faa`)

---

### 2. **Feature Extraction**
Physicochemical and structural descriptors were extracted from amino acid sequences using **feature_extract.R**.  
Feature categories and dimensions:

| Feature Type | Dimension |
|---------------|------------|
| Amino acid composition | 20D |
| Dipeptide composition | 400D |
| Pseudo amino acid composition | 30D |
| Amphiphilic pseudo amino acid composition | 40D |
| Composition (C) | 21D |
| Transition (T) | 21D |
| Distribution (D) | 105D |
| Secondary structure | 6D |
| Sequence-order-coupling number | 20D |
| Quasi-sequence-order descriptors | 40D |
| Position-specific scoring matrix-based bigrams | 400D |

---

### 3. **Feature Filtering (Pearson Correlation Coefficient)**
Highly correlated features (r ≥ 0.90) were removed using **pearsons_corr.R**.  
**Output:** `pearsoncorrelation.csv`

---

### 4. **Feature Evaluation Approaches**

#### 🧮 Cross-Validation-Based Feature Selection (CVFS)
Implements cross-validation-driven ranking and selection of the most predictive features.  
- Files: `cvfs.csv`, `cvfs.py`  
- Run using:
  ```bash
  cvfs.py -i cvfs.csv -o C2P5P0.4 -c 2 -e 5 -p 0.4
  ```

#### 🧠 Hypergraph-Based Feature Evaluation (HFE)
Employs a hypergraph representation of feature interrelationships to evaluate multi-feature dependencies.  
- Notebook: `hfe.ipynb` (executed in Jupyter)

---

### 5. **Machine Learning Model Implementation**
Final model training and interpretation were conducted using **XGBoost**.

- **Script:** `shap.R`  
- Performs:
  - Training and validation on selected features  
  - SHAP-based feature importance analysis  
  - Grouped SHAP aggregation (CTD: Composition, Transition, Distribution)  
  - Model evaluation via ROC–AUC, MCC, and confusion matrix

---

## 📊 Outputs

| File | Description |
|------|--------------|
| `pearsoncorrelation.csv` | Filtered feature set after correlation thresholding |
| `selected_training_merged_file.csv` | Processed training dataset |
| `selected_validation_merged_file.csv` | Processed validation dataset |

---

## 🌐 Web Application

Explore and run predictions using the interactive web app:  
🔗 **[Bacteriocin Prediction Web Application](https://shiny.tricities.wsu.edu/bacteriocin-prediction/)**

If you find our web application useful, please consider citing the following works:

> **Akhter, S. and Miller, J.H.**  
> *BPAGS: A web application for bacteriocin prediction via feature evaluation using alternating decision tree, genetic algorithm, and linear support vector classifier.*  
> *Frontiers in Bioinformatics*, 3, p.1284705.

> **Akhter, S., & Miller, J.H.**  
> *BaPreS: a software tool for predicting bacteriocins using an optimal set of features.*  
> *BMC Bioinformatics*, 24(1), 313.

---

## ⚙️ Requirements

### **R (≥ 4.0)**
Required packages:  
`seqinr`, `protr`, `Biostrings`, `xgboost`, `SHAPforxgboost`, `caret`, `pROC`, `ggplot2`

### **Python (≥ 3.8)**
Required libraries:  
`pandas`, `numpy`, `scikit-learn`, `networkx`

---

## 📂 Repository Structure

```
├── feature_extract.R
├── pearsons_corr.R
├── shap.R
├── cvfs.py
├── cvfs.csv
├── hfe.ipynb
├── positiveBac.faa
├── negativeBac.faa
├── pearsoncorrelation.csv
├── selected_training_merged_file.csv
└── selected_validation_merged_file.csv
```

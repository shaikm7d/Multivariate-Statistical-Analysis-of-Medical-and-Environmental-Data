# 🧮 Multivariate Statistical Analysis in R  
**Author:** Shaik Mohammed

## 📘 Overview  
This project performs a comprehensive **multivariate statistical analysis** in R using two datasets:  
1. `data_treatment.csv` - a medical dataset for analyzing treatment effects on cholesterol and glucose.  
2. `data_pollution.csv` - an environmental dataset for clustering cities by pollution and meteorological characteristics.

It applies **Hotelling’s T² test**, **MANOVA**, and **K-Means clustering** after systematic data cleaning and visualization, providing both inferential and exploratory insights.

---

## 🎯 Problem Statement  

### **Part 1 - Medical Dataset (`data_treatment.csv`)**
Determine whether:
- The population mean vector of two health indicators (Total Cholesterol, Blood Glucose) equals a specified target vector [5, 5] claimed by a researcher.
- There is a multivariate difference in these indicators across three treatment groups (T1, T2, T3).

### **Part 2 - Pollution Dataset (`data_pollution.csv`)**
Identify meaningful clusters of cities based on environmental variables (SO₂, temperature, precipitation, population, etc.) using correlation analysis and K-means clustering.

---

## 🧠 Methodology  

### **Part 1:**
### **Hotelling’s One-Sample T² Test**
**Objective:** Test whether the true mean vector equals [5, 5].  
**Steps:**
1. Compute sample mean and covariance matrix.  
2. Calculate `T² = n (X̄ − μ)' S⁻¹ (X̄ − μ)`  
3. Convert T² to an approximate F-statistic:  
   `F = ((n - p) / (p * (n - 1))) * T²`  
   where `p` = number of variables, `n` = sample size.  
4. Obtain p-value from the F-distribution.

**Result Summary:**  
- `T² = 76.9971`, `F = 37.171`, `df₁ = 2`, `df₂ = 28`  
- `p = 1.317×10⁻⁸`  
✅ **Reject H₀** → The true mean vector ≠ [5, 5].  

---

### **Boxplot Visualization**
![Cholesterol & Glucose by Treatment](./1.(iii)_boxplots.png)

Visualized distributions of:
- **V1 (Total Cholesterol)**  
- **V2 (Blood Glucose)**  

**Insights:**
- Treatment 2 (T2) shows highest medians for both variables.
- Treatment 1 (T1) yields the lowest readings.
- Treatment 3 (T3) lies in between, with moderate variability.

---

### **One-Way MANOVA**
**Goal:** Test whether mean vectors of (V1, V2) differ across treatments.  
**Test statistic:** Wilks’ Lambda `Λ = |W| / |W + B|`

**Computation results:**  
- `W = [[32.741, 4.592], [4.592, 22.569]]`  
- `B = [[37.948, 13.784], [13.784, 19.465]]`  
- `Λ = 0.2725`  
- `F = 11.901`, `df₁ = 4`, `df₂ = 52`  
- `p = 6.15×10⁻⁷`

✅ **Reject H₀** → There is a significant multivariate difference between treatment means.  
**Conclusion:** The treatments significantly affect cholesterol and glucose jointly.

---

### **Part 2:**
### **Pollution Dataset Analysis**

#### **Data Preprocessing**
- Outlier capping using **IQR method**: replaced values beyond `Q1 − 1.5 × IQR` or `Q3 + 1.5 × IQR`.  
- Imputed missing values with **median** (robust to outliers).  
- Standardized all numeric variables (`mean = 0`, `SD = 1`) for fair distance-based clustering.

#### **Exploratory Correlation**
- Computed **Pearson correlation matrix** and visualized via heatmap.
![Correlation Matrix Heatmap](./2.(ii)_corrplot.png)
- Observations:  
  - Strong +ve correlation: Manufacturing (`y₃`) ↔ Population (`y₄`) = 0.83.  
  - Moderate +ve: SO₂ ↔ Precipitation days (`y₇`) = 0.49.  
  - Moderate −ve: Temperature (`y₂`) ↔ SO₂ (`y₁`) = −0.42.  

---

### **K-Means Clustering**
- Used **Elbow method** → optimal `k = 3`.
![Elbow Method for Choosing k](./2.(ii)_elbow.png)
- Performed PCA for 2D visualization of clusters.
![K-means Clusters (PCA Projection)](./2.(ii)_pca_clusters.png)
- Evaluated mean pollution and climate features per cluster.

#### **Cluster Interpretation**
| Cluster | Description | Key Traits |
|:--:|:--|:--|
| **1** | Urban-industrial hubs | High SO₂ (45.9), High population (967k), High manufacturing (814). |
| **2** | Medium, wetter cities | Moderate pollution, highest precipitation (43.9 in), warmer temps. |
| **3** | Cleaner, drier cities | Lowest SO₂ (16.0), moderate population, low rainfall. |

**Conclusion:** Industrialization and population density strongly influence pollution levels. Meteorological variables (wind, rain, temperature) modulate these effects.

---

## ⚙️ Technologies & Packages  
- **R base packages:** `stats`, `graphics`, `grDevices`  
- **Libraries:** `MASS`, `corrplot`, `ggplot2`  
- **Methods Used:**  
  - Hotelling’s T² Test  
  - MANOVA (Wilks’ Lambda)  
  - K-Means Clustering  
  - Principal Component Analysis  

---

## 📈 Key Findings
| Analysis | Result Summary | Inference |
|-----------|----------------|------------|
| **Hotelling’s T² Test** | `p ≈ 1.3e−8` | Significant deviation from [5, 5] mean vector |
| **MANOVA** | `Λ = 0.2725`, `F = 11.90`, `p ≈ 6.15e−7` | Treatment effect significant |
| **Clustering** | `k = 3` optimal; clusters clearly separated | Industrialization drives SO₂ variation |

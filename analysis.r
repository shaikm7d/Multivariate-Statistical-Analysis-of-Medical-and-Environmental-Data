# Read the dataset

treatment_data <- read.csv("data_treatment.csv", header = TRUE, stringsAsFactors = FALSE)
pollution_data <- read.csv("data_pollution.csv", header = TRUE, stringsAsFactors = FALSE)

# Create a dedicated folder to save plots
plots_dir <- "plots"
if (!dir.exists(plots_dir)) dir.create(plots_dir)

# 1.(i) Hotelling’s T² Test: One Population

if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
library(MASS)

hotelling_test_one_population <- function(data, mu) {
  if (!is.matrix(data)) stop("Data must be a matrix.")
  if (ncol(data) != length(mu)) stop("Mean vector length must match number of columns in data.")
  
  n <- nrow(data)
  p <- ncol(data)
  x_bar <- colMeans(data)
  S <- cov(data)
  S_inv <- solve(S)
  
  diff <- x_bar - mu
  T2 <- n * t(diff) %*% S_inv %*% diff
  
  F_stat <- ((n - p) / (p * (n - 1))) * T2
  df1 <- p
  df2 <- n - p
  p_value <- 1 - pf(F_stat, df1, df2)
  
  return(list(
    T2_stat = as.numeric(T2),
    p_value = as.numeric(p_value),
    F_stat = as.numeric(F_stat),
    df1 = df1,
    df2 = df2
  ))
}

# 1.(ii) Run Hotelling’s Test

data_matrix <- as.matrix(treatment_data[, c("V1", "V2")])
mu0 <- c(5, 5)
result <- hotelling_test_one_population(data_matrix, mu0)

cat("Hotelling's T² statistic:", round(result$T2_stat, 4), "\n")
cat("p-value:", signif(result$p_value, 4), "\n")
cat("F-statistic:", round(result$F_stat, 4), 
    "with df1 =", result$df1, "and df2 =", result$df2, "\n")

# 1.(iii) Boxplots by Treatment

png(file.path(plots_dir, "1.(iii)_boxplots.png"), width = 1000, height = 450)
par(mfrow = c(1, 2))

boxplot(V1 ~ Treatment, data = treatment_data,
        main = "V1: Total Cholesterol",
        xlab = "Treatment",
        ylab = "mmol/L",
        col = "lightblue", border = "darkblue", 
        outline = TRUE)

boxplot(V2 ~ Treatment, data = treatment_data,
        main = "V2: Blood Glucose",
        xlab = "Treatment",
        ylab = "mmol/L",
        col = "lightgreen", border = "darkgreen", 
        outline = TRUE)

par(mfrow = c(1, 1))
dev.off()

# 1.(iv) MANOVA Components

p  <- 2
g  <- length(unique(treatment_data$Treatment))
N  <- nrow(treatment_data)

y_bar <- colMeans(treatment_data[, 1:2])
print(y_bar)

W <- matrix(0, nrow = p, ncol = p)
B <- matrix(0, nrow = p, ncol = p)

for(tr in unique(treatment_data$Treatment)) {
  grp     <- treatment_data[treatment_data$Treatment == tr, 1:2]
  n_i     <- nrow(grp)
  y_bar_i <- colMeans(grp)
  
  for(j in 1:n_i) {
    diff <- matrix(as.numeric(grp[j, ]) - y_bar_i, ncol = 1)
    W <- W + diff %*% t(diff)
  }
  
  diff_bar <- matrix(y_bar_i - y_bar, ncol = 1)
  B <- B + n_i * diff_bar %*% t(diff_bar)
}

print("Within-group sum of squares and cross-products (W):")
print(W)
print("Between-group sum of squares and cross-products (B):")
print(B)

wilklambda  <- det(W) / det(W + B)
F_stat  <- ((N - g - 1)/2) * ((1 - sqrt(wilklambda)) / sqrt(wilklambda))
df1     <- p * (g - 1)
df2     <- 2 * (N - 4)
p_val   <- pf(F_stat, df1, df2, lower.tail = FALSE)

wilklambda; F_stat; df1; df2; p_val

manova_tab <- list(
  Treatment = list(SSP = B, df = g - 1),
  Residual  = list(SSP = W, df = N - g),
  Total     = list(SSP = B + W, df = N - 1)
)
print(manova_tab)

manova_model <- manova(cbind(V1, V2) ~ Treatment, data = treatment_data)
summary_manova <- summary(manova_model, test = "Wilks")
print(summary_manova)

#  2.(i) Outlier Capping, Cleaning, and Scaling

rownames(pollution_data) <- pollution_data[, 1]
pollution_data <- pollution_data[, -1]

cap_outliers_dataset <- function(data, numeric_cols, factor = 1.5) {
  detect_outliers_iqr <- function(x) {
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1
    lower <- Q1 - factor * IQR_val
    upper <- Q3 + factor * IQR_val
    x < lower | x > upper
  }
  
  cap_outliers_iqr <- function(x) {
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1
    lower <- Q1 - factor * IQR_val
    upper <- Q3 + factor * IQR_val
    x[x < lower] <- lower
    x[x > upper] <- upper
    x
  }
  
  outlier_flags <- as.data.frame(
    sapply(numeric_cols, function(col) detect_outliers_iqr(data[[col]]))
  )
  outlier_counts <- colSums(outlier_flags)
  cat("Number of outliers in each variable:\n")
  print(outlier_counts)
  
  data_capped <- data
  for (col in numeric_cols) {
    data_capped[[col]] <- cap_outliers_iqr(data_capped[[col]])
  }
  
  return(data_capped)
}

numeric_cols <- c("y_1", "y_2", "y_3", "y_4", "y_5", "y_6", "y_7")
pollution_data_capped <- cap_outliers_dataset(pollution_data, numeric_cols)

colSums(is.na(pollution_data_capped))

for (v in numeric_cols) {
  pollution_data_capped[[v]][is.na(pollution_data_capped[[v]])] <- 
    median(pollution_data_capped[[v]], na.rm = TRUE)
}

pollution_data_scaled <- pollution_data_capped
for (col in numeric_cols) {
  mean_val <- mean(pollution_data_capped[[col]], na.rm = TRUE)
  sd_val <- sd(pollution_data_capped[[col]], na.rm = TRUE)
  pollution_data_scaled[[col]] <- (pollution_data_capped[[col]] - mean_val) / sd_val
}

# 2.(ii) Correlation and Clustering

summary(pollution_data_scaled)
cor_matrix <- cor(pollution_data_scaled, method = "pearson")

png(file.path(plots_dir, "2.(ii)_corrplot.png"), width = 800, height = 700)
library(corrplot)
corrplot(cor_matrix, 
         method = "color",
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black",
         col = colorRampPalette(c("#6BAED6", "#F7F7F7", "#FB6A4A"))(200)
)
dev.off()

set.seed(123)
max_k <- 10
wss <- numeric(max_k)
for (k in 1:max_k) {
  km_model <- kmeans(pollution_data_scaled, centers = k, nstart = 100)
  wss[k] <- km_model$tot.withinss
}

png(file.path(plots_dir, "2.(ii)_elbow.png"), width = 800, height = 500)
plot(1:max_k, wss, type="b",
     xlab="Number of clusters (k)",
     ylab="Total within-cluster sum of squares (SSE)",
     main="Elbow Method for Choosing Optimal k")
dev.off()

optimal_k <- 3
set.seed(123)
final_kmeans <- kmeans(pollution_data_scaled, centers = optimal_k, nstart = 100)

cluster_assignments <- final_kmeans$cluster
print(cluster_assignments)
table(cluster_assignments)

library(ggplot2)
pca_result <- prcomp(pollution_data_scaled)
pca_df <- data.frame(pca_result$x[, 1:2])
pca_df$Cluster <- as.factor(cluster_assignments)
pca_df$City <- rownames(pollution_data_scaled)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = City)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text(vjust = -0.5, size = 3, check_overlap = TRUE) +
  labs(title = "K-means Clusters with City Labels (PCA Projection)",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey70", fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  scale_color_brewer(palette = "Set2")

ggplot2::ggsave(
  filename = file.path(plots_dir, "2.(ii)_pca_clusters.png"),
  plot = p,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

aggregate(pollution_data_capped[, 1:7], by = list(final_kmeans$cluster), FUN = mean)

##########################################
#             Load Libraries             #
##########################################
library(dplyr)
library(tidyverse)
library(cluster)
library(mclust)
library(mixtools)
library(poLCA)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(softImpute)
library(patchwork)
library(clustMixType)
library(cowplot)


##########################################
#          Load and Preprocess Data      #
##########################################
# Load data
Dash_Data <- read.csv("C:/Users/loisr/OneDrive/Documents/JKim_Project/merged_df.csv") %>%
  arrange(rand_id)

# Filter for treatment group
Dash_Data_treat <- Dash_Data %>%
  filter(treatmnt == 2)

# Remove columns with >20% missing
missing_fraction <- colMeans(is.na(Dash_Data_treat))
columns_to_keep <- names(missing_fraction[missing_fraction <= 0.2])
Dash_Data_filtered <- Dash_Data_treat[, columns_to_keep]


##########################################
#     Handle Data Types Separately       #
##########################################
# Convert character to factor
Dash_Data_filtered <- Dash_Data_filtered %>%
  mutate(across(where(is.character), as.factor))

# Identify types
numeric_cols <- names(select(Dash_Data_filtered, where(is.numeric)))
factor_cols <- names(select(Dash_Data_filtered, where(is.factor)))

# Convert binary numeric to factor
binary_numeric_cols <- sapply(Dash_Data_filtered[numeric_cols], function(x) length(unique(na.omit(x))) == 2)
Dash_Data_filtered[numeric_cols[binary_numeric_cols]] <- lapply(Dash_Data_filtered[numeric_cols[binary_numeric_cols]], as.factor)

# Remove near-zero variance columns
Dash_Data_filtered <- Dash_Data_filtered[, sapply(Dash_Data_filtered, function(x) length(unique(x)) > 1)]


##########################################
#     Imputation for Missing Values      #
##########################################
# 1. Numeric → Median Imputation
num_cols <- names(select(Dash_Data_filtered, where(is.numeric)))
Dash_Data_filtered[num_cols] <- lapply(Dash_Data_filtered[num_cols], function(x) {
  x[is.nan(x)] <- NA
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(scale(x))
})

# 2. Factor → Mode Imputation
get_mode <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}
cat_cols <- names(select(Dash_Data_filtered, where(is.factor)))
Dash_Data_filtered[cat_cols] <- lapply(Dash_Data_filtered[cat_cols], function(x) {
  x[is.na(x)] <- get_mode(x)
  return(x)
})


##########################################
#     Gower Distance + PAM Clustering    #
##########################################
gower_dist <- daisy(Dash_Data_filtered, metric = "gower")

# Optimal k using silhouette
silhouette_scores <- sapply(2:10, function(k) {
  pam_fit <- pam(gower_dist, k = k)
  mean(pam_fit$silinfo$avg.width)
})
optimal_k <- which.max(silhouette_scores) + 1

# PAM clustering
pam_fit <- pam(gower_dist, k = optimal_k)
Dash_Data_filtered$pam_cluster <- factor(pam_fit$clustering)

# Visualize PAM
mds_fit <- cmdscale(gower_dist, k = 2)
mds_df <- data.frame(Dim1 = mds_fit[,1], Dim2 = mds_fit[,2], Cluster = Dash_Data_filtered$pam_cluster)

pam_plot <- ggplot(mds_df, aes(Dim1, Dim2, color = Cluster)) +
  geom_point(size = 2) +
  stat_ellipse(aes(fill = Cluster), geom = "polygon", alpha = 0.2) +
  labs(title = "PAM Clustering (MDS Visualization)") +
  theme_classic()


##########################################
#       K-Prototypes Clustering          #
##########################################
set.seed(123)
kp_fit <- kproto(Dash_Data_filtered, k = optimal_k)
Dash_Data_filtered$kproto_cluster <- factor(kp_fit$cluster)

mds_df_kproto <- mds_df
mds_df_kproto$Cluster <- Dash_Data_filtered$kproto_cluster

kproto_plot <- ggplot(mds_df_kproto, aes(Dim1, Dim2, color = Cluster)) +
  geom_point(size = 2) +
  stat_ellipse(aes(fill = Cluster), geom = "polygon", alpha = 0.2) +
  labs(title = "K-Prototypes Clustering (MDS Visualization)") +
  theme_classic()


##########################################
#      Model-Based Clustering (mclust)   #
##########################################
numeric_data <- Dash_Data_filtered %>% select(where(is.numeric))
mclust_fit <- Mclust(numeric_data)
Dash_Data_filtered$mclust_cluster <- factor(mclust_fit$classification)

pca_mclust <- prcomp(numeric_data)
pca_df <- data.frame(PC1 = pca_mclust$x[,1], PC2 = pca_mclust$x[,2],
                     Cluster = Dash_Data_filtered$mclust_cluster)

mclust_plot <- ggplot(pca_df, aes(PC1, PC2, color = Cluster)) +
  geom_point(size = 2) +
  labs(title = "Model-Based Clustering (PCA Visualization)") +
  theme_minimal()


##########################################
#      Latent Class Analysis (LCA)       #
##########################################
cat_data <- Dash_Data_filtered %>% select(where(is.factor))

f <- as.formula(paste("cbind(", paste(names(cat_data), collapse = ", "), ") ~ 1"))
lca_fit <- poLCA(f, data = cat_data, nclass = 3, maxiter = 1000, na.rm = FALSE)
Dash_Data_filtered$Cluster_LCA <- factor(lca_fit$predclass)

mca_res <- MCA(cat_data, graph = FALSE)
mca_df <- as.data.frame(mca_res$ind$coord[, 1:2])
mca_df$Cluster <- Dash_Data_filtered$Cluster_LCA

lca_plot <- ggplot(mca_df, aes(Dim.1, Dim.2, color = Cluster)) +
  geom_point(size = 2) +
  labs(title = "Latent Class Analysis (MCA Visualization)") +
  theme_minimal()


##########################################
# SoftImpute + KMeans Clustering (Num)   #
##########################################
num_matrix <- as.matrix(numeric_data)
fit_si <- softImpute(num_matrix, rank.max = 50, lambda = 10)
completed_matrix <- complete(num_matrix, fit_si)

kmeans_fit <- kmeans(completed_matrix, centers = 3)
Dash_Data_filtered$Cluster_SoftImpute <- factor(kmeans_fit$cluster)

pca_res <- prcomp(completed_matrix)
pca_si <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Cluster = Dash_Data_filtered$Cluster_SoftImpute)

softimpute_plot <- ggplot(pca_si, aes(PC1, PC2, color = Cluster)) +
  geom_point(size = 2) +
  labs(title = "SoftImpute + KMeans Clustering (PCA)") +
  theme_minimal()


##########################################
#         Combine and Save Plots         #
##########################################
combined_plot <- pam_plot / kproto_plot / mclust_plot / lca_plot / softimpute_plot +
  plot_layout(ncol = 1)

ggsave("C:/Users/loisr/OneDrive/Documents/Clustering_Results/All_Clusters.pdf",
       plot = combined_plot, width = 10, height = 20)

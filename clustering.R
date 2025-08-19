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
library(stats)
library(patchwork)
library(clustMixType)
library(cowplot)



############################# Loading data + data processing #############################
# Reading csv files
Dash_Data <- read.csv(file = "C:/Users/loisr/OneDrive/Documents/JKim_Project/merged_df.csv")
Dash_Data <- Dash_Data %>% 
  dplyr::arrange(rand_id)
View(Dash_Data)


# Dash_Data_4 <- read.csv(file = "C:/Users/loisr/OneDrive/Documents/JKim_Project/not_filtered_dataframes/dataframe_4.csv")
# # Removing days 3 and 5
# Dash_Data_4 <- Dash_Data_4[Dash_Data_4$HC4DAY%in%c(1,7,14), ]
# # Sorting in ascending order by rand_id
# Dash_Data_4 <- Dash_Data_4 %>%
#   dplyr::arrange(HC4DAY, rand_id)
# View(Dash_Data_4)


# Filtering for 'treatment' == 2 
Dash_Data_treat <- Dash_Data %>%
  dplyr::filter(treatmnt == 2)
#View(Dash_Data_treat)


# Filtering features with 20% missing data
missing_fraction <- colMeans(is.na(Dash_Data_treat))
columns_to_keep <- names(missing_fraction[missing_fraction <= 0.2])
Dash_Data_filtered <- Dash_Data_treat[, columns_to_keep]
# View(Dash_Data_filtered)


# Apply scaling and transformations
Dash_Data_filtered[setdiff(names(Dash_Data_filtered)[sapply(Dash_Data_filtered, is.numeric)], "treatmnt")] <- 
  scale(Dash_Data_filtered[setdiff(names(Dash_Data_filtered)[sapply(Dash_Data_filtered, is.numeric)], "treatmnt")])

Dash_Data_filtered <- Dash_Data_filtered %>%
  mutate(across(where(is.character), as.factor))

#str(Dash_Data_filtered)

numeric_cols <- setdiff(names(Dash_Data_filtered)[sapply(Dash_Data_filtered, is.numeric)], "treatmnt")

# Median data imputation 
Dash_Data_filtered[numeric_cols] <- lapply(Dash_Data_filtered[numeric_cols], function(x) {
  x[is.nan(x)] <- NA                  # convert NaN to NA if any
  x[is.na(x)] <- median(x, na.rm = TRUE)  # replace NAs with median
  return(x)
})
############################# Loading data + data processing #############################



############################# Removing 0 variance + converting binary to factor #############################
# Removing constant column (treatment)
constant_cols <- sapply(Dash_Data_filtered, function(x) length(unique(x)) == 1)
constant_cols
Dash_Data_filtered_no_const <- Dash_Data_filtered[, !constant_cols]


# Identify binary numeric columns (only 2 unique values, numeric)
binary_numeric_cols <- sapply(Dash_Data_filtered_no_const, function(x) {
  is.numeric(x) && length(unique(na.omit(x))) == 2
})


# Convert binary numeric columns to factor
Dash_Data_filtered_no_const[binary_numeric_cols] <- lapply(Dash_Data_filtered_no_const[binary_numeric_cols], as.factor)
############################# Removing 0 variance + converting binary to factor #############################




############################# Gower Distance + Silhouette Width #############################
gower_dist <- daisy(Dash_Data_filtered_no_const, metric = "gower")

# Finding optimal k clusters
avg_sil <- sapply(2:10, function(k) {
  pam_fit <- pam(gower_dist, k = k)
  mean(pam_fit$silinfo$avg.width)
})


# Plotting average silhouette width
k <- 2:10
silhouette_gg <- ggplot(data = data.frame(k, avg_sil), aes(x = k, y = avg_sil)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = which.max(avg_sil) + 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.10876340, lty = "dashed", color = "red") +
  annotate("text", x = which.max(avg_sil) + 1 + 0.3, 
           y = max(avg_sil) + 0.02, 
           label = paste("Optimal k =", which.max(avg_sil) + 1), hjust = 0) +
  labs(
    x = "Number of Clusters (k)",
    y = "Average Silhouette Width",
    title = "Silhouette Method for Optimal k (Mixed Data)"
  ) +
  theme_minimal() +
  theme(plot.margin = unit(c(1,1,1,1), "inches"),
        plot.title = element_text(face = "bold"))

silhouette_gg
############################# Gower Distance + Silhouette Width #############################



############################# PAM Clustering #############################
# Apply PAM clustering
pam_fit <- pam(gower_dist, k = 2)

# View cluster assignment
Dash_Data_filtered$pam_cluster <- pam_fit$clustering
#View(Dash_Data_filtered)


# Visualizing Gower clustering
# Perform MDS (cmdscale) on the Gower distance
mds_fit <- cmdscale(gower_dist, k = 2)  # 2D coordinates

# Create a data frame for plotting
mds_df <- data.frame(
  Dim1 = mds_fit[,1],
  Dim2 = mds_fit[,2],
  Cluster = factor(pam_fit$clustering)  # PAM cluster assignment
)

# Plot with ggplot2
G_clustering <- ggplot(mds_df, aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  stat_ellipse(aes(fill = Cluster), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  labs(title = "PAM Clusters visualized via MDS", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() +
  theme(plot.margin = unit(c(1,1,1,1), "inches"),
        plot.title = element_text(face = "bold"))

G_clustering
############################# PAM Clustering #############################



############################# k-prototypes clustering #############################
#Dash_Data_filtered_no_const <- Dash_Data_filtered_no_const %>% dplyr::select(-HC2CONSX)


# Run k-prototypes clustering (k = 2)
set.seed(123)
kpres <- kproto(Dash_Data_filtered_no_const, k = 2)
Dash_Data_filtered$k_proto_cluster <- factor(kpres$cluster)


mds_fit <- cmdscale(gower_dist, k = 2)

mds_df <- data.frame(
  Dim1 = mds_fit[,1],
  Dim2 = mds_fit[,2],
  cluster = Dash_Data_filtered$k_proto_cluster
)

# Plot clusters
k_prototype_clust <- ggplot(mds_df, aes(x = Dim1, y = Dim2, color = cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  stat_ellipse(aes(fill = cluster), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  labs(title = "K-Prototypes Clustering Visualization via MDS",
       x = "Dimension 1",
       y = "Dimension 2",
       color = "Cluster") +
  theme_classic() +
  theme(plot.margin = unit(c(1,1,1,1), "inches"),
        plot.title = element_text(face = "bold"))

k_prototype_clust
############################# k-prototypes clustering #############################




############################# Hierarchical Clustering #############################
# dist_matrix <- dist(numeric_data_imputed, method = "euclidean")
# hc <- hclust(dist_matrix, method = "ward.D2")
# plot(hc, labels = FALSE, hang = -1, main = "Hierarchical Clustering Dendrogram")
# rect.hclust(hc, k = 4, border = "red")
# 
# clusters_hc_5 <- cutree(hc, k = 4)
# pca <- prcomp(numeric_data_imputed[,-53], scale. = FALSE)
# pca_df <- as.data.frame(pca$x[, 1:2])
# pca_df$Cluster <- factor(clusters_hc_5)
# 
# hieararchical_clust <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
#   geom_point(size = 2) +
#   labs(title = "Hierarchical Clustering") +
#   theme_classic() +
#   theme(plot.margin = unit(c(1,1,1,1), "inches"),
#         plot.title = element_text(face = "bold"))
############################# Hierarchical Clustering #############################






############################# Merging all plots together #############################
combined_plot <- silhouette_gg / G_clustering / k_prototype_clust  # vertically stack the 3 plots
combined_plot_clusters_3 <- G_clustering + k_prototype_clust


ggsave("C:/Users/loisr/OneDrive/Documents/Pediatrics Neonatology - Dr. M/Projects/Clinical_hydrocortisone/combined_plot.pdf",
       plot = combined_plot, width = 8, height = 12)
ggsave("C:/Users/loisr/OneDrive/Documents/Pediatrics Neonatology - Dr. M/Projects/Clinical_hydrocortisone/combined_plot_clusters_3.pdf",
       plot = combined_plot_clusters_3, width = 12, height = 6)
############################# Merging all plots together #############################




#################################################################################################
#################################################################################################
############################# CANNOT HANDLE MISSING DATA INTERNALLY #############################
# Model-Based Clustering with mclust
numeric_data <- Dash_Data_filtered %>%
  dplyr::select(where(is.numeric))

# Fitting model-based clustering
mclust_fit <- Mclust(numeric_data)
summary(mclust_fit)
############################# CANNOT HANDLE MISSING DATA INTERNALLY #############################



############################# Latent Class Analysis (LCA) with poLCA #############################
# ONLY works with categorical data
cat_data <- Dash_Data_filtered %>%
  dplyr::select(where(is.factor))


f <- as.formula(paste("cbind(", paste(names(cat_data), collapse = ", "), ") ~ 1"))
set.seed(123)
lca_fit <- poLCA(f, data = cat_data, nclass = 3, na.rm = FALSE, maxiter = 1000)
Dash_Data_filtered$Cluster_LCA <- lca_fit$predclass


mca_res <- MCA(cat_data, graph = FALSE)
mca_coords <- as.data.frame(mca_res$ind$coord[, 1:2])
colnames(mca_coords) <- c("Dim1", "Dim2")
mca_coords$Cluster_LCA <- factor(Dash_Data_filtered$Cluster_LCA)


LCA_clustering <- ggplot(mca_coords, aes(x = Dim1, y = Dim2, color = Cluster_LCA)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "MCA Scatterplot colored by LCA clusters", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")
############################# Latent Class Analysis (LCA) with poLCA #############################




############################# Soft Imputation + Clustering #############################
# ONLY works for numeric data
#Convert to matrix and apply softImpute
numeric_matrix <- as.matrix(numeric_data)

# Running softImpute
fit_softimpute <- softImpute(numeric_matrix, rank.max = 50, lambda = 10)

#Impute missing values
completed_matrix <- complete(numeric_matrix, fit_softimpute)

# Run clustering on completed data, e.g. K-means
set.seed(123)
kmeans_fit <- kmeans(completed_matrix, centers = 3)

# Add cluster assignment to original data
Dash_Data_filtered$Cluster_SoftImpute <- kmeans_fit$cluster

# Perform PCA on the completed matrix
pca_res <- prcomp(completed_matrix, scale. = FALSE)

plot_data <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Cluster = factor(kmeans_fit$cluster)
)

softImpute_clustering <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.7, size = 3) +
  labs(title = "SoftImpute + K-means Clustering Visualization (PCA)",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")
############################# Soft Imputation + Clustering #############################
#################################################################################################
#################################################################################################

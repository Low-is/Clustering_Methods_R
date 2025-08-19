# Clustering Methods in R
A collection of clustering techniques in R for numerical, categorical, and mixed data. Includes methods for data preprocessing, handling missing values through imputation, and applying various clustering algorithms with practical examples.


Clustering is an unsupervised learning technique where you do not need labels. Instead, it finds natural groupings in the data, that is useful because:
- Clustering reveals hidden structure(s) such as subgroups.
    - Example: In patient data, clustering may reveal subtypes of disease progression or response to treatment.
- Before building predictive models, clustering helps you to explore patterns and generate hypotheses.
    - Example: Cluster patients by lab results and then investigate why certain cluster have better survival outcomes.
- Dimensionality reduction, reduces complex datasets into simpler groups and makes it easier to interpret and visualize large datasets.


## Types of data:
- Numeric data: Captures quantitative variation and helps to identify gradient, trends, and magnitudes.
- Categorical data: Captures qualitative distinctions and reveals group membership or structural differences.
- Binary data: cateogrical data with two outcomes. 

  
Working with mixed data can offer lots of advantages vs just only with numeric data or only working with categorical data, limiting analyses to only numeric or only categorical variables throws away valuable info. Mixed data makes clustering more powerful because it respects the complexity of real populations, captures complementary types of information, and reveals patterns that would remain hidden if you analyzed only one data type. Fortunately, there are clustering methods that can handle heterogenuous data. 

### The following clustering methods will be used in this repo:
- PAM (Partitioning Around Medoids) + Gower Distance: A clustering algorithm similar to K-means but more robust to noise and outliers. K-means is a clustering algorithm that computes the mean of a centroid (the mean of all points in the cluster), but instead of using the centroid, the medoid is used which represents the most centrally located actual data point. Additionally, this algortihm can handle mixed data types by apply different similarity rules:
      - Numeric -> scaled absolute difference
      - Categorical -> 0 if same, 1 if different
      - Binary -> similar logic to categorical variables

- K-prototypes Clustering: is a variant of k-means that can handle mixed data. It extends k-means (for numeric data) and k-modes (for categorical data). This clustering algortihm uses Euclidean distance for numeric features, and simple matching (count of mismatches) for categorical features.

- Hieararchical Clustering: clustering method that builds a hiearachy (tree structure) of clusters. Works only with numeric/continuous data. 
      - Agglomerative (bottom-up): start with each data point as its own cluster, then it iteratively merge the closest clusters until one big cluster remains.
      - Divisive (top-bottom): start with all points in one cluster, then iteratively split clusters.

# Load necessary libraries
library(readr)
library(ggplot2)
library(factoextra)
library(dplyr)
library(ggalt)  # For encircling clusters

# Create a function to save plots automatically
save_plot <- function(plot, filename) {
  ggsave(filename, plot = plot, width = 8, height = 6, dpi = 300)
}

# Define output folder
output_folder <- "enhanced_results/"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Load dataset
brightest_stars <- read_csv("C:/Users/User/Desktop/Debashis 2024/Astronomy Dataset/MAIN/brightest_stars.csv")

# Rename columns for consistency
colnames(brightest_stars) <- c("Traditional", "Bayer", "Luminosity", "Distance", "SpectralType", 
                               "ProperMotionRA", "ProperMotionDec", "RA_h", "RA_m", "Dec_deg", "Dec_min")

# Convert Spectral Type to numerical format for analysis
brightest_stars$SpectralTypeNum <- as.numeric(factor(brightest_stars$SpectralType))
brightest_stars$Distance_pc <- brightest_stars$Distance * 0.3066  # Convert light-years to parsecs

# Select numeric variables for clustering
numeric_data <- brightest_stars %>%
  select(Luminosity, Distance_pc, ProperMotionRA, ProperMotionDec)

# Set star names as row names for numeric_data so they appear as labels in dendrogram
rownames(numeric_data) <- brightest_stars$Traditional

#######################
# Hierarchical Clustering using Ward's method
dist_matrix <- dist(numeric_data)  # Compute distance matrix
hclust_result <- hclust(dist_matrix, method = "ward.D2")  # Apply Ward's method for clustering

# Plot dendrogram with star labels
dendrogram_plot <- fviz_dend(hclust_result, k = 3, rect = TRUE, rect_fill = TRUE, 
                             show_labels = TRUE, cex = 0.5, labels_track_height = 0.8) +
  ggtitle("Hierarchical Clustering Dendrogram with Star Labels") +
  xlab("Stars") +
  ylab("Height") +
  theme_minimal()

# Save the dendrogram plot
save_plot(dendrogram_plot, paste0(output_folder, "hierarchical_clustering_with_labels.png"))

#######################
# K-Means Clustering with Star Names
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(numeric_data, centers = 3)

# Add cluster assignment to the dataset
brightest_stars$Cluster <- factor(kmeans_result$cluster)

# Plot K-Means clustering results with ellipses and star names beside dots
kmeans_plot <- ggplot(brightest_stars, aes(x = Luminosity, y = SpectralTypeNum, color = Cluster)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_text(aes(label = Traditional), hjust = -0.2, vjust = 0.5, size = 3) +  # Add star names beside points
  geom_encircle(aes(group = Cluster), alpha = 0.2, color = "black", size = 1.5) +  # Ellipses around clusters
  ggtitle("K-Means Clustering: Luminosity vs Spectral Type with Star Names") +
  xlab("Luminosity (Sun = 1)") +
  ylab("Spectral Type (Numerical)") +
  theme_minimal()

# Save the K-Means plot
save_plot(kmeans_plot, paste0(output_folder, "kmeans_clusters_with_names.png"))

#######################
# Summary statistics and cluster table
# Create cluster summary table
cluster_summary <- brightest_stars %>%
  group_by(Cluster) %>%
  summarise(
    Mean_Luminosity = mean(Luminosity),
    SD_Luminosity = sd(Luminosity),
    Mean_Distance = mean(Distance_pc),
    SD_Distance = sd(Distance_pc),
    Mean_ProperMotionRA = mean(ProperMotionRA),
    SD_ProperMotionRA = sd(ProperMotionRA),
    Mean_ProperMotionDec = mean(ProperMotionDec),
    SD_ProperMotionDec = sd(ProperMotionDec)
  )

# Save the cluster summary table
write.csv(cluster_summary, file = paste0(output_folder, "cluster_summary.csv"), row.names = FALSE)

#######################
# Plotting Summary Statistics
# Summary plot of mean Luminosity and Distance for clusters
summary_plot <- ggplot(cluster_summary, aes(x = Mean_Luminosity, y = Mean_Distance, label = Cluster)) +
  geom_point(color = "blue", size = 5) +
  geom_text(vjust = -1, size = 4) +
  ggtitle("Cluster Summary: Mean Luminosity vs. Mean Distance") +
  xlab("Mean Luminosity") +
  ylab("Mean Distance (pc)") +
  theme_minimal()

# Save the summary plot
save_plot(summary_plot, paste0(output_folder, "cluster_summary_plot.png"))

#######################
# Display the plots
print(dendrogram_plot)
print(kmeans_plot)
print(summary_plot)


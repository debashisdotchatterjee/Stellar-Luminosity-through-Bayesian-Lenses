# Load necessary libraries
library(readr)
library(ggplot2)
library(factoextra)
library(dplyr)
library(ggalt)  # For encircling clusters
library(broom)  # For summarizing regression models

# Create a function to save plots automatically
save_plot <- function(plot, filename) {
  ggsave(filename, plot = plot, width = 10, height = 8, dpi = 300)
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

# Set star names as row names for numeric_data
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
# K-Means Clustering Plot
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(numeric_data, centers = 3)

# Add cluster assignment to the dataset
brightest_stars$Cluster <- factor(kmeans_result$cluster)

# Plot K-Means clustering results with ellipses and star names beside dots
kmeans_plot <- ggplot(brightest_stars, aes(x = Luminosity, y = SpectralTypeNum, color = Cluster)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_text(aes(label = Traditional), hjust = -0.1, vjust = 0.5, size = 3, check_overlap = TRUE) +  # Add star names
  geom_encircle(aes(group = Cluster), alpha = 0.2, color = "black", size = 1.5) +  # Ellipses around clusters
  ggtitle("K-Means Clustering: Luminosity vs Spectral Type with Star Names") +
  xlab("Luminosity (Sun = 1)") +
  ylab("Spectral Type (Numerical)") +
  theme_minimal()

# Save the K-Means plot
save_plot(kmeans_plot, paste0(output_folder, "kmeans_clusters_with_names.png"))

#######################
# Perform Nonlinear Regression Cluster-wise and Add Regression Line
# Create a list to store model summaries for each cluster
model_summaries <- list()

# Loop through each cluster and fit a nonlinear regression model
for (cluster in unique(brightest_stars$Cluster)) {
  
  # Subset data for the current cluster
  cluster_data <- brightest_stars %>% filter(Cluster == cluster)
  
  # Estimate reasonable starting values based on data
  a_start <- mean(cluster_data$Luminosity)
  b_start <- 0.5
  c_start <- 0.5
  
  # Ensure SpectralTypeNum has a meaningful value, filter out invalid data if necessary
  if(length(unique(cluster_data$SpectralTypeNum)) > 1) {
    # Nonlinear regression model: Luminosity as a function of Distance and SpectralTypeNum
    model <- tryCatch({
      nls(Luminosity ~ a * Distance_pc^b * SpectralTypeNum^c, 
          data = cluster_data, 
          start = list(a = a_start, b = b_start, c = c_start),
          control = nls.control(maxiter = 500))  # Increase max iterations
    }, error = function(e) {
      message(paste("Nonlinear model did not converge for cluster:", cluster))
      NULL
    })
  } else {
    message(paste("Insufficient variation in SpectralTypeNum for nonlinear model in cluster:", cluster))
    model <- NULL
  }
  
  # Check if the model converged
  if (!is.null(model)) {
    
    # Store the model summary
    model_summaries[[cluster]] <- glance(model)
    
    # Plot nonlinear regression fit for this cluster
    cluster_plot <- ggplot(cluster_data, aes(x = Distance_pc, y = Luminosity, color = factor(SpectralTypeNum))) +
      geom_point() +
      geom_smooth(method = "nls", formula = y ~ a * x^b * SpectralTypeNum^c, 
                  method.args = list(start = list(a = a_start, b = b_start, c = c_start)), 
                  se = FALSE) +
      ggtitle(paste("Nonlinear Regression Fit for Cluster", cluster)) +
      xlab("Distance (pc)") +
      ylab("Luminosity") +
      theme_minimal()
    
    # Save each cluster's regression plot
    save_plot(cluster_plot, paste0(output_folder, "cluster_", cluster, "_regression_plot.png"))
    
    # Display the plot
    print(cluster_plot)
    
  } else {
    # Fallback to linear regression if nonlinear fails
    message(paste("Falling back to linear regression for cluster:", cluster))
    
    # Fit linear regression: Luminosity ~ Distance + SpectralTypeNum
    model <- lm(Luminosity ~ Distance_pc + SpectralTypeNum, data = cluster_data)
    
    # Store the model summary
    model_summaries[[cluster]] <- glance(model)
    
    # Plot linear regression fit
    cluster_plot <- ggplot(cluster_data, aes(x = Distance_pc, y = Luminosity, color = factor(SpectralTypeNum))) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      ggtitle(paste("Linear Regression Fit for Cluster", cluster)) +
      xlab("Distance (pc)") +
      ylab("Luminosity") +
      theme_minimal()
    
    # Save the linear regression plot
    save_plot(cluster_plot, paste0(output_folder, "cluster_", cluster, "_linear_regression_plot.png"))
    
    # Display the plot
    print(cluster_plot)
  }
}

#######################
# Summarize regression goodness-of-fit for each cluster
regression_summary_table <- bind_rows(model_summaries, .id = "Cluster")

# Save the regression summary table as CSV
write.csv(regression_summary_table, file = paste0(output_folder, "regression_goodness_of_fit.csv"), row.names = FALSE)

#######################
# Display Final K-Means Plot and Dendrogram
print(kmeans_plot)
print(dendrogram_plot)

################################

# Loop through each cluster and fit a nonlinear regression model
for (cluster in unique(brightest_stars$Cluster)) {
  
  # Subset data for the current cluster
  cluster_data <- brightest_stars %>% filter(Cluster == cluster)
  
  # Estimate reasonable starting values based on data
  a_start <- mean(cluster_data$Luminosity, na.rm = TRUE)
  b_start <- 0.5
  c_start <- 0.5
  
  # Nonlinear regression model: Luminosity as a function of Distance and SpectralTypeNum
  model <- tryCatch({
    nls(Luminosity ~ a * Distance_pc^b * SpectralTypeNum^c, 
        data = cluster_data, 
        start = list(a = a_start, b = b_start, c = c_start),
        control = nls.control(maxiter = 500))  # Increase max iterations
  }, error = function(e) {
    message(paste("Nonlinear model did not converge for cluster:", cluster))
    NULL
  })
  
  # Check if the model converged
  if (!is.null(model)) {
    
    # Store the model summary
    model_summaries[[cluster]] <- glance(model)
    
    # Get model-predicted values for plotting the curve
    cluster_data$predicted <- predict(model)
    
    # Plot nonlinear regression fit for this cluster
    cluster_plot <- ggplot(cluster_data, aes(x = Distance_pc, y = Luminosity, color = SpectralTypeNum)) +
      geom_point(size = 4) +  # Use larger points
      geom_line(aes(y = predicted), color = "red", size = 1, linetype = "dotted") +  # Plot dotted regression curve
      ggtitle(paste("Nonlinear Regression Fit for Cluster", cluster)) +
      xlab("Distance (pc)") +
      ylab("Luminosity") +
      theme_minimal()
    
    # Save each cluster's regression plot
    save_plot(cluster_plot, paste0(output_folder, "cluster_", cluster, "_regression_plot.png"))
    
    # Display the plot
    print(cluster_plot)
    
  } else {
    # Fallback to linear regression if nonlinear fails
    message(paste("Falling back to linear regression for cluster:", cluster))
    
    # Fit linear regression: Luminosity ~ Distance + SpectralTypeNum
    model <- lm(Luminosity ~ Distance_pc + SpectralTypeNum, data = cluster_data)
    
    # Store the model summary
    model_summaries[[cluster]] <- glance(model)
    
    # Get model-predicted values for plotting the line
    cluster_data$predicted <- predict(model)
    
    # Plot linear regression fit
    cluster_plot <- ggplot(cluster_data, aes(x = Distance_pc, y = Luminosity, color = SpectralTypeNum)) +
      geom_point(size = 4) +  # Use larger points
      geom_line(aes(y = predicted), color = "blue", size = 1, linetype = "dotted") +  # Plot dotted linear regression line
      ggtitle(paste("Linear Regression Fit for Cluster", cluster)) +
      xlab("Distance (pc)") +
      ylab("Luminosity") +
      theme_minimal()
    
    # Save the linear regression plot
    save_plot(cluster_plot, paste0(output_folder, "cluster_", cluster, "_linear_regression_plot.png"))
    
    # Display the plot
    print(cluster_plot)
  }
}

# Summarize regression goodness-of-fit for each cluster
regression_summary_table <- bind_rows(model_summaries, .id = "Cluster")

# Save the regression summary table as CSV
write.csv(regression_summary_table, file = paste0(output_folder, "regression_goodness_of_fit.csv"), row.names = FALSE)


##########################
# Define consistent color palette for K-means clusters
kmeans_colors <- c("red", "blue", "green")  # Make sure to match the K-means plot color scheme

# 1. Hierarchical Clustering using Ward's method
dist_matrix <- dist(numeric_data)  # Compute distance matrix
hclust_result <- hclust(dist_matrix, method = "ward.D2")  # Apply Ward's method for clustering

# Cut the hierarchical dendrogram into 3 clusters
hclust_clusters <- cutree(hclust_result, k = 3)  # Create 3 clusters based on dendrogram cut

# Match hierarchical clustering colors with new color (different from k-means)
hclust_line_colors <- c("purple", "orange", "cyan")  # Different colors for the dendrogram lines

# Plot dendrogram with K-means consistent colors for points and separate color for dendrogram lines
dendrogram_plot <- fviz_dend(hclust_result, k = 3, rect = TRUE, rect_fill = TRUE, 
                             show_labels = TRUE, cex = 0.5,
                             k_colors = hclust_line_colors) +  # Different colors for lines
  ggtitle("Hierarchical Clustering Dendrogram with Separate Line Colors") +
  xlab("Stars") +
  ylab("Height") +
  theme_minimal()

# Save the hierarchical dendrogram plot
save_plot(dendrogram_plot, paste0(output_folder, "hierarchical_clustering_separate_line_colors.png"))

# 2. K-Means Clustering plot (points with consistent K-means colors)
kmeans_plot <- ggplot(brightest_stars, aes(x = Luminosity, y = SpectralTypeNum, color = Cluster)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_text(aes(label = Traditional), hjust = -0.1, vjust = 0.5, size = 3, check_overlap = TRUE) +  # Add star names
  geom_encircle(aes(group = Cluster), alpha = 0.2, color = "black", size = 1.5) +  # Ellipses around clusters
  ggtitle("K-Means Clustering: Luminosity vs Spectral Type with Star Names") +
  xlab("Luminosity (Sun = 1)") +
  ylab("Spectral Type (Numerical)") +
  scale_color_manual(values = kmeans_colors) +  # Apply consistent K-means colors
  theme_minimal()

# Save the K-Means plot
save_plot(kmeans_plot, paste0(output_folder, "kmeans_clusters_consistent_colors.png"))

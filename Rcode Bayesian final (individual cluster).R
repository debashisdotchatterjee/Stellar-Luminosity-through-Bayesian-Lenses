# Load necessary libraries
library(tidyverse)
library(brms)       # For Bayesian modeling
library(ggplot2)
library(bayesplot) # For posterior predictive checks

library(brms)   
library(dplyr)

# Define output folder to save plots and tables
output_folder <- "cluster_bayesian_results/"
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

# Clustering (assuming K-Means clusters for separation)
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(brightest_stars[,c("Luminosity", "Distance_pc", "ProperMotionRA", "ProperMotionDec")], centers = 3)

# Add cluster information to the dataset
brightest_stars$Cluster <- as.factor(kmeans_result$cluster)

# Function to save plots
save_plot <- function(plot, filename) {
  ggsave(filename, plot = plot, width = 10, height = 8, dpi = 300)
}

##################################
# Cluster 1: Bayesian Model, Plots and Summaries
##################################
cluster_1_data <- brightest_stars %>% filter(Cluster == 1)

# Define the Bayesian hierarchical model for Cluster 1
bayesian_model_1 <- brm(
  formula = Luminosity ~ (1 | SpectralTypeNum) + Distance_pc + ProperMotionRA + ProperMotionDec,
  data = cluster_1_data,
  family = gaussian(),
  prior = c(
    prior(normal(0, 10), class = "b"),
    prior(cauchy(0, 2), class = "sd")  # Cauchy prior for the group-level standard deviation
  ),
  iter = 10000,  # Increase iterations for better sampling
  chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  cores = 4
)

# Posterior predictive check for Cluster 1
ppc_plot_1 <- pp_check(bayesian_model_1, ndraws = 100) + 
  ggtitle("Posterior Predictive Check for Cluster 1") +
  theme_minimal()

# Save and print the posterior predictive plot for Cluster 1
save_plot(ppc_plot_1, paste0(output_folder, "cluster_1_ppc.png"))

# Fitted vs Actual for Cluster 1
fitted_luminosity_1 <- apply(posterior_epred(bayesian_model_1), 2, mean)

fitted_actual_plot_1 <- ggplot(cluster_1_data, aes(x = Luminosity, y = fitted_luminosity_1)) +
  geom_point(color = 'blue', alpha = 0.7, size = 3) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dotted") +
  ggtitle("Fitted vs Actual Luminosity for Cluster 1") +
  xlab("Actual Luminosity") +
  ylab("Fitted Luminosity") +
  theme_minimal()

# Save and print fitted vs actual plot for Cluster 1
save_plot(fitted_actual_plot_1, paste0(output_folder, "cluster_1_fitted_actual.png"))

# Print model summary for Cluster 1
print(summary(bayesian_model_1))

##################################
# Cluster 2: Bayesian Model, Plots and Summaries
##################################
cluster_2_data <- brightest_stars %>% filter(Cluster == 2)

# Define the Bayesian hierarchical model for Cluster 2
bayesian_model_2 <- brm(
  formula = Luminosity ~ (1 | SpectralTypeNum) + Distance_pc + ProperMotionRA + ProperMotionDec,
  data = cluster_2_data,
  family = gaussian(),
  prior = c(
    prior(normal(0, 10), class = "b"),
    prior(cauchy(0, 2), class = "sd")
  ),
  iter = 10000,  # Increase iterations for better sampling
  chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  cores = 4
)

# Posterior predictive check for Cluster 2
ppc_plot_2 <- pp_check(bayesian_model_2, ndraws = 100) + 
  ggtitle("Posterior Predictive Check for Cluster 2") +
  theme_minimal()

# Save and print the posterior predictive plot for Cluster 2
save_plot(ppc_plot_2, paste0(output_folder, "cluster_2_ppc.png"))

# Fitted vs Actual for Cluster 2
fitted_luminosity_2 <- apply(posterior_epred(bayesian_model_2), 2, mean)

fitted_actual_plot_2 <- ggplot(cluster_2_data, aes(x = Luminosity, y = fitted_luminosity_2)) +
  geom_point(color = 'blue', alpha = 0.7, size = 3) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dotted") +
  ggtitle("Fitted vs Actual Luminosity for Cluster 2") +
  xlab("Actual Luminosity") +
  ylab("Fitted Luminosity") +
  theme_minimal()

# Save and print fitted vs actual plot for Cluster 2
save_plot(fitted_actual_plot_2, paste0(output_folder, "cluster_2_fitted_actual.png"))

# Print model summary for Cluster 2
print(summary(bayesian_model_2))

##################################
# Cluster 3: Bayesian Model, Plots and Summaries
##################################
cluster_3_data <- brightest_stars %>% filter(Cluster == 3)

# Define the Bayesian hierarchical model for Cluster 3
bayesian_model_3 <- brm(
  formula = Luminosity ~ (1 | SpectralTypeNum) + Distance_pc + ProperMotionRA + ProperMotionDec,
  data = cluster_3_data,
  family = gaussian(),
  prior = c(
    prior(normal(0, 10), class = "b"),
    prior(cauchy(0, 2), class = "sd")
  ),
  iter = 10000,  # Increase iterations for better sampling
  chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  cores = 4
)

# Posterior predictive check for Cluster 3
ppc_plot_3 <- pp_check(bayesian_model_3, ndraws = 100) + 
  ggtitle("Posterior Predictive Check for Cluster 3") +
  theme_minimal()

# Save and print the posterior predictive plot for Cluster 3
save_plot(ppc_plot_3, paste0(output_folder, "cluster_3_ppc.png"))

# Fitted vs Actual for Cluster 3
fitted_luminosity_3 <- apply(posterior_epred(bayesian_model_3), 2, mean)

fitted_actual_plot_3 <- ggplot(cluster_3_data, aes(x = Luminosity, y = fitted_luminosity_3)) +
  geom_point(color = 'blue', alpha = 0.7, size = 3) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dotted") +
  ggtitle("Fitted vs Actual Luminosity for Cluster 3") +
  xlab("Actual Luminosity") +
  ylab("Fitted Luminosity") +
  theme_minimal()

# Save and print fitted vs actual plot for Cluster 3
save_plot(fitted_actual_plot_3, paste0(output_folder, "cluster_3_fitted_actual.png"))

# Print model summary for Cluster 3
print(summary(bayesian_model_3))

###################################

# Additional diagnostics: Trace plots and ESS measures for each cluster
diagnostics_plots <- list()

for (cluster in unique(brightest_stars$Cluster)) {
  
  # Filter the cluster-specific Bayesian model summaries
  bayesian_model <- get(paste0("bayesian_model_", cluster))
  
  # Trace plot for each parameter in the model
  trace_plot <- mcmc_trace(as.array(bayesian_model), regex_pars = "b_") +
    ggtitle(paste("Trace Plot for Cluster", cluster)) +
    theme_minimal()
  
  # Save trace plot
  save_plot(trace_plot, paste0(output_folder, "cluster_", cluster, "_trace_plot.png"))
  diagnostics_plots[[paste0("cluster_", cluster, "_trace_plot")]] <- trace_plot
  
  # Print trace plot
  print(trace_plot)
  
  # Effective Sample Size (ESS) diagnostic
  ess_plot <- mcmc_areas(as.array(bayesian_model), regex_pars = "b_", prob_outer = 0.95) +
    ggtitle(paste("Effective Sample Size (ESS) for Cluster", cluster)) +
    theme_minimal()
  
  # Save ESS plot
  save_plot(ess_plot, paste0(output_folder, "cluster_", cluster, "_ess_plot.png"))
  diagnostics_plots[[paste0("cluster_", cluster, "_ess_plot")]] <- ess_plot
  
  # Print ESS plot
  print(ess_plot)
  
  # Summary ESS diagnostic printed in console
  ess_summary <- summary(bayesian_model)$summary
  cat("\n\nEffective Sample Size (ESS) Summary for Cluster", cluster, ":\n")
  print(ess_summary)
}

# Save all the diagnostics summary and plots
for (plot_name in names(diagnostics_plots)) {
  print(diagnostics_plots[[plot_name]])
}

######################
# Load necessary libraries
library(ggplot2)
library(brms)
library(dplyr)

# Create a function to calculate and display the correlation in the plot
add_correlation <- function(data, fitted, actual) {
  cor_value <- cor(data[[actual]], data[[fitted]], use = "complete.obs")
  return(paste("Correlation:", round(cor_value, 2)))
}

# Loop through clusters and generate fitted vs actual plots
for (cluster in unique(brightest_stars$Cluster)) {
  
  # Filter data for the current cluster
  cluster_data <- brightest_stars %>% filter(Cluster == cluster)
  
  # Define the Bayesian model for this cluster
  bayesian_model <- brm(
    formula = Luminosity ~ (1 | SpectralTypeNum) + Distance_pc + ProperMotionRA + ProperMotionDec,
    data = cluster_data,
    family = gaussian(),
    prior = c(
      prior(normal(0, 10), class = "b"),
      prior(cauchy(0, 2), class = "sd")
    ),
    iter = 4000,
    chains = 4,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    cores = 4
  )
  
  # Extract fitted values from the Bayesian model
  fitted_luminosity <- posterior_epred(bayesian_model) %>% apply(2, mean)
  
  # Add the fitted values to the cluster_data
  cluster_data$fitted_luminosity <- fitted_luminosity
  
  # Calculate the correlation value
  cor_text <- add_correlation(cluster_data, "fitted_luminosity", "Luminosity")
  
  # Create the plot
  fitted_actual_plot <- ggplot(cluster_data, aes(x = Luminosity, y = fitted_luminosity)) +
    geom_point(color = 'blue', alpha = 0.7, size = 4) +  # Larger points
    geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dotted", size = 1.2) +  # Prominent dotted line
    ggtitle(paste("Fitted vs Actual Luminosity for Cluster", cluster)) +
    xlab("Actual Luminosity") +
    ylab("Fitted Luminosity") +
    annotate("text", x = Inf, y = -Inf, label = cor_text, vjust = -1.5, hjust = 1.1, size = 5, color = "black") +  # Append correlation value
    theme_minimal(base_size = 16)
  
  # Save the plot
  save_plot(fitted_actual_plot, paste0(output_folder, "cluster_", cluster, "_fitted_actual_with_correlation.png"))
  
  # Print the plot
  print(fitted_actual_plot)
}


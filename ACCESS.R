#### ACCESS ######
raster_data <- rast("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ACCESS-ESM1-5_1950-1990_regridded.nc")

spi_12_matrix <- matrix(NA, nrow = ncell(raster_data), ncol = 492)


for (i in 1:ncell(raster_data)) {
  cell_values <- extract(raster_data, i)  
  
  spi_12_result <- spi(as.numeric(cell_values), scale = 12, 
                       distribution = "Gamma", fit = "ub-pwm", 
                       na.rm = TRUE, verbose = FALSE)
  
  spi_12_matrix[i, ] <- spi_12_result$fitted
}

spi_raster <- raster_data  
values(spi_raster) <- spi_12_matrix  


writeCDF(spi_raster, "C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/SPI12_ACCESS-ESM1-5_1950-1990.nc", 
         overwrite = TRUE, varname = "spi12", 
         longname = "SPI-12", zname = "time")

################ Plot

spi_access <- rast("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/SPI12_ACCESS-ESM1-5_1950-1990.nc")
print(spi_access)

plot(spi_access[[30]])



##### PCA, SVD and Varimax Rotation #####

# Note: The spi_12_matrix currently has dimensions: (ncell, 492),
# meaning each row is a spatial cell and each column is a time lag.
# For PCA analogous to the Python script (observations = time steps),
# we need to transpose the matrix so that rows represent time and columns represent spatial cells.
spi_data <- t(spi_12_matrix)  # New dimensions: 492 (time steps) x ncell (spatial cells)

# Remove rows with any NA values
spi_data_clean <- spi_data[complete.cases(spi_data), ]

# Check dimensions after removal
dim(spi_data_clean)


# Now perform PCA on the cleaned data
pca_res <- prcomp(spi_data_clean, center = TRUE, scale. = FALSE)

# Set the maximum number of components to retain (we took 60 components)
max_comps <- min(100, ncol(pca_res$rotation))
pca_loadings <- pca_res$rotation[, 1:max_comps, drop = FALSE]  # Spatial loadings (ncell x max_comps)
pca_scores   <- pca_res$x[, 1:max_comps, drop = FALSE]           # Temporal scores (time steps x max_comps)

# Apply Varimax rotation on the PCA loadings.
varimax_res <- varimax(pca_loadings)
rotated_loadings <- varimax_res$loadings             # Rotated loadings (ncell x max_comps)
rotated_scores   <- pca_scores %*% varimax_res$rotmat   # Rotated time series (time steps x max_comps)

# Calculate the variance of each rotated component and the fraction of total variance explained.
rotated_variance <- apply(rotated_scores, 2, var)
total_variance <- sum(pca_res$sdev[1:max_comps]^2)
rotated_explained <- rotated_variance / total_variance

# For each rotated component, find the cell with the maximum absolute loading.
coords <- xyFromCell(raster_data, 1:ncell(raster_data))
max_idx <- apply(abs(rotated_loadings), 2, which.max)
max_coords <- coords[max_idx, ]

# Print summary results.
print("Explained variance of rotated components:")
print(rotated_explained)
print("Coordinates of maximum absolute loadings for each component:")
print(max_coords)

# Suppose you have already computed rotated_explained and cumulative variance:
cumulative_variance <- cumsum(rotated_explained)
print(cumulative_variance)

# Decide on a cutoff (e.g., 90% of the variance)
cutoff <- 0.9
n_components <- which(cumulative_variance >= cutoff)[1]
cat("Number of components selected:", n_components, "\n")

# Subset the rotated scores (rotated_scores has dimensions: time steps x components)
subset_rotated_scores <- rotated_scores[, 1:n_components, drop = FALSE]

# Assuming 'subset_rotated_scores' is your matrix (or data frame) of selected components.
# If you want to include a time index, ensure it's part of the data frame:
subset_data <- as.data.frame(subset_rotated_scores)
subset_data$time <- seq_len(nrow(subset_data))

# Write the data to a CSV file.
write.csv(subset_data, "ACCESS_subset_rotated_scores.csv", row.names = FALSE)

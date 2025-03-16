library(terra)
library(SPEI)
library(data.table)
library(arrow)
library(maps)
library(ggplot2)
library(ggrepel)
library(cowplot)

##### DATA ####
era5 <- rast("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5/1950_1990_monthly_mm.nc")


common_extent <- ext(
  max(c(ext(era5)[1], ext(access)[1])),   
  min(c(ext(era5)[2], ext(access)[2])),   
  max(c(ext(era5)[3], ext(access)[3])),   
  min(c(ext(era5)[4], ext(access)[4]))    
)

target_res <- rast(common_extent, resolution = c(2, 2))
print(target_res)

era5_regridded <- resample(era5, target_res, method = "bilinear")


writeCDF(era5_regridded, "C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5/ERA5_1950-1990_regridded.nc", overwrite = TRUE)

plot(era5_regridded[[1]])
print(era5_regridded)

##### SPI #####
raster_data <- rast("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5/ERA5_1950-1990_regridded.nc")
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

print(spi_12_matrix)
plot(spi_12_result)
writeCDF(spi_raster, "C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5/SPI12_ERA5_1950-1990.nc", 
         overwrite = TRUE, varname = "spi12", 
         longname = "SPI-12", zname = "time")
spi_era5 <- rast("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5/SPI12_ERA5_1950-1990.nc")
print(spi_era5)

plot(spi_era5[[20]])
# Save SPI raster plot
png(file.path("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5", "spi_era5_plot.png"), width = 800, height = 600)
plot(spi_era5[[20]], main = "SPI-12 ERA5 (Time Step 20)")
dev.off()
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
write.csv(subset_data, "C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5/ERA5_subset_rotated_scores.csv", row.names = FALSE)





####### BARPLOT ######
# Multiply by 100 for percentage representation
barplot(rotated_explained * 100, 
        names.arg = 1:length(rotated_explained),
        xlab = "Rotated Component", ylab = "Variance Explained (%)", 
        main = "Scree Plot of Rotated Components", col = "lightblue")
png(file.path("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5", "rotated_components_barplot.png"), width = 800, height = 600)
barplot(rotated_explained * 100, 
        names.arg = 1:length(rotated_explained),
        xlab = "Rotated Component", ylab = "Variance Explained (%)", 
        main = "Scree Plot of Rotated Components", col = "lightblue")
dev.off()



####### FIRST COMPONENTS MAP ######
# For each cell (row), find the index of the component with the largest absolute loading
component_assignments <- apply(rotated_loadings, 1, function(x) which.max(abs(x)))

# Create a new SpatRaster with these integer assignments
assign_raster <- raster_data
values(assign_raster) <- component_assignments
assign_raster <- rotate(assign_raster)



# Convert the SpatRaster to a data frame
df_assign <- as.data.frame(assign_raster, xy = TRUE)
# The default column name for the raster values might be something like "layer"
# Rename it to "component"
names(df_assign)[3] <- "component"


# 1) Convert 'component' to numeric instead of factor
df_assign$component <- as.numeric(df_assign$component)

# 2) Plot with a continuous viridis color scale
area_components_map <- ggplot() +
  geom_raster(data = df_assign, aes(x = x, y = y, fill = component)) +
  geom_path(data = world, aes(x = long, y = lat, group = group), 
            color = "black", size = 0.2) +
  scale_fill_viridis_c(
    name = "Component",
    # Optional: set limits & breaks if desired
    # limits = c(1, max(df_assign$component)), 
    # breaks = 1:max(df_assign$component),
    guide = guide_colorbar(barwidth = 1, barheight = 10)
  ) +
  coord_fixed(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  theme_minimal() +
  labs(
    title = "Dominant Rotated Component by Cell",
    x = "Longitude", 
    y = "Latitude"
  )
ggsave(file.path("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5", "area_dominant_components_map.png"), area_components_map, width = 10, height = 6)

####### SECOND COMPONENTS MAP #########

# max_idx[i] = the cell index with the largest absolute loading for component i
max_idx <- apply(abs(rotated_loadings), 2, which.max)
coords <- xyFromCell(raster_data, 1:ncell(raster_data))
max_coords <- coords[max_idx, ]  # lat/lon for each component

# for nicer text labels
# Shift negative longitudes from [-180, 0) to [180, 360)
df_nodes$x_360 <- ifelse(df_nodes$x < 0, df_nodes$x + 360, df_nodes$x)
components_map <-ggplot() +
  geom_polygon(data = world_360, 
               aes(x = long_360, y = lat, group = group),
               fill = "white", color = "black", size = 0.2) +
  geom_point(data = df_nodes, 
             aes(x = x_360, y = y),
             color = "red", size = 3) +
  geom_text_repel(data = df_nodes, 
                  aes(x = x_360, y = y, label = component),
                  color = "black", size = 4) +
  coord_fixed(xlim = c(0, 360), ylim = c(-90, 90)) +
  theme_minimal() +
  labs(title = "Primary Node for Each Rotated Component (0–360)",
       x = "Longitude (0–360)", y = "Latitude")
ggsave(file.path("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5", "dominant_components_map.png"), components_map, width = 10, height = 6)

######## Time series ########
time <- seq_len(nrow(rotated_scores))  # each row in rotated_scores is a time step
df_ts <- data.frame(
  time = time,
  comp1 = rotated_scores[, 1],
  comp2 = rotated_scores[, 2],
  comp3 = rotated_scores[, 3]
  # etc. if you want more
)

p_left <- ggplot(df_ts, aes(x = time)) +
  geom_line(aes(y = comp1, color = "Component 1")) +
  geom_line(aes(y = comp2, color = "Component 2")) +
  geom_line(aes(y = comp3, color = "Component 3")) +
  scale_color_manual(values = c("blue", "red", "darkgreen")) +
  theme_minimal() +
  labs(x = "Time (index)", y = "Rotated Score", color = "Series",
       title = "Time Series of Selected Components")

# 5B. Map showing maximum-loading points for these 3 components
df_nodes_subset <- data.frame(
  component = 1:3,
  x = max_coords[1:3, 1],
  y = max_coords[1:3, 2]
)


# Combine side by side
plot(p_left)
ggsave(file.path("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5", "time_series_plot.png"), p_left, width = 10, height = 6)




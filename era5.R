library(terra)
library(SPEI)
library(data.table)
library(arrow)

##### DATA ####
era5 <- rast("1950_1990_monthly_mm.nc")
access <- rast("pr_monthly_sums_ACCESS-ESM1-5_1950-1990.nc")

common_extent <- ext(
  max(c(ext(era5)[1], ext(access)[1])),   
  min(c(ext(era5)[2], ext(access)[2])),   
  max(c(ext(era5)[3], ext(access)[3])),   
  min(c(ext(era5)[4], ext(access)[4]))    
)

target_res <- rast(common_extent, resolution = c(2, 2))
print(target_res)

era5_regridded <- resample(era5, target_res, method = "bilinear")
access_regridded <- resample(access, target_res, method = "bilinear")

writeCDF(access_regridded, "C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ACCESS-ESM1-5_1950-1990_regridded.nc", overwrite = TRUE)
writeCDF(era5_regridded, "C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5_1950-1990_regridded.nc", overwrite = TRUE)

plot(access_regridded[[1]])
plot(era5_regridded[[1]])
print(access_regridded)
print(era5_regridded)

##### SPI #####
raster_data <- rast("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5_1950-1990_regridded.nc")
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
writeCDF(spi_raster, "C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/SPI12_ERA5_1950-1990.nc", 
         overwrite = TRUE, varname = "spi12", 
         longname = "SPI-12", zname = "time")
spi_era5 <- rast("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/SPI12_ERA5_1950-1990.nc")
print(spi_era5)

plot(spi_era5[[20]])
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
max_comps <- min(60, ncol(pca_res$rotation))
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
write.csv(subset_data, "ERA5_subset_rotated_scores.csv", row.names = FALSE)





########## visualisations ###########


####### BARPLOT ######
# Multiply by 100 for percentage representation
barplot(rotated_explained * 100, 
        names.arg = 1:length(rotated_explained),
        xlab = "Rotated Component", ylab = "Variance Explained (%)", 
        main = "Scree Plot of Rotated Components", col = "lightblue")


####### DOES NOT WORK YET ######
# Install required packages if not already installed:
#install.packages(c("leaflet", "terra", "raster", "viridisLite"))

library(leaflet)
library(terra)
library(raster)  # For conversion compatibility with leaflet
library(viridisLite)

# Select the component to visualize (e.g., first component)
component_index <- 1

# Create a raster for the selected component using your template raster
# Assuming 'raster_data' is a SpatRaster from terra and rotated_loadings is a matrix
component_raster <- raster_data  # Using your existing terra object
terra::values(component_raster) <- rotated_loadings[, component_index]

print(raster_data)

# (Optional) If your data are not in EPSG:4326 (WGS84), you should transform them.
# For example, if they are in another CRS:
# component_raster <- project(component_raster, "EPSG:4326")

# Convert terra raster to a raster object from the 'raster' package for compatibility with leaflet.
component_raster_sp <- raster(component_raster)

# Create a color palette based on the raster values.
pal <- colorNumeric(palette = viridis(100), domain = values(component_raster_sp), na.color = "transparent")

# Create the interactive map.
leaflet() %>%
  addTiles() %>%  # Base map layer
  addRasterImage(component_raster_sp, colors = pal, opacity = 0.7) %>%
  addLegend(pal = pal, values = values(component_raster_sp), title = paste("Component", component_index, "Loadings")) %>%
  # Overlay a marker at the maximum absolute loading point for the selected component.
  addMarkers(
    lng = max_coords[component_index, "x"],
    lat = max_coords[component_index, "y"],
    popup = paste("Max loading for Component", component_index)
  )


########## PLOTLY ######### (WORKS)
# Install plotly if needed:
#install.packages("plotly")

library(plotly)
library(ggplot2)
# Create a data frame containing time (e.g., time steps) and rotated scores for two components.
time <- 1:nrow(rotated_scores)  # Assuming rows correspond to time steps
data_ts <- data.frame(
  Time = time,
  Component1 = rotated_scores[, 1],
  Component2 = rotated_scores[, 2]
)

# Build the interactive time series plot.
p <- plot_ly(data_ts, x = ~Time) %>%
  add_lines(y = ~Component1, name = "Component 1", line = list(width = 2)) %>%
  add_lines(y = ~Component2, name = "Component 2", line = list(width = 2)) %>%
  layout(
    title = "Interactive Time Series of Rotated Scores",
    xaxis = list(title = "Time Step"),
    yaxis = list(title = "Rotated Score")
  )

# Display the plot.
p






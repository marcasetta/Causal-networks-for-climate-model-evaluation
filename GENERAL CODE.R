###### GENERAL FILE ######

# PACKAGES 
#install.packages("arrow")
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


writeCDF(access_regridded, "C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ACCESS-ESM1-5_1950-1990_regridded.nc", overwrite=TRUE)
writeCDF(era5_regridded, "C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5_1950-1990_regridded.nc", overwrite=TRUE)

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


writeCDF(spi_raster, "C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/SPI12_ERA5_1950-1990.nc", 
         overwrite = TRUE, varname = "spi12", 
         longname = "SPI-12", zname = "time")

################ Plot

spi_era5 <- rast("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/SPI12_ERA5_1950-1990.nc")
print(spi_era5)

plot(spi_era5[[20]])

########################### access

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

######## change data structure for fitting python ##########


#### ERA DATA ####

# Load the original NetCDF file
raster_data <- rast("C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/ERA5_1950-1990_regridded.nc")

# Get raster dimensions
n_cells <- ncell(raster_data)  # Number of spatial locations
n_time <- 492  # Assuming monthly data from 1950 to 1990

# Initialize an empty matrix (T, N) for PCMCI
spi_12_matrix <- matrix(NA, nrow = n_time, ncol = n_cells)

# Compute SPI-12 for each cell
for (i in 1:n_cells) {
  cell_values <- extract(raster_data, i)  # Extract time series for cell
  
  spi_12_result <- spi(as.numeric(cell_values), scale = 12, 
                       distribution = "Gamma", fit = "ub-pwm", 
                       na.rm = TRUE, verbose = FALSE)
  
  spi_12_matrix[, i] <- spi_12_result$fitted  # Store result in (T, N) format
}

plot(spi_12_result)
# Convert matrix to a DataTable (efficient format)
spi_12_dt <- as.data.table(spi_12_matrix)

# Save as CSV (easiest format for Python)
write.csv(spi_12_dt, "C:/Users/marta/OneDrive/Documenti/UNI/Climate change seminar/processed_spi12.csv", row.names = FALSE)


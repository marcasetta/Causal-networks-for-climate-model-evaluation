library(terra)
library(SPEI)


raster_data <- rast("C:/Users/Jonas/Desktop/StatistikSeminar/ERA5_1950-1990_regridded.nc")

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


writeCDF(spi_raster, "C:/Users/Jonas/Desktop/StatistikSeminar/SPI12_ERA5_1950-1990.nc", 
         overwrite = TRUE, varname = "spi12", 
         longname = "SPI-12", zname = "time")

################ Plot

spi_era5 <- rast("C:/Users/Jonas/Desktop/StatistikSeminar/SPI12_ERA5_1950-1990.nc")
print(spi_era5)

plot(spi_era5[[20]])

########################### access

raster_data <- rast("C:/Users/Jonas/Desktop/StatistikSeminar/ACCESS-ESM1-5_1950-1990_regridded.nc")

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


writeCDF(spi_raster, "C:/Users/Jonas/Desktop/StatistikSeminar/SPI12_ACCESS-ESM1-5_1950-1990.nc", 
         overwrite = TRUE, varname = "spi12", 
         longname = "SPI-12", zname = "time")

################ Plot

spi_access <- rast("C:/Users/Jonas/Desktop/StatistikSeminar/SPI12_ACCESS-ESM1-5_1950-1990.nc")
print(spi_access)

plot(spi_access[[30]])


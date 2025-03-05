library(terra)


era5 <- rast("C:/Users/Jonas/Desktop/StatistikSeminar/1950_1990_monthly_mm.nc")
access <- rast("C:/Users/Jonas/Desktop/StatistikSeminar/pr_monthly_sums_ACCESS-ESM1-5_1950-1990.nc")


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


writeCDF(access_regridded, "C:/Users/Jonas/Desktop/StatistikSeminar/ACCESS-ESM1-5_1950-1990_regridded.nc", overwrite=TRUE)
writeCDF(era5_regridded, "C:/Users/Jonas/Desktop/StatistikSeminar/ERA5_1950-1990_regridded.nc", overwrite=TRUE)

plot(access_regridded[[1]])
plot(era5_regridded[[1]])

print(access_regridded)
print(era5_regridded)





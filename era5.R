library(terra)

# Setze den Pfad zu deinen NetCDF-Dateien
data_path <- "C:/Users/Jonas/Desktop/StatistikSeminar/"

# Liste der Jahre
years <- 1950:1990

# Funktion zur Bestimmung der Tage pro Monat (Schaltjahr berücksichtigen)
days_in_month <- function(year) {
  if ((year %% 4 == 0 && year %% 100 != 0) || (year %% 400 == 0)) {
    return(c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))  # Schaltjahr
  } else {
    return(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))  # Normales Jahr
  }
}

# Schleife über alle Jahre
for (year in years) {
  file_path <- paste0(data_path, year, ".nc")
  
  if (!file.exists(file_path)) {
    warning(paste("Datei fehlt:", file_path))
    next
  }
  
  cat("Verarbeite Jahr:", year, "\n")  # Fortschritt ausgeben
  
  daily_raster <- rast(file_path)
  
  if (is.na(crs(daily_raster))) {
    crs(daily_raster) <- "EPSG:4326"
  }
  
  days_per_month <- days_in_month(year)
  
  monthly_sums <- list()
  start_layer <- 1
  
  for (month in 1:12) {
    end_layer <- start_layer + days_per_month[month] - 1
    monthly_sums[[month]] <- app(daily_raster[[start_layer:end_layer]], sum)
    start_layer <- end_layer + 1
  }
  
  yearly_raster <- rast(monthly_sums)
  
  # Speichere die Monatsdaten direkt als NetCDF für jedes Jahr
  output_file_year <- paste0(data_path, year, "_monthly.nc")
  writeCDF(yearly_raster, output_file_year, varname="precip", overwrite=TRUE)
  
  cat("Gespeichert:", output_file_year, "\n")  # Bestätigung ausgeben
  
  # Speicher freigeben
  rm(daily_raster, yearly_raster, monthly_sums)
  gc()
}

cat("Alle Jahre erfolgreich verarbeitet!\n")

############################################################################
# Setze den Pfad zu den NetCDF-Dateien
data_path <- "C:/Users/Jonas/Desktop/StatistikSeminar/"

# Liste der Jahre
years <- 1950:1990

# Liste für alle Monatsdaten
all_years_raster <- list()

# Alle NetCDF-Dateien einladen und zusammenfügen
for (year in years) {
  file_path <- paste0(data_path, year, "_monthly.nc")
  
  if (!file.exists(file_path)) {
    warning(paste("Datei fehlt:", file_path))
    next  # Falls eine Datei fehlt, einfach weitermachen
  }
  
  cat("Lade Datei:", file_path, "\n")
  yearly_raster <- rast(file_path)
  
  # Füge das Raster zur Liste hinzu
  all_years_raster <- c(all_years_raster, yearly_raster)
}

# Erstelle ein SpatRaster aus der Liste
final_raster <- rast(all_years_raster)

# Speichere das gesamte Raster als eine einzige NetCDF-Datei
output_file <- paste0(data_path, "1950_1990_monthly.nc")
writeCDF(final_raster, output_file, varname="precip", overwrite=TRUE)

cat("Gesamtdaten erfolgreich gespeichert in:", output_file, "\n")

data_path <- "C:/Users/Jonas/Desktop/StatistikSeminar/"
final_file <- paste0(data_path, "1950_1990_monthly.nc")

# Prüfen, ob die Datei existiert
if (!file.exists(final_file)) {
  stop("Die Datei existiert nicht: ", final_file)
}

# NetCDF-Datei laden
final_raster <- rast(final_file)

# Datei-Infos ausgeben
print(final_raster)

# Test: Ein Monat aus einem bestimmten Jahr plotten (z. B. August 1980)
year <- 1980
month <- 8
layer_index <- ((year - 1950) * 12) + month  # Richtigen Layer berechnen

cat("Plotte Jahr:", year, "Monat:", month, "Layer:", layer_index, "\n")

# Plot des ausgewählten Monats
plot(final_raster[[layer_index]], main=paste("Niederschlag:", year, "-", month))


# in mm umrechnen


input_file <- "C:/Users/Jonas/Desktop/StatistikSeminar/1950_1990_monthly.nc"

# Datei einladen
monthly_raster <- rast(input_file)

# Umrechnung von Metern in Millimeter (Faktor 1000)
monthly_raster_mm <- monthly_raster * 1000

# **Speichern der umgerechneten Datei**
output_file_mm <- "C:/Users/Jonas/Desktop/StatistikSeminar/1950_1990_monthly_mm.nc"
writeCDF(monthly_raster_mm, output_file_mm, varname="pr_sum")

# **Überprüfung der Ausgabe**
print(monthly_raster_mm)
plot(monthly_raster_mm[[1]], main = "Monatliche Niederschlagssumme - Januar 1950 (mm)")
plot(monthly_raster_mm[[492]], main = "Monatliche Niederschlagssumme - Dezember 1990 (mm)")

file_path <- "C:/Users/Jonas/Desktop/StatistikSeminar/1950_1990_monthly_mm.nc"

# NetCDF-Datei als Raster einladen
pr_data <- rast(file_path)

# Überprüfen, ob die Daten geladen wurden
print(pr_data)
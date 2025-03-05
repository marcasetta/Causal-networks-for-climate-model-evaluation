library(terra)
library(ncdf4)

# Dateipfad zur NetCDF-Datei
file_path <- "C:/Users/Jonas/Desktop/StatistikSeminar/pr_day_ACCESS-ESM1-5_historical_r16i1p1f1_gn_19500101-19991231.nc"

#NetCDF-Datei als RasterStack einladen
pr_data <- rast(file_path)

#Überblick über die Daten
print(pr_data)

# Metadaten auslesen (Variable, Dimensionen, Zeit)
nc <- nc_open(file_path)  # Datei öffnen
print(nc)                 # Alle Metadaten anzeigen
nc_close(nc)              # Datei wieder schließen

# Einzelne Tage plotten (z. B. 1. Januar 1950 und 1000. Tag)
plot(pr_data[[1]], main = "Niederschlag - 1. Januar 1950")  
plot(pr_data[[1000]], main = "Niederschlag - Tag 1000")  

# Zeitkoordinaten auslesen
time_values <- time(pr_data)  
print(time_values)  # Zeigt die Zeitstempel

# Falls die Zeit in "Tagen seit ..." gespeichert ist, Umrechnung in Datum
# time_converted <- as.Date(time_values, origin="1950-01-01")
# print(time_converted)

file_path <- "C:/Users/Jonas/Desktop/StatistikSeminar/pr_day_ACCESS-ESM1-5_historical_r16i1p1f1_gn_19500101-19991231.nc"

# NetCDF-Datei als Raster einladen
pr_data <- rast(file_path)

# Überprüfen, ob die Daten geladen wurden
print(pr_data)

# Funktion zur Bestimmung der Tage pro Monat (Schaltjahr beachten)
days_in_month <- function(year) {
  if ((year %% 4 == 0 && year %% 100 != 0) || (year %% 400 == 0)) {
    return(c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))  # Schaltjahr
  } else {
    return(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))  # Normales Jahr
  }
}

# Start- und Endjahr aus Dateiname entnehmen (1950–1999)
start_year <- 1950
end_year <- 1999

# Liste zur Speicherung der monatlichen Summen
monthly_sums <- list()

# Index für die Layer
start_layer <- 1

# Loop über alle Jahre und Monate
for (year in start_year:end_year) {
  days_per_month <- days_in_month(year)
  
  for (month in 1:12) {
    end_layer <- start_layer + days_per_month[month] - 1
    
    # Umrechnung: pr (kg/m²/s) * Sekunden pro Tag * Tage pro Monat
    monthly_sums[[paste(year, month, sep = "_")]] <- app(pr_data[[start_layer:end_layer]], sum) * 86400
    
    # Nächster Startlayer
    start_layer <- end_layer + 1
  }
}

# In SpatRaster umwandeln
monthly_raster <- rast(monthly_sums)

# Plot eines Monats als Beispiel (Januar 1950)
plot(monthly_raster[[1]], main = "Monatliche Niederschlagssumme - Januar 1950")

# Speichern der monatlichen Summen als NetCDF
output_file <- "C:/Users/Jonas/Desktop/StatistikSeminar/pr_monthly_sums_ACCESS-ESM1-5.nc"
writeCDF(monthly_raster, output_file, varname="pr_sum")

# Kontrolle der Ausgabe
print(monthly_raster)

output_file <- "C:/Users/Jonas/Desktop/StatistikSeminar/pr_monthly_sums_ACCESS-ESM1-5.nc"

# NetCDF-Datei als Raster einladen
monthly_raster <- rast(output_file)

# Überblick über die Daten
print(monthly_raster)

# **Plot der monatlichen Niederschlagssummen für verschiedene Monate**
# Januar 1950 (erste Schicht)
plot(monthly_raster[[1]], main = "Monatliche Niederschlagssumme - Januar 1950")

# Juli 1950 (7. Schicht)
plot(monthly_raster[[7]], main = "Monatliche Niederschlagssumme - Juli 1950")

# Dezember 1999 (letzte Schicht)
plot(monthly_raster[[nlyr(monthly_raster)]], main = "Monatliche Niederschlagssumme - Dezember 1999")

#auf 1950-1990 zuschneiden 

# Original-Datei mit allen Jahren
input_file <- "C:/Users/Jonas/Desktop/StatistikSeminar/pr_monthly_sums_ACCESS-ESM1-5.nc"

# NetCDF als Raster laden
monthly_raster <- rast(input_file)

# Anzahl der Monate insgesamt
num_months <- nlyr(monthly_raster)
print(paste("Gesamtzahl der Monate:", num_months))

# **Schichten auswählen:**
# 1950–1990 = 41 Jahre * 12 Monate = 492 Schichten (von 1 bis 492)
filtered_raster <- monthly_raster[[1:492]]

# **Neue Datei speichern**
output_file <- "C:/Users/Jonas/Desktop/StatistikSeminar/pr_monthly_sums_ACCESS-ESM1-5_1950-1990.nc"
writeCDF(filtered_raster, output_file, varname="pr_sum")

# **Ergebnis überprüfen**
print(filtered_raster)
plot(filtered_raster[[1]], main = "Monatliche Niederschlagssumme - Januar 1950")
plot(filtered_raster[[492]], main = "Monatliche Niederschlagssumme - Dezember 1990")
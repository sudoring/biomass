# Load necessary libraries
library(sf)

# Define the center point (lat, lon)
lat_center <- 37.460749
lon_center <- 128.516865

# Create a data frame for the center point
center_point <- data.frame(lon = lon_center, lat = lat_center)

# Convert the data frame to an sf object
center_sf <- st_as_sf(center_point, coords = c("lon", "lat"), crs = 4326)

# Reproject the center point to a CRS that uses meters (UTM zone 52N for South Korea)
center_utm <- st_transform(center_sf, crs = 32652)

# Define half of the side of the square in meters (11 km / 2 = 5.5 km)
half_side_meters <- 5500

# Create a bounding box (square) in UTM coordinates
utm_square <- st_buffer(center_utm, dist = half_side_meters)

# Convert the UTM square back to lat/lon (CRS 4326)
square_latlon <- st_transform(utm_square, crs = 4326)

# Write the GeoJSON file using the sf package's st_write function
st_write(square_latlon, "11km_square.geojson", driver = "GeoJSON")

# Print a success message
cat("GeoJSON file '11km_square.geojson' created successfully.")

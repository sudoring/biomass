# Load necessary libraries
library(sp)
library(terra)  # For projection transformations
library(jsonlite)  # For writing GeoJSON

# Site name
siteName <- 'gariwang'

# Define the center point (lat, lon)
lat_center <- 37.460749
lon_center <- 128.516865

# Create a matrix for the center point in WGS84 (EPSG:4326)
center_coords <- cbind(lon_center, lat_center)

# Create a SpatVector (terra) object for the center point in WGS84
center_point <- vect(center_coords, crs = "+proj=longlat +datum=WGS84")

# Reproject the center point to UTM zone 52N (EPSG:32652)
center_point_utm <- project(center_point, "+proj=utm +zone=52 +datum=WGS84")

# Define half of the side of the square in meters (11 km / 2 = 5.5 km)
half_side_meters <- 5500

# Extract the UTM coordinates of the center point
utm_x <- crds(center_point_utm)[1, 1]
utm_y <- crds(center_point_utm)[1, 2]

# Calculate the UTM coordinates for the square corners
square_coords_utm <- matrix(c(
  utm_x - half_side_meters, utm_y - half_side_meters,
  utm_x + half_side_meters, utm_y - half_side_meters,
  utm_x + half_side_meters, utm_y + half_side_meters,
  utm_x - half_side_meters, utm_y + half_side_meters,
  utm_x - half_side_meters, utm_y - half_side_meters
), ncol = 2, byrow = TRUE)

# Create a polygon from the UTM square coordinates
square_polygon_utm <- vect(list(square_coords_utm), type = "polygons", crs = "+proj=utm +zone=52 +datum=WGS84")

# Reproject the square back to WGS84 (EPSG:4326)
square_polygon_wgs84 <- project(square_polygon_utm, "+proj=longlat +datum=WGS84")

# Extract the coordinates in WGS84 for GeoJSON
square_coords_wgs84 <- crds(square_polygon_wgs84)

# Create GeoJSON structure manually
geojson_list <- list(
  type = "FeatureCollection",
  name = siteName,
  features = list(
    list(
      type = "Feature",
      properties = list('f',siteName),  # Empty properties
      geometry = list(
        type = "Polygon",
        coordinates = list(square_coords_wgs84)
      )
    )
  )
)

# Write the GeoJSON to a file using jsonlite
path <- '/projectnb/modislc/users/mkmoon/biomass/geojson/'
write_json(geojson_list, paste0(path,siteName,'.geojson'), auto_unbox = TRUE)

# comments
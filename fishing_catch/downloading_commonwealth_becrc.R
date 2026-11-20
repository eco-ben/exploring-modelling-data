library(sf)
url <- "https://www.cmar.csiro.au/geoserver/becrc/ows?service=WFS&version=2.0.0&request=GetFeature&typeNames=becrc:LAYERNAME&outputFormat=application/json"
x <- st_read(url)

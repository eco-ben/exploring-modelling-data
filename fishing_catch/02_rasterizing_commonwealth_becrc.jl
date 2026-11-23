import GeoDataFrames as GDF
import GeometryOps as GO
import GeoInterface as GI

using Rasters
using NetCDF
using YAXArrays
using DimensionalData

using Statistics
using DataFrames
using GLMakie

data_dir = "C:/Users/grier/OneDrive - The University of Queensland/Documents/Data/Commonwealth_Spatial_Fisheries_Catch/"

files = readdir(data_dir)
files = files[contains.(files, [".gpkg"])]

test_polygons = GDF.read(joinpath(data_dir, files[1]))
test_polygons = test_polygons[test_polygons.catch_kg_total .> 0, :]

poly_fig = poly(test_polygons.geometry, color=test_polygons.catch_kg_total)
save("first_test_polygons.png", poly_fig)

test_raster = Rasters.rasterize(sum, test_polygons; fill=:catch_kg_total, res=0.1)
missed_geoms = Rasters.metadata(test_raster)["missed_geometries"]

missed_geoms = test_polygons[missed_geoms, :]
missed_points = DataFrame(geometry = GI.Point.(GO.centroid.(missed_geoms.geometry)), catch_kg_total = missed_geoms.catch_kg_total)
GDF.setcrs!(missed_points, crs(test_polygons))
missed_raster = Rasters.rasterize(sum, missed_points; fill=:catch_kg_total, to=test_raster)

# Combine the data from original geometries and from the missing geometries
total_raster = replace_missing(test_raster, 0) .+ replace_missing(missed_raster, 0)

total_raster = ifelse.(total_raster .== 0.0, missing, total_raster)
total_raster = rebuild(total_raster, missingval=-9999.0f0)
total_raster = reverse(total_raster; dims=:Y)

write(joinpath(data_dir, replace(files[1], ".gpkg" => ".tif")), total_raster, force=true)
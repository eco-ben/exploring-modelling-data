import GeoDataFrames as GDF
import GeoInterface as GI
import GeometryOps as GO
using GLMakie
import ArchGDAL as AG

using DataFrames

wfs_url = "WFS:https://www.cmar.csiro.au/geoserver/becrc/wfs"
dataset = AG.read(wfs_url)

target_data_dir = "C:/Users/uqbgrier/OneDrive - The University of Queensland/Documents/Data/Commonwealth_Spatial_Fisheries_Catch/"

nlayers = AG.nlayer(dataset)
layer_names = []
for i in 1:nlayers
    try
        push!(layer_names, AG.getname(AG.getlayer(dataset, i)))    
    catch
        print("Error getting layer number: $(i)")
    end
end

fisheries_layers = layer_names[contains.(layer_names, ["fisher"])]
errored_layers = []

for i_layer in 1:length(fisheries_layers)
    layer_name = fisheries_layers[i_layer]
    filename = target_data_dir * "/" * last(split(layer_name, ":")) * ".gpkg"
    if isfile(filename) continue end
    try
        layer_df = GDF.read(wfs_url, layer = layer_name)

        layer_df.geometry = AG.lineargeom.(layer_df.the_geom)
        layer_df = layer_df[layer_df.catch_kg .> -1, :]

        if nrow(layer_df) == 0 continue end
        
        select!(layer_df, Not(:the_geom))
        GDF.setgeometrycolumn!(layer_df, :geometry)
        
        GDF.write(filename, layer_df)
    catch
        print("Error loading or converting fisheries_layer: $(i_layer)")
        push!(errored_layers, i_layer)
    end
end

csv_files = readdir(target_data_dir)
csv_files = csv_files[contains.(csv_files, ".csv")]

base_df = GDF.read(wfs_url, layer=1)
base_df.geometry = AG.lineargeom.(base_df.the_geom)

for file in csv_files
    print("Loading file $(file) \n")

    out_filename = target_data_dir * "/" * replace(file, ".csv" => ".gpkg")

    layer_df = GDF.read(
        GDF.ArchGDALDriver(),
        target_data_dir * "/" * file,
        options = ["GEOM_POSSIBLE_NAMES=the_geom", "KEEP_GEOM_COLUMNS=NO"]
    )

    layer_df.catch_kg_total = parse.([Float64], layer_df.catch_kg_total)
    layer_df = layer_df[layer_df.catch_kg_total .> -1, :]

    
    layer_df.geometry = AG.lineargeom.(layer_df.the_geom)
    layer_df = layer_df[.!GI.isempty.(layer_df.geometry), :]
    layer_df.geometry = GO.flip(layer_df.geometry)

    if nrow(layer_df) == 0 continue end
    
    select!(layer_df, Not(:the_geom))
    GDF.setgeometrycolumn!(layer_df, :geometry)
    
    GDF.write(filename, layer_df[:, [:geometry, :time_period, :catch_kg_total]])
end
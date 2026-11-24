import GeoDataFrames as GDF
import GeometryOps as GO
import GeoInterface as GI
import DimensionalData as DD

using Rasters
using NetCDF
using YAXArrays
using DimensionalData

using Statistics
using DataFrames
using GLMakie

data_dir = "C:/Users/uqbgrier/OneDrive - The University of Queensland/Documents/Data/Commonwealth_Spatial_Fisheries_Catch/"

files = readdir(data_dir)
files = files[contains.(files, [".gpkg"])]

timespans = Dict(
    "200304_200708" => "200304-200708",
    "200809_201213" => "200809-201213",
    "201314_201718" => "201314-201718",
    "201819_202223" => "201819-202223",
    "2022_2023" => "2022-2023",
    "20222023" => "2022-2023",
    "202223_geom_20222023" => "2022-2023"
)
fisheries = [
    "alldata_allfisheries",
    "afma_Bass_Strait_Central_Zone_Scallop",
    "afma_Coral_Sea",
    "afma_North_West_Slope_Trawl",
    "afma_Northern_Prawn",
    "afma_Southern_and_Eastern_Scalefish_and_Shark",
    "afma_Southern_Squid_Jig",
    "afma_Torres_Strait_Finfish",
    "afma_Torres_Strait_Prawn",
    "afma_Torres_Strait_Rock_Lobster",
    "afma_Tuna_and_Billfish",
    "afma_Western_Deepwater_Trawl",
    "nsw_Abalone",
    "nsw_Estuary_General",
    "nsw_Estuary_Prawn_Trawl",
    "nsw_Lobster",
    "nsw_Ocean_Hauling",
    "nsw_Ocean_Trawl",
    "nsw_s37_Permit",
    "nsw_Southern_Fish_Trawl",
    "qld_East_Coast_Inshore_Finfish",
    "qld_East_Coast_Otter_Trawl",
    "qld_Gulf_Inshore_Finfish",
    "qld_Gulf_Line",
    "qld_Other",
    "qld_Queensland_Crab",
    "qld_Reef_Line",
    "qld_River__Inshore_Beam_Trawl",
    "qld_Tropical_Rock_Lobster",
    "sa_abalone",
    "sa_giant_crab",
    "sa_goolwa_cockle",
    "sa_lakes_and_coorong",
    "sa_lcf_prod",
    "sa_marine_scale_prod",
    "sa_msf_prod",
    "sa_prawn",
    "sa_river_prod",
    "sa_rock_lobster",
    "sa_sardine",
    "tas_abalone",
    "tas_giant_crab",
    "tas_scalefish",
    "tas_scallop",
    "tas_southern_rock_lobster",
    "wa_Abalone",
    "wa_Abrolhos_Islands_and_Mid_West_Trawl",
    "wa_Broome_Prawn",
    "wa_c127_and_L02_Joint_Authority_Northern_Shark",
    "wa_Cockburn_Sound_Crab",
    "wa_Cockburn_Sound_Line_and_Pot",
    "wa_Exmouth_Gulf_Developing_Crab",
    "wa_Exmouth_Gulf_Prawn",
    "wa_FBL_condition_19_Estuaries",
    "wa_FBL_condition_42_herring",
    "wa_FBL_condition_73_South_Coast_Trawl",
    "wa_FBL_condition_74_Fish_Trapping",
    "wa_FBL_condition_84_SW_Coast_Estuaries",
    "wa_FBL_condition_93_Purse_Seine_Development_Zone",
    "wa_FBL_condition_192_King_George_Sound_Leatherjackets",
    "wa_Gascoyne_Demersal_Scalefish",
    "wa_Joint_Authority_Northern_Shark",
    "wa_Kimberley_Crab",
    "wa_Kimberley_Gillnet_and_Barramundi",
    "wa_Kimberley_Prawn",
    "wa_Mackerel",
    "wa_Mandurah_to_Bunbury_Developing_Crab",
    "wa_Nickol_Bay_Prawn",
    "wa_Northern_Demersal_Scalefish",
    "wa_Octopus_Interim",
    "wa_Onslow_Prawn",
    "wa_Open_access",
    "wa_Open_Access_in_the_North_Coast_Gascoyne_Coast_and_West_Coast_Bioregions",
    "wa_Open_Access_in_the_South_Coast_Bioregion",
    "wa_Pearl_Oyster",
    "wa_Pilbara_Crab",
    "wa_Pilbara_Fish_Trawl_Interim",
    "wa_Pilbara_Line",
    "wa_Pilbara_Trap",
    "wa_Shark_Bay_Beach_Seine_and_Mesh_Net",
    "wa_Shark_Bay_Crab",
    "wa_Shark_Bay_Pink_Snapper",
    "wa_Shark_Bay_Prawn",
    "wa_Shark_Bay_Scallop",
    "wa_South_Coast_Crustacean",
    "wa_South_Coast_Estuarine",
    "wa_South_Coast_Purse-Seine",
    "wa_South_Coast_Salmon",
    "wa_South_West_Coast_Beach_Net",
    "wa_South_West_Coast_Salmon",
    "wa_South_West_Trawl",
    "wa_Southern_Demersal_Gillnet_and_Demersal_Longline",
    "wa_West_Australian_Sea_Cucumber",
    "wa_West_Coast_Deep_Sea_Crustacean",
    "wa_West_Coast_Demersal_Gillnet_and_Demersal_Longline",
    "wa_West_Coast_Demersal_Scalefish",
    "wa_West_Coast_Estuarine",
    "wa_West_Coast_Purse_Seine",
    "wa_West_Coast_Rock_Lobster",
    "wa_Western_Australian_North_Coast_Shark"
]

template_polygons = GDF.read(joinpath(data_dir, files[1]))
template_polygons = template_polygons[template_polygons.catch_kg_total .> 0, :]

template_raster = Rasters.rasterize(
    sum, 
    template_polygons; 
    fill=:catch_kg_total, 
    res=0.1, 
    crs=crs(template_polygons)
)

fishery = fisheries[1]

fishery_files = files[contains.(files, [fishery])]

fishery_timespans = sort([
    collect(keys(timespans))[findfirst(contains.([f], keys(timespans)))] 
    for f in fishery_files
])

function yaxarray_from_gpkg(fishery_files, input_timespan; data_dir=data_dir, timespans=timespans, template_raster=template_raster, catch_col=:catch_kg)
    input_fn = first(fishery_files[contains.(fishery_files, input_timespan)])

    polygons = GDF.read(joinpath(data_dir, input_fn))
    polygons = polygons[polygons[:, catch_col] .> 0, :]

    base_raster = Rasters.rasterize(
        sum,
        polygons;
        fill=catch_col,
        to=template_raster
    )
    
    if !isempty(Rasters.metadata(base_raster))
        missed_geoms = Rasters.metadata(base_raster)["missed_geometries"]
        if any(missed_geoms)
            missed_polygons = polygons[missed_geoms, :]
            missed_points = DataFrame(
                geometry = GI.Point.(GO.centroid.(missed_polygons.geometry))
            )
            missed_points[!, catch_col] = missed_polygons[:, catch_col]
            GDF.setcrs!(missed_points, crs(polygons))

            missed_raster = Rasters.rasterize(
                sum, 
                missed_points; 
                fill=catch_col, 
                to=template_raster
            )

            if !isempty(Rasters.metadata(missed_raster)) error("Rasterized missed_points still containing missed geometries.") end
            
            base_raster = replace_missing(base_raster, 0) .+ 
                replace_missing(missed_raster, 0)
        end

        base_raster = rebuild(
                base_raster, 
                metadata=delete!(
                    DD.val(Rasters.metadata(base_raster)), 
                    "missed_geometries"
            )
        )
        base_raster = replace_missing(base_raster, 0)
    end

    if isa(base_raster.data, Matrix{Int64})
        base_raster = rebuild(base_raster, data=Float64.(base_raster.data))
    end

    base_raster = reverse(base_raster; dims=:Y)
    base_raster = base_raster ./ 1000 # Convert from kg catch to tonnes
    base_raster = ifelse.(base_raster .== 0.0, -9999.0f0, base_raster)
    base_raster = rebuild(base_raster, missingval=-9999.0f0)

    base_raster = rebuild(
        base_raster,
        metadata=merge(
            Rasters.metadata(base_raster),
            Dict("Units" => "tonnes")
        )
    )
    base_array = YAXArray(
        dims(base_raster), 
        Float64.(parent(base_raster)), 
        Rasters.metadata(base_raster)
    )
    # base_array = rebuild(
    #     base_array,
    #     data = reshape(base_array.data, size(base_array)..., 1),
    #     dims = (dims(base_array)..., Dim{:time}([output_timespan]))
    # )

    return base_array
end

for fishery in fisheries[fisheries .!= "alldata_allfisheries"]
    print("Converting $(fishery) \n")

    fishery_files = files[contains.(files, [fishery])]

    fishery_timespans = sort([
        collect(keys(timespans))[findfirst(contains.([f], keys(timespans)))] 
        for f in fishery_files
    ])
    
    if isfile(joinpath(data_dir, "$(fishery).nc")) continue end

    fishery_arrays = [yaxarray_from_gpkg(fishery_files, t) for t in fishery_timespans]
    fishery_array = concatenatecubes(fishery_arrays, Ti(fishery_timespans))
    savecube(fishery_array, joinpath(data_dir, "$(fishery).nc"); overwrite=true)
end

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

total_raster = total_raster ./ 1000

write(joinpath(data_dir, replace(files[1], ".gpkg" => ".tif")), total_raster, force=true)
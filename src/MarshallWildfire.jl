module MarshallWildfire

using ArchGDAL, CSV, DataFrames, Dates, Downloads, Extents, GeoJSON, GeoSurrogates, JSON3, Tyler,
    NCDatasets, ModelingToolkit, GLMakie, GeoMakie, OSMGeocoder, ProgressMeter, URIs,
    RapidRefreshData, Rasters, Statistics, TimeZones, Tyler, TileProviders
using ModelingToolkit: t_nounits as t, D_nounits as D
using GeoJSON: GeoJSON

import Landfire, WebAssets, Wildfires

import GeoInterface as GI
import GeoFormatTypes as GFT
import GeometryOps as GO

#-----------------------------------------------------------------------------# constants
const extent = Extent(X = (-105.23332f0, -105.131424f0), Y = (39.92924f0, 39.98638f0))
const extent2 = Extents.grow(extent, 0.1f0)

const ignition_point = (lon = -105.231f0, lat = 39.955f0)

const unique_fire_id = "2021-COBLX-000995"
const irwin_id = "C63FC371-BC70-4615-841B-B0838C21064F"

const start_time = ZonedDateTime(DateTime(2021, 12, 30, 10), tz"America/Denver")
const stop_time = ZonedDateTime(DateTime(2022, 1, 1, 6), tz"America/Denver")
const start_time_utc = DateTime(start_time, UTC)
const stop_time_utc = DateTime(stop_time, UTC)


#-----------------------------------------------------------------------------#
include("data.jl")
include("surrogates.jl")
include("rothermel.jl")
include("level_set.jl")
include("model.jl")
include("plots.jl")
include("map_viewer.jl")
include("wildfires_sim.jl")

#-----------------------------------------------------------------------------# stats
function stats()
    perim = get_perimeter()
    area_km2 = GO.area(GO.Geodesic(), perim.geometry) / 1000 ^ 2
end


end

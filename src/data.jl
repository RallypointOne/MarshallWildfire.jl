#-----------------------------------------------------------------------------# download helper with retry
"""
    download_with_retry(url, path; max_retries=3, base_delay=2.0)

Download a URL with exponential backoff retry logic.
Returns `true` on success, `false` on failure after all retries.
"""
function download_with_retry(url::String, path::String; max_retries::Int=3, base_delay::Float64=2.0)
    for attempt in 1:max_retries
        try
            Downloads.download(url, path)
            return true
        catch e
            if attempt < max_retries
                delay = base_delay * (2 ^ (attempt - 1))
                @warn "Download failed (attempt $attempt/$max_retries), retrying in $(delay)s..." exception=e
                sleep(delay)
            else
                @warn "Download failed after $max_retries attempts" exception=e
                return false
            end
        end
    end
    return false
end

const EMPTY_FEATURE_COLLECTION = """{"type":"FeatureCollection","features":[]}"""

#-----------------------------------------------------------------------------# perimeter
const PERIMETER_URL = "https://services3.arcgis.com/0jWpHMuhmHsukKE3/arcgis/rest/services/bcpos_Marshall_Fire_Perimeter/FeatureServer/0/query?where=1%3D1&outFields=*&returnGeometry=true&outSR=4326&f=geojson"

function get_perimeter()
    path = joinpath(@__DIR__, "..", "data", "perimeter.geojson")
    if !isfile(path)
        mkpath(dirname(path))
        Downloads.download(PERIMETER_URL, path)
    end
    GeoJSON.read(path)
end

#-----------------------------------------------------------------------------# landfire
function get_landfire_data()
    prods = Landfire.products(conus=true)
    data = Landfire.Dataset(prods, extent2; output_projection="4326")
    tif = Landfire.get(data)
    r = replace_missing(Raster(tif))
    (; (Symbol(b) => r[Band = At(b)] for b in lookup(r, Rasters.Band))...)
end

#-----------------------------------------------------------------------------# building footprints
"""
    get_building_footprints()

Download building footprints from OpenStreetMap for the extent2 area.
Uses the Overpass API to query building polygons.
Returns an empty FeatureCollection if the download fails after retries.
"""
function get_building_footprints()
    path = joinpath(@__DIR__, "..", "data", "buildings.geojson")

    if !isfile(path)
        mkpath(dirname(path))

        # Overpass API query for buildings in extent2
        xmin, xmax = extent2.X
        ymin, ymax = extent2.Y
        bbox = "$(ymin),$(xmin),$(ymax),$(xmax)"

        query = """
        [out:json][timeout:120];
        (
          way["building"]($(bbox));
          relation["building"]($(bbox));
        );
        out body;
        >;
        out skel qt;
        """

        # URL encode the query
        encoded_query = URIs.escapeuri(query)
        url = "https://overpass-api.de/api/interpreter?data=$(encoded_query)"

        # Download OSM JSON with retry logic
        osm_path = joinpath(@__DIR__, "..", "data", "buildings_osm.json")
        if !download_with_retry(url, osm_path; max_retries=3, base_delay=5.0)
            @warn "Failed to download building footprints from Overpass API, using empty dataset"
            write(path, EMPTY_FEATURE_COLLECTION)
            return GeoJSON.read(path)
        end

        # Convert to GeoJSON using osmtogeojson or manual parsing
        # For simplicity, use the raw Overpass JSON and parse it
        osm_data = JSON3.read(read(osm_path, String))

        # Build GeoJSON features from OSM data
        features = _osm_to_geojson_buildings(osm_data)

        # Write GeoJSON
        geojson = Dict(
            "type" => "FeatureCollection",
            "features" => features
        )
        write(path, JSON3.write(geojson))

        # Clean up temp file
        rm(osm_path, force=true)
    end

    GeoJSON.read(path)
end

function _osm_to_geojson_buildings(osm_data)
    # Build node lookup
    nodes = Dict{Int64, Tuple{Float64, Float64}}()
    for element in osm_data.elements
        if element.type == "node"
            nodes[element.id] = (element.lon, element.lat)
        end
    end

    # Build features from ways
    features = []
    for element in osm_data.elements
        if element.type == "way" && hasproperty(element, :nodes)
            coords = [get(nodes, nid, nothing) for nid in element.nodes]
            filter!(!isnothing, coords)

            if length(coords) >= 4
                # Close the polygon if needed
                if coords[1] != coords[end]
                    push!(coords, coords[1])
                end

                feature = Dict(
                    "type" => "Feature",
                    "geometry" => Dict(
                        "type" => "Polygon",
                        "coordinates" => [coords]
                    ),
                    "properties" => hasproperty(element, :tags) ? Dict(pairs(element.tags)) : Dict()
                )
                push!(features, feature)
            end
        end
    end

    return features
end

#-----------------------------------------------------------------------------# power lines
"""
    get_power_lines()

Download power line locations from OpenStreetMap for the extent2 area.
Uses the Overpass API to query power lines and minor power lines.
Returns an empty FeatureCollection if the download fails after retries.
"""
function get_power_lines()
    path = joinpath(@__DIR__, "..", "data", "powerlines.geojson")

    if !isfile(path)
        mkpath(dirname(path))

        # Overpass API query for power lines in extent2
        xmin, xmax = extent2.X
        ymin, ymax = extent2.Y
        bbox = "$(ymin),$(xmin),$(ymax),$(xmax)"

        query = """
        [out:json][timeout:120];
        (
          way["power"="line"]($(bbox));
          way["power"="minor_line"]($(bbox));
          way["power"="cable"]($(bbox));
        );
        out body;
        >;
        out skel qt;
        """

        # URL encode the query
        encoded_query = URIs.escapeuri(query)
        url = "https://overpass-api.de/api/interpreter?data=$(encoded_query)"

        # Download OSM JSON with retry logic
        osm_path = joinpath(@__DIR__, "..", "data", "powerlines_osm.json")
        if !download_with_retry(url, osm_path; max_retries=3, base_delay=5.0)
            @warn "Failed to download power lines from Overpass API, using empty dataset"
            write(path, EMPTY_FEATURE_COLLECTION)
            return GeoJSON.read(path)
        end

        osm_data = JSON3.read(read(osm_path, String))

        # Build GeoJSON features from OSM data
        features = _osm_to_geojson_lines(osm_data)

        # Write GeoJSON
        geojson = Dict(
            "type" => "FeatureCollection",
            "features" => features
        )
        write(path, JSON3.write(geojson))

        # Clean up temp file
        rm(osm_path, force=true)
    end

    GeoJSON.read(path)
end

function _osm_to_geojson_lines(osm_data)
    # Build node lookup
    nodes = Dict{Int64, Tuple{Float64, Float64}}()
    for element in osm_data.elements
        if element.type == "node"
            nodes[element.id] = (element.lon, element.lat)
        end
    end

    # Build features from ways
    features = []
    for element in osm_data.elements
        if element.type == "way" && hasproperty(element, :nodes)
            coords = [get(nodes, nid, nothing) for nid in element.nodes]
            filter!(!isnothing, coords)

            if length(coords) >= 2
                feature = Dict(
                    "type" => "Feature",
                    "geometry" => Dict(
                        "type" => "LineString",
                        "coordinates" => coords
                    ),
                    "properties" => hasproperty(element, :tags) ? Dict(pairs(element.tags)) : Dict()
                )
                push!(features, feature)
            end
        end
    end

    return features
end

#-----------------------------------------------------------------------------# Damage Assessment Data
# From Boulder County's official damage assessment web app:
# https://www.arcgis.com/apps/webappviewer/index.html?id=9f3314c39ad64fac925101aae0bdd62c

const RDA_URL = "https://services3.arcgis.com/0jWpHMuhmHsukKE3/arcgis/rest/services/RDA_Public_View/FeatureServer/0/query?where=1%3D1&outFields=*&returnGeometry=true&outSR=4326&f=geojson"
const CDA_URL = "https://services3.arcgis.com/0jWpHMuhmHsukKE3/arcgis/rest/services/CDA_Public_View/FeatureServer/0/query?where=1%3D1&outFields=*&returnGeometry=true&outSR=4326&f=geojson"

"""
    get_residential_damage_assessment()

Download residential damage assessment points from Boulder County's official damage assessment.
Returns a GeoJSON FeatureCollection with point geometries.

# Attributes:
- `addr`: Address
- `damage`: Damage category ("Destroyed", "Major", "Minor", "Affected")
- `jurisdiction`: City/jurisdiction
"""
function get_residential_damage_assessment()
    path = joinpath(@__DIR__, "..", "data", "residential_damage_assessment.geojson")
    if !isfile(path)
        mkpath(dirname(path))
        Downloads.download(RDA_URL, path)
    end
    GeoJSON.read(path)
end

"""
    get_commercial_damage_assessment()

Download commercial damage assessment points from Boulder County's official damage assessment.
Returns a GeoJSON FeatureCollection with point geometries.

# Attributes:
- `addr`: Address
- `damage`: Damage category ("Destroyed", "Major", "Minor", "Affected")
- `jurisdiction`: City/jurisdiction
"""
function get_commercial_damage_assessment()
    path = joinpath(@__DIR__, "..", "data", "commercial_damage_assessment.geojson")
    if !isfile(path)
        mkpath(dirname(path))
        Downloads.download(CDA_URL, path)
    end
    GeoJSON.read(path)
end

#-----------------------------------------------------------------------------# NREL WTK-LED wind data
"""
    get_nrel_wind_data(; year=2020)

Download 5-minute wind data from the NREL Wind Toolkit Long-term Ensemble Dataset (WTK-LED)
for the Marshall Fire ignition point.

Returns attributes: windspeed_10m, winddirection_10m, turbulent_kinetic_energy_20m.

Requires `ENV["NREL_API_KEY"]` to be set. Get a free key at https://developer.nrel.gov/signup/

The data is cached at `data/nrel/wtk_led_5min_{year}.csv`.
"""
function get_nrel_wind_data(; year::Int=2020)
    path = joinpath(@__DIR__, "..", "data", "nrel", "wtk_led_5min_$(year).csv")

    if !isfile(path)
        mkpath(dirname(path))

        api_key = get(ENV, "NREL_API_KEY", "")
        isempty(api_key) && error("Set ENV[\"NREL_API_KEY\"] to download NREL data. Get a free key at https://developer.nrel.gov/signup/")

        lon, lat = ignition_point.lon, ignition_point.lat
        url = string(
            "https://developer.nrel.gov/api/wind-toolkit/v2/wind/wtk-conus-5min-v1-0-0-download.csv",
            "?api_key=", api_key,
            "&wkt=POINT($(lon)%20$(lat))",
            "&attributes=windspeed_10m,winddirection_10m,turbulent_kinetic_energy_20m",
            "&names=", year,
            "&interval=5",
            "&utc=true",
            "&email=emailjoshday@gmail.com",
        )

        @info "Downloading NREL WTK-LED 5-min wind data for $year..."
        Downloads.download(url, path)
        @info "Cached NREL data at $path"
    end

    return path
end

#-----------------------------------------------------------------------------# HRRR wind data
function get_hrrr_data()
    start_date = Date(start_time_utc)
    start_hour = hour(start_time_utc)
    stop_date = Date(stop_time_utc)
    stop_hour = hour(stop_time_utc)

    datasets = [
        HRRRDataset(date = d, cycle = lpad(h, 2, '0'), forecast = "f00")
        for d in start_date:Day(1):stop_date
        for h in (d == start_date ? start_hour : 0):(d == stop_date ? stop_hour : 23)
    ]

    wind_bands = filter(RapidRefreshData.bands(datasets[1])) do band
        band.variable in ["UGRD", "VGRD"] && band.level == "10 m above ground"
    end

    rasters = map(datasets) do data
        path = get(data, wind_bands)
        r = Raster(path)
        r = resample(r, crs = GFT.EPSG(4326))
        r = crop(r; to = Extents.grow(extent, 3.0f0))
        r = modify(x -> Float64.(x), r)  # hack to remove Missing
        # Rename bands from UGRD/VGRD to u/v
        set(r, Rasters.Band => [:u, :v])
    end

    # Create time dimension from dataset dates/hours
    times = [DateTime(d.date) + Hour(parse(Int, d.cycle)) for d in datasets]

    cat(rasters...; dims=Ti(times))
end

#-----------------------------------------------------------------------------# gridMET fuel moisture data
const GRIDMET_BASE_URL = "http://www.northwestknowledge.net/metdata/data"

"""
    get_gridmet_fuel_moisture(; year=2021, variable=:fm100)

Download gridMET dead fuel moisture data for the specified year and variable.
Data is cropped to the Marshall Fire extent.

# Arguments
- `year`: Year to download (default: 2021 for Marshall Fire)
- `variable`: Fuel moisture variable - `:fm100` (100-hour) or `:fm1000` (1000-hour)

# Returns
A Raster with daily fuel moisture values (%) cropped to the fire extent.

# Data Source
gridMET is a ~4km daily gridded dataset covering CONUS from 1979-present.
See: https://www.climatologylab.org/gridmet.html
"""
function get_gridmet_fuel_moisture(; year::Int=2021, variable::Symbol=:fm100)
    @assert variable in (:fm100, :fm1000) "Variable must be :fm100 or :fm1000"

    filename = "$(variable)_$(year).nc"
    path = joinpath(@__DIR__, "..", "data", "gridmet", filename)

    if !isfile(path)
        mkpath(dirname(path))
        url = "$(GRIDMET_BASE_URL)/$(filename)"
        @info "Downloading gridMET $variable data for $year..."
        Downloads.download(url, path)
    end

    # Load and crop to fire extent (with buffer for context)
    # Use 1.0 degree buffer to show multiple gridMET pixels (~4km resolution)
    r = Raster(path; checkmem=false)
    r = crop(r; to=Extents.grow(extent, 1.0f0))

    return r
end

"""
    get_fuel_moisture_at_ignition()

Get the 100-hour fuel moisture at the Marshall Fire ignition point on December 30, 2021.
Returns fuel moisture as a percentage.
"""
function get_fuel_moisture_at_ignition()
    fm = get_gridmet_fuel_moisture(; year=2021, variable=:fm100)

    # December 30, 2021 is day 364 of the year
    fire_date = Date(2021, 12, 30)

    # Extract value at ignition point for fire date
    fm_day = fm[Ti=At(fire_date)]
    val = fm_day[X=Near(ignition_point.lon), Y=Near(ignition_point.lat)]

    return val
end

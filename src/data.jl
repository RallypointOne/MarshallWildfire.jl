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

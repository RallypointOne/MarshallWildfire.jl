function get_elevation1()
    data = get_landfire_data()
    elev = data.US_ELEV2020
    return GeoSurrogates.RasterWrap(elev)
end

#-----------------------------------------------------------------------------# NRELWindData
"""
    NRELWindData

Parsed NREL WTK-LED 5-minute wind data with precomputed hourly variability statistics.

Provides sub-hourly wind variability metrics (speed std dev, direction std dev, TKE) to
complement the HRRR hourly mean wind used in the simulation.

# Fields
- `times`, `windspeed`, `winddirection`, `tke`: Raw 5-minute data
- `hourly_times`: Hour-aligned timestamps for precomputed stats
- `hourly_speed_std`: Std dev of wind speed within each hour (m/s)
- `hourly_dir_std`: Std dev of wind direction within each hour (degrees)
- `hourly_tke_mean`: Mean TKE within each hour (m²/s²)

# Usage
    nrel = NRELWindData(get_nrel_wind_data())
    speed_std, dir_std, tke = predict(nrel, datetime)
"""
struct NRELWindData
    # Raw 5-minute data
    times::Vector{DateTime}
    windspeed::Vector{Float64}      # m/s at 10m
    winddirection::Vector{Float64}  # degrees
    tke::Vector{Float64}            # m²/s² at 20m

    # Precomputed hourly stats
    hourly_times::Vector{DateTime}
    hourly_speed_std::Vector{Float64}
    hourly_dir_std::Vector{Float64}
    hourly_tke_mean::Vector{Float64}
end

"""
    NRELWindData(csv_path::String)

Parse a WTK-LED CSV file and compute hourly variability statistics.
The CSV has 1 header row of metadata then the actual column headers.
"""
function NRELWindData(csv_path::String)
    # WTK-LED CSVs have 1 metadata row, then column headers, then data
    df = CSV.read(csv_path, DataFrame; header=2, skipto=3)

    # Parse timestamps from Year, Month, Day, Hour, Minute columns
    times = DateTime.(df.Year, df.Month, df.Day, df.Hour, df.Minute)

    # Extract wind columns (CSV column names from NREL API)
    windspeed = Float64.(df[!, "wind speed at 10m (m/s)"])
    winddirection = Float64.(df[!, "wind direction at 10m (deg)"])
    tke = Float64.(df[!, "turbulent kinetic energy at 20m (m^2/s^2)"])

    # Compute hourly statistics
    hourly_times = DateTime[]
    hourly_speed_std = Float64[]
    hourly_dir_std = Float64[]
    hourly_tke_mean = Float64[]

    # Group by hour (floor to hour)
    hour_groups = Dict{DateTime, Vector{Int}}()
    for (i, t) in enumerate(times)
        h = floor(t, Hour)
        push!(get!(Vector{Int}, hour_groups, h), i)
    end

    for h in sort!(collect(keys(hour_groups)))
        idxs = hour_groups[h]
        push!(hourly_times, h)
        push!(hourly_speed_std, length(idxs) > 1 ? std(windspeed[idxs]) : 0.0)
        push!(hourly_dir_std, length(idxs) > 1 ? _circular_std(winddirection[idxs]) : 0.0)
        push!(hourly_tke_mean, mean(tke[idxs]))
    end

    NRELWindData(times, windspeed, winddirection, tke,
                 hourly_times, hourly_speed_std, hourly_dir_std, hourly_tke_mean)
end

"""Circular standard deviation for wind direction (degrees)."""
function _circular_std(angles::AbstractVector{<:Real})
    s = sum(sind, angles)
    c = sum(cosd, angles)
    n = length(angles)
    R = sqrt(s^2 + c^2) / n  # mean resultant length
    # Circular std dev in degrees (Mardia & Jupp, 2000)
    rad2deg(sqrt(-2 * log(max(R, 1e-10))))
end

"""
    predict(nrel::NRELWindData, dt::DateTime)

Return `(speed_std, dir_std, tke)` at the given time via linear interpolation
of the precomputed hourly statistics.
"""
function GeoSurrogates.predict(nrel::NRELWindData, dt::DateTime)
    times = nrel.hourly_times

    if dt <= first(times)
        return (speed_std=nrel.hourly_speed_std[1],
                dir_std=nrel.hourly_dir_std[1],
                tke=nrel.hourly_tke_mean[1])
    elseif dt >= last(times)
        return (speed_std=nrel.hourly_speed_std[end],
                dir_std=nrel.hourly_dir_std[end],
                tke=nrel.hourly_tke_mean[end])
    end

    i = searchsortedlast(times, dt)
    t1, t2 = times[i], times[i + 1]
    α = Dates.value(dt - t1) / Dates.value(t2 - t1)

    speed_std = (1 - α) * nrel.hourly_speed_std[i] + α * nrel.hourly_speed_std[i + 1]
    dir_std   = (1 - α) * nrel.hourly_dir_std[i]   + α * nrel.hourly_dir_std[i + 1]
    tke       = (1 - α) * nrel.hourly_tke_mean[i]   + α * nrel.hourly_tke_mean[i + 1]

    return (speed_std=speed_std, dir_std=dir_std, tke=tke)
end



#-----------------------------------------------------------------------------# WindField
"""
    WindField(hrrr_data::Raster)

A wrapper around HRRR wind data that provides interpolation in space and time.
Creates RasterWrap objects for each time step's u and v components.

# Usage
    wf = WindField(get_hrrr_data())
    u, v = predict(wf, (lon, lat, datetime))
"""
struct WindField
    times::Vector{DateTime}
    u_wraps::Vector{GeoSurrogates.RasterWrap}
    v_wraps::Vector{GeoSurrogates.RasterWrap}
end

function WindField(hrrr_data::Raster)
    times = collect(lookup(hrrr_data, Ti))
    u_data = hrrr_data[Band=At(:u)]
    v_data = hrrr_data[Band=At(:v)]

    u_wraps = [GeoSurrogates.RasterWrap(u_data[Ti=At(t)]) for t in times]
    v_wraps = [GeoSurrogates.RasterWrap(v_data[Ti=At(t)]) for t in times]

    WindField(times, u_wraps, v_wraps)
end

function GeoSurrogates.predict(wf::WindField, coords::Tuple{<:Real, <:Real, <:DateTime})
    lon, lat, dt = coords
    times = wf.times

    # Find bracketing time indices
    if dt <= first(times)
        u = GeoSurrogates.predict(wf.u_wraps[1], (lon, lat))
        v = GeoSurrogates.predict(wf.v_wraps[1], (lon, lat))
        return (u = u, v = v)
    elseif dt >= last(times)
        u = GeoSurrogates.predict(wf.u_wraps[end], (lon, lat))
        v = GeoSurrogates.predict(wf.v_wraps[end], (lon, lat))
        return (u = u, v = v)
    end

    # Linear interpolation in time
    i = searchsortedlast(times, dt)
    t1, t2 = times[i], times[i + 1]
    α = Dates.value(dt - t1) / Dates.value(t2 - t1)

    u1 = GeoSurrogates.predict(wf.u_wraps[i], (lon, lat))
    u2 = GeoSurrogates.predict(wf.u_wraps[i + 1], (lon, lat))
    v1 = GeoSurrogates.predict(wf.v_wraps[i], (lon, lat))
    v2 = GeoSurrogates.predict(wf.v_wraps[i + 1], (lon, lat))

    u = (1 - α) * u1 + α * u2
    v = (1 - α) * v1 + α * v2

    return (u = u, v = v)
end

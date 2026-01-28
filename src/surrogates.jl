function get_elevation1()
    data = get_landfire_data()
    elev = data.US_ELEV2020
    return GeoSurrogates.RasterWrap(elev)
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

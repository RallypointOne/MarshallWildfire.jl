#-----------------------------------------------------------------------------# FireModel
# Level-set fire propagation model based on Rothermel's equations
# Inspired by: https://github.com/jClugstor/WildfireModels

"""
    FireModel

A fire propagation model using the level-set method with Rothermel fire spread.

# Fields
- `wind_field::WindField` - Spatiotemporal wind interpolation
- `landfire::NamedTuple` - Landfire raster data (fuel, slope, aspect, etc.)
- `start_time::DateTime` - Reference time for the simulation
- `domain::Extent` - Spatial domain
"""
struct FireModel
    wind_field::WindField
    landfire::NamedTuple
    start_time::DateTime
    domain::Extent
end

function FireModel(; wind_data = get_hrrr_data(), landfire_data = get_landfire_data())
    wf = WindField(wind_data)
    FireModel(wf, landfire_data, start_time_utc, extent)
end

#-----------------------------------------------------------------------------# Model Interface Functions
"""
    get_wind(model::FireModel, x, y, t_seconds)

Get wind vector (u, v) at position (x, y) and time t_seconds from start.
"""
function get_wind(model::FireModel, x, y, t_seconds)
    dt = model.start_time + Second(round(Int, t_seconds))
    w = GeoSurrogates.predict(model.wind_field, (x, y, dt))
    return (w.u, w.v)
end

"""
    get_terrain(model::FireModel, x, y)

Get terrain slope and aspect at position (x, y).
Returns (slope_rad, aspect_rad).
"""
function get_terrain(model::FireModel, x, y)
    lf = model.landfire

    if haskey(lf, :SLP) && haskey(lf, :ASP)
        slp = GeoSurrogates.predict(GeoSurrogates.RasterWrap(lf.SLP), (x, y))
        asp = GeoSurrogates.predict(GeoSurrogates.RasterWrap(lf.ASP), (x, y))
        return (deg2rad(slp), deg2rad(asp))
    else
        return (0.0, 0.0)
    end
end

"""
    get_fuel(model::FireModel, x, y)

Get fuel model properties at position (x, y).
Returns a FuelModel from rothermel.jl.
"""
function get_fuel(model::FireModel, x, y)
    lf = model.landfire

    # Look for LANDFIRE fuel model layer (may have US_250 prefix)
    fuel_key = nothing
    for k in keys(lf)
        ks = string(k)
        if occursin("FBFM13", ks) || occursin("FBFM40", ks)
            fuel_key = k
            break
        end
    end

    if !isnothing(fuel_key)
        fuel_code = round(Int, GeoSurrogates.predict(GeoSurrogates.RasterWrap(lf[fuel_key]), (x, y)))
        fuel = get_fuel_model(fuel_code)
        return isnothing(fuel) ? DEFAULT_FUEL : fuel
    else
        return DEFAULT_FUEL
    end
end

"""
    spread_rate(model::FireModel, x, y, t; Mf=0.05)

Calculate fire spread rate and direction at (x, y, t).
"""
function spread_rate(model::FireModel, x, y, t; Mf=0.05)
    u, v = get_wind(model, x, y, t)

    # Guard against NaN wind values
    if !isfinite(u) || !isfinite(v)
        return (rate = 0.0, direction = 0.0, R0 = 0.0)
    end

    wind_speed = sqrt(u^2 + v^2)
    wind_dir = atan(v, u)

    slope, aspect = get_terrain(model, x, y)

    # Guard against NaN terrain values
    if !isfinite(slope) || !isfinite(aspect)
        slope, aspect = 0.0, 0.0
    end

    fuel = get_fuel(model, x, y)

    result = rothermel_spread_rate(fuel, wind_speed, wind_dir, slope, aspect; Mf)

    # Guard against NaN results
    if !isfinite(result.rate)
        return (rate = 0.0, direction = result.direction, R0 = result.R0)
    end

    return result
end

#-----------------------------------------------------------------------------# Fire Simulation using Level Set

"""
    init_fire!(ls::LevelSet, model::FireModel)

Initialize a level set for fire simulation from the ignition point.
"""
function init_fire!(ls::LevelSet)
    init_signed_distance!(ls, ignition_point.lon, ignition_point.lat)
end

"""
    create_speed_function(model::FireModel; Mf=0.05)

Create a speed function for use with LevelSet simulation.
Returns a function `speed(x, y, t) -> Real`.
"""
function create_speed_function(model::FireModel; Mf::Real=0.05)
    return (x, y, t) -> spread_rate(model, x, y, t; Mf).rate
end

"""
    simulate(model::FireModel, duration_seconds; dt=60.0, nx=100, ny=100, Mf=0.05)

Run fire simulation for given duration using the level set method.

# Arguments
- `model`: FireModel containing wind, terrain, and fuel data
- `duration_seconds`: Total simulation time in seconds

# Keyword Arguments
- `dt`: Time step in seconds (default: 60.0)
- `nx`, `ny`: Grid resolution (default: 100×100)
- `Mf`: Fuel moisture content (default: 0.05)

# Returns
Vector of `(time, LevelSet)` snapshots.
"""
function simulate(model::FireModel, duration_seconds::Real;
                  dt::Real=60.0, nx::Int=100, ny::Int=100, Mf::Real=0.05)
    # Create level set on model domain
    ls = LevelSet(model.domain; nx, ny)

    # Initialize from ignition point
    init_fire!(ls)

    # Create speed function from model
    speed_fn = create_speed_function(model; Mf)

    # Estimate max speed for CFL (use a reasonable upper bound)
    max_speed = 10.0  # m/s, conservative estimate

    # Run simulation
    simulate(ls, speed_fn, duration_seconds;
             dt=dt, save_interval=dt, max_speed=max_speed)
end

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

#-----------------------------------------------------------------------------# Fuel Models
# Anderson 13 fuel model properties
# (w0: fuel load kg/m², σ: SAV ratio 1/m, Mx: moisture of extinction, h: heat content kJ/kg, δ: fuel depth m)
const FUEL_MODELS = Dict(
    1  => (w0=0.166, σ=11483.0, Mx=0.12, h=18622.0, δ=0.305),   # Short grass
    2  => (w0=0.896, σ=9843.0,  Mx=0.15, h=18622.0, δ=0.305),   # Timber grass
    3  => (w0=1.345, σ=4921.0,  Mx=0.25, h=18622.0, δ=0.762),   # Tall grass
    4  => (w0=2.242, σ=6562.0,  Mx=0.20, h=18622.0, δ=1.829),   # Chaparral
    5  => (w0=0.448, σ=6562.0,  Mx=0.20, h=18622.0, δ=0.610),   # Brush
    6  => (w0=0.673, σ=5741.0,  Mx=0.25, h=18622.0, δ=0.762),   # Dormant brush
    7  => (w0=0.507, σ=5741.0,  Mx=0.40, h=18622.0, δ=0.762),   # Southern rough
    8  => (w0=0.673, σ=6562.0,  Mx=0.30, h=18622.0, δ=0.061),   # Compact timber litter
    9  => (w0=1.307, σ=8203.0,  Mx=0.25, h=18622.0, δ=0.061),   # Hardwood litter
    10 => (w0=1.345, σ=6562.0,  Mx=0.25, h=18622.0, δ=0.305),   # Timber understory
    11 => (w0=0.673, σ=4921.0,  Mx=0.15, h=18622.0, δ=0.305),   # Light logging slash
    12 => (w0=1.794, σ=4921.0,  Mx=0.20, h=18622.0, δ=0.701),   # Medium logging slash
    13 => (w0=3.140, σ=4921.0,  Mx=0.25, h=18622.0, δ=0.914),   # Heavy logging slash
)

const DEFAULT_FUEL = FUEL_MODELS[1]
const ρp = 513.0  # Ovendry particle density (kg/m³)

#-----------------------------------------------------------------------------# Rothermel Fire Spread Rate
"""
    rothermel_spread_rate(fuel, wind_speed, wind_dir, slope, aspect; Mf=0.05)

Calculate fire spread rate using Rothermel's model.

# Arguments
- `fuel`: Fuel model properties (w0, σ, Mx, h, δ)
- `wind_speed`: Wind speed (m/s)
- `wind_dir`: Wind direction (radians, direction wind is coming FROM)
- `slope`: Terrain slope (radians)
- `aspect`: Terrain aspect (radians, downslope direction)
- `Mf`: Fuel moisture content (fraction, default 0.05)

# Returns
- Spread rate (m/s) and direction (radians)
"""
function rothermel_spread_rate(fuel, wind_speed, wind_dir, slope, aspect; Mf=0.05)
    (; w0, σ, Mx, h, δ) = fuel

    # Packing ratio
    β = w0 / (ρp * δ)
    βop = 3.348 * σ^(-0.8189)  # Optimum packing ratio
    βratio = β / βop

    # Reaction intensity
    Γmax = σ^1.5 / (495.0 + 0.0594 * σ^1.5)
    A = 133.0 * σ^(-0.7913)
    Γ = Γmax * (βratio)^A * exp(A * (1.0 - βratio))

    # Moisture damping
    ηM = clamp(1.0 - 2.59 * (Mf / Mx) + 5.11 * (Mf / Mx)^2 - 3.52 * (Mf / Mx)^3, 0.0, 1.0)

    # Mineral damping (constant for most fuels)
    ηs = 0.174 * 0.01^(-0.19)  # Se = 0.01 typical

    # Reaction intensity (kW/m²)
    IR = Γ * w0 * h * ηM * ηs

    # Propagating flux ratio
    ξ = exp((0.792 + 0.681 * σ^0.5) * (β + 0.1)) / (192.0 + 0.2595 * σ)

    # Heat sink (effective heating number × heat of preignition)
    ε = exp(-138.0 / σ)  # Effective heating number
    Qig = 250.0 + 1116.0 * Mf  # Heat of preignition (kJ/kg)
    ρb = w0 / δ  # Bulk density

    # No-wind, no-slope spread rate (m/min -> m/s)
    R0 = (IR * ξ) / (ρb * ε * Qig) / 60.0

    # Wind factor
    C = 7.47 * exp(-0.133 * σ^0.55)
    B = 0.02526 * σ^0.54
    E = 0.715 * exp(-3.59e-4 * σ)
    U = wind_speed * 60.0  # Convert to m/min for formula
    ϕw = C * U^B * βratio^(-E)

    # Slope factor
    ϕs = 5.275 * β^(-0.3) * tan(slope)^2

    # Combined spread rate
    # Wind pushes fire in direction it's blowing TO (opposite of wind_dir)
    # Slope pushes fire uphill (opposite of aspect which points downhill)
    wind_push_dir = wind_dir + π
    slope_push_dir = aspect + π

    # Vector combination of wind and slope effects
    Rw = R0 * ϕw
    Rs = R0 * ϕs

    Rx = R0 + Rw * cos(wind_push_dir) + Rs * cos(slope_push_dir)
    Ry = Rw * sin(wind_push_dir) + Rs * sin(slope_push_dir)

    R = sqrt(Rx^2 + Ry^2)
    θ = atan(Ry, Rx)

    return (rate = R, direction = θ)
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
"""
function get_fuel(model::FireModel, x, y)
    lf = model.landfire

    fuel_key = haskey(lf, :FBFM13) ? :FBFM13 : (haskey(lf, :FBFM40) ? :FBFM40 : nothing)

    if !isnothing(fuel_key)
        fuel_code = round(Int, GeoSurrogates.predict(GeoSurrogates.RasterWrap(lf[fuel_key]), (x, y)))
        return get(FUEL_MODELS, fuel_code, DEFAULT_FUEL)
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
    wind_speed = sqrt(u^2 + v^2)
    wind_dir = atan(v, u)

    slope, aspect = get_terrain(model, x, y)
    fuel = get_fuel(model, x, y)

    rothermel_spread_rate(fuel, wind_speed, wind_dir, slope, aspect; Mf)
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

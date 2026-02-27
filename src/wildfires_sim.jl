#=
    Marshall Fire simulation using the Wildfires.jl propagation engine.

    Bridges MarshallWildfire data (HRRR wind, LANDFIRE terrain/fuel) to
    Wildfires.jl's level set method and Rothermel spread rate model.
=#

# Coordinate conversion constants at ~40°N
const M_PER_DEG_LAT = 111_320.0
const M_PER_DEG_LON = 111_320.0 * cosd(40.0)  # ≈ 85,280 m

# Map Anderson 13 fuel codes to Wildfires.Rothermel fuel models
const WF_FUEL_MAP = Dict{Int, Wildfires.Rothermel.Rothermel}(
    1  => Wildfires.Rothermel.SHORT_GRASS,
    2  => Wildfires.Rothermel.TIMBER_GRASS,
    3  => Wildfires.Rothermel.TALL_GRASS,
    4  => Wildfires.Rothermel.CHAPARRAL,
    5  => Wildfires.Rothermel.BRUSH,
    6  => Wildfires.Rothermel.DORMANT_BRUSH,
    7  => Wildfires.Rothermel.SOUTHERN_ROUGH,
    8  => Wildfires.Rothermel.CLOSED_TIMBER_LITTER,
    9  => Wildfires.Rothermel.HARDWOOD_LITTER,
    10 => Wildfires.Rothermel.TIMBER_UNDERSTORY,
    11 => Wildfires.Rothermel.LIGHT_SLASH,
    12 => Wildfires.Rothermel.MEDIUM_SLASH,
    13 => Wildfires.Rothermel.HEAVY_SLASH,
)

# Key Marshall Fire events for timeline annotation.
# Times in minutes from simulation start (10:00 AM MST, Dec 30, 2021).
# Sources: NWS Boulder, Colorado Sun, Boulder County After-Action Report.
const MARSHALL_FIRE_EVENTS = [
    (t = 81.0,   label = "Fire discovered"),        # 11:21 AM — Engine 2209 finds grass fire
    (t = 135.0,  label = "Sagamore burns"),          # 12:15 PM — first neighborhood destroyed
    (t = 210.0,  label = "Louisville reached"),      # 1:30 PM — fire enters Louisville
    (t = 315.0,  label = "Emergency declared"),      # 3:15 PM — Governor declares state of emergency
    (t = 420.0,  label = "Winds subside"),           # 5:00 PM — forward spread largely stops
    (t = 840.0,  label = "Cold front"),              # 12:00 AM Dec 31 — wind shift ends fire threat
    (t = 1680.0, label = "Snow begins"),             # 2:00 PM Dec 31 — snow blankets fire area
]

#-----------------------------------------------------------------------------# MarshallSpreadModel
"""
    MarshallSpreadModel

Spread rate model that bridges MarshallWildfire data to Wildfires.jl's Rothermel
implementation.  Callable as `model(t, x, y) → spread_rate [m/min]` where `(x, y)`
are meters relative to the ignition point and `t` is simulation time in minutes.

Used with `Wildfires.SpreadModel.simulate!` to drive the level set propagation.
"""
struct MarshallSpreadModel{T<:AbstractFloat}
    wind_field::WindField
    slp_wrap::GeoSurrogates.RasterWrap
    asp_wrap::GeoSurrogates.RasterWrap
    fuel_wrap::GeoSurrogates.RasterWrap
    start_time::DateTime
    lon0::T
    lat0::T
    moisture::Wildfires.Rothermel.FuelClasses{T}
    lon_range::Tuple{T, T}
    lat_range::Tuple{T, T}
    wind_adjustment::T  # 10m-to-midflame wind speed multiplier
    nrel_wind::Union{Nothing, NRELWindData}  # optional NREL WTK-LED wind variability data
    lb_formula::Symbol  # length-to-breadth formula for elliptical spread (:anderson or :green)
end

"""
    MarshallSpreadModel(; wind_data, landfire_data, Mf=0.03, wind_adjustment=0.36, nrel_wind=nothing)

# Keyword Arguments
- `wind_adjustment`: Multiplier converting 10m wind speed to midflame wind speed.
  HRRR provides wind at 10m above ground, but the Rothermel model requires the
  midflame wind speed (at ~1.5× fuel bed depth).  Standard adjustment factors
  from Anderson (1982):
    - 0.36 for open grassland (NFFL 1, unsheltered)
    - 0.28 for brush/shrub (NFFL 4–7, partially sheltered)
    - 0.10 for closed timber (NFFL 8–10, fully sheltered)
  Default 0.36 is appropriate since the Marshall Fire area is predominantly grass.
- `nrel_wind`: Optional `NRELWindData` providing sub-hourly wind variability from
  NREL WTK-LED. When present, wind speed std dev (σ_u ≈ √(2/3 · TKE)) is available
  for uncertainty quantification and ensemble runs.
"""
function MarshallSpreadModel(;
    wind_data = get_hrrr_data(),
    landfire_data = get_landfire_data(),
    Mf::Real = 0.03,
    wind_adjustment::Real = 0.36,
    nrel_wind::Union{Nothing, NRELWindData} = nothing,
    T::Type{<:AbstractFloat} = Float64,
    lb_formula::Symbol = :anderson,
)
    wf = WindField(wind_data)
    moisture = Wildfires.Rothermel.FuelClasses(
        d1=T(Mf), d10=T(Mf)+T(0.01), d100=T(Mf)+T(0.02),
        herb=zero(T), wood=zero(T),
    )

    # Find LANDFIRE layers by substring matching on key names
    function _find_layer(data, patterns...)
        for k in keys(data), p in patterns
            occursin(p, string(k)) && return data[k]
        end
        error("No LANDFIRE layer matching $(patterns) found. Available: $(keys(data))")
    end

    slp_raster = _find_layer(landfire_data, "SLP")
    slp_wrap = GeoSurrogates.RasterWrap(slp_raster)
    asp_wrap = GeoSurrogates.RasterWrap(_find_layer(landfire_data, "ASP"))
    fuel_wrap = GeoSurrogates.RasterWrap(_find_layer(landfire_data, "FBFM13", "FBFM40"))

    # Extract data bounds from raster for out-of-bounds checking
    xs = lookup(slp_raster, X)
    ys = lookup(slp_raster, Y)
    lon_range = (T(minimum(xs)), T(maximum(xs)))
    lat_range = (T(minimum(ys)), T(maximum(ys)))

    MarshallSpreadModel(
        wf, slp_wrap, asp_wrap, fuel_wrap,
        start_time_utc,
        T(ignition_point.lon), T(ignition_point.lat),
        moisture, lon_range, lat_range, T(wind_adjustment),
        nrel_wind, lb_formula,
    )
end

# No dynamic state to update between time steps
Wildfires.SpreadModel.update!(::MarshallSpreadModel, ::Wildfires.LevelSet.LevelSetGrid, dt) = nothing

"""
    wind_variability(model::MarshallSpreadModel, t_min)

Return wind variability metrics at simulation time `t_min` (minutes from start).
Uses NREL WTK-LED data if available.

Returns a named tuple `(σ_u, dir_std, tke)` where:
- `σ_u`: Wind speed std dev (m/s), derived from TKE via `σ_u ≈ √(2/3 · TKE)`
- `dir_std`: Wind direction std dev (degrees)
- `tke`: Mean turbulent kinetic energy (m²/s²)

Returns `nothing` if no NREL data is attached.
"""
function wind_variability(model::MarshallSpreadModel, t_min::Real)
    model.nrel_wind === nothing && return nothing
    dt = model.start_time + Second(round(Int, t_min * 60))
    v = GeoSurrogates.predict(model.nrel_wind, dt)
    σ_u = sqrt(2/3 * v.tke)  # TKE → wind speed std dev
    return (σ_u=σ_u, dir_std=v.dir_std, tke=v.tke)
end

#-----------------------------------------------------------------------------# Directional spread_rate_field!
# Override the generic spread_rate_field! to produce directional (anisotropic) spread.
#
# The generic fallback applies the head-fire rate (maximum, wind-aligned) isotropically
# to all cells, which dramatically overestimates flanking and backing fire spread.
# This override mirrors the approach in Wildfires.SpreadModel._directional_rate:
#
#   F = R_base + (R_head - R_base) · max(0, n̂ · d̂)
#
# where n̂ = ∇φ/|∇φ| is the front propagation direction and d̂ is the combined
# wind + slope push direction.  Cells propagating into the wind get only R_base
# (the no-wind, no-slope rate), while cells propagating downwind get R_head.
function Wildfires.SpreadModel.spread_rate_field!(
    F::AbstractMatrix{FT}, model::MarshallSpreadModel, grid::Wildfires.LevelSet.LevelSetGrid,
) where {FT}
    xs = Wildfires.LevelSet.xcoords(grid)
    ys = Wildfires.LevelSet.ycoords(grid)
    t  = grid.t
    φ  = grid.φ
    nrows, ncols = size(φ)
    dxg, dyg = grid.dx, grid.dy

    for j in eachindex(xs), i in eachindex(ys)
        x, y = xs[j], ys[i]
        lon, lat = _to_lonlat(model, x, y)

        # Out of LANDFIRE bounds → no spread
        if !(model.lon_range[1] <= lon <= model.lon_range[2] &&
             model.lat_range[1] <= lat <= model.lat_range[2])
            F[i, j] = zero(FT)
            continue
        end

        # Fuel lookup
        code = round(Int, GeoSurrogates.predict(model.fuel_wrap, (lon, lat)))
        if code in NON_BURNABLE_CODES
            F[i, j] = zero(FT)
            continue
        end
        fuel = get(WF_FUEL_MAP, code, Wildfires.Rothermel.SHORT_GRASS)
        moist = model.moisture

        # Wind (10m → midflame adjusted)
        dt_time = model.start_time + Second(round(Int, t * 60))
        w = GeoSurrogates.predict(model.wind_field, (lon, lat, dt_time))
        u, v = w.u, w.v
        if !isfinite(u) || !isfinite(v)
            F[i, j] = zero(FT)
            continue
        end
        wind_kmh = sqrt(u^2 + v^2) * 3.6 * model.wind_adjustment

        # Wind direction: HRRR (u,v) points where wind blows TO.
        # Rothermel convention: direction wind blows FROM.
        wind_from = atan(v, u) + π

        # Terrain
        slp_deg = GeoSurrogates.predict(model.slp_wrap, (lon, lat))
        slope = isfinite(slp_deg) ? tan(deg2rad(clamp(slp_deg, 0.0, 80.0))) : 0.0
        asp_deg = GeoSurrogates.predict(model.asp_wrap, (lon, lat))
        aspect = isfinite(asp_deg) ? deg2rad(asp_deg) : 0.0

        R_base = Wildfires.Rothermel.rate_of_spread(fuel; moisture=moist, wind=0.0, slope=0.0)

        # Front normal from ∇φ (central differences)
        dφdx = if j == 1
            (φ[i, 2] - φ[i, 1]) / dxg
        elseif j == ncols
            (φ[i, ncols] - φ[i, ncols-1]) / dxg
        else
            (φ[i, j+1] - φ[i, j-1]) / (2dxg)
        end
        dφdy = if i == 1
            (φ[2, j] - φ[1, j]) / dyg
        elseif i == nrows
            (φ[nrows, j] - φ[nrows-1, j]) / dyg
        else
            (φ[i+1, j] - φ[i-1, j]) / (2dyg)
        end
        grad = hypot(dφdx, dφdy)

        if grad == 0
            R_head = Wildfires.Rothermel.rate_of_spread(fuel; moisture=moist, wind=wind_kmh, slope=slope)
            F[i, j] = iszero(R_head) ? zero(FT) : R_head
            continue
        end

        nx, ny = dφdx / grad, dφdy / grad

        # Directional slope: project terrain slope onto front propagation direction.
        # Only apply slope boost when fire propagates uphill; downhill gets slope=0.
        # Aspect = downhill direction, so uphill = (-cos(aspect), -sin(aspect)).
        uphill_dot = nx * (-cos(aspect)) + ny * (-sin(aspect))
        slope_dir = slope * max(0.0, uphill_dot)

        # Head-fire rate uses directional slope (not raw terrain steepness)
        R_head = Wildfires.Rothermel.rate_of_spread(fuel; moisture=moist, wind=wind_kmh, slope=slope_dir)

        if iszero(R_head) || R_head ≈ R_base
            F[i, j] = R_head
            continue
        end

        # Push direction: weighted combination of wind push and slope push.
        # Wind pushes opposite to FROM direction; slope pushes uphill (opposite aspect).
        # Push direction uses full slope to define where fire "wants" to go.
        R_w = Wildfires.Rothermel.rate_of_spread(fuel; moisture=moist, wind=wind_kmh, slope=0.0)
        R_s = Wildfires.Rothermel.rate_of_spread(fuel; moisture=moist, wind=0.0, slope=slope)
        w_wind  = R_w - R_base
        w_slope = R_s - R_base

        px = w_wind * (-cos(wind_from)) + w_slope * (-cos(aspect))
        py = w_wind * (-sin(wind_from)) + w_slope * (-sin(aspect))
        pmag = hypot(px, py)

        if pmag == 0
            F[i, j] = R_head
            continue
        end

        cos_theta = (nx * px + ny * py) / pmag

        # Elliptical spread model (Anderson 1983)
        U_ms = wind_kmh / 3.6  # midflame km/h → m/s
        LB = Wildfires.SpreadModel.length_to_breadth(U_ms; formula=model.lb_formula)
        ε = Wildfires.SpreadModel.fire_eccentricity(LB)
        F[i, j] = R_head * (1 - ε) / (1 - ε * cos_theta)
    end
    F
end

# Convert grid coordinates (meters from ignition) to lon/lat
@inline function _to_lonlat(m::MarshallSpreadModel, x, y)
    (m.lon0 + x / M_PER_DEG_LON, m.lat0 + y / M_PER_DEG_LAT)
end

function (m::MarshallSpreadModel{T})(t, x, y) where {T}
    lon, lat = _to_lonlat(m, x, y)

    # Return 0 for cells outside the LANDFIRE data extent
    (m.lon_range[1] <= lon <= m.lon_range[2] && m.lat_range[1] <= lat <= m.lat_range[2]) || return zero(T)

    # Fuel lookup
    code = round(Int, GeoSurrogates.predict(m.fuel_wrap, (lon, lat)))
    code in NON_BURNABLE_CODES && return zero(T)
    fuel = get(WF_FUEL_MAP, code, Wildfires.Rothermel.SHORT_GRASS)

    # Wind: HRRR returns 10m (u, v) in m/s; Wildfires.Rothermel wants midflame km/h.
    # Apply wind_adjustment to convert 10m → midflame height (see Anderson 1982,
    # Table 3).  Without this correction, the Rothermel wind factor—which scales as
    # U^B—dramatically overestimates spread rate (e.g. 1400 vs 170 m/min for grass).
    dt = m.start_time + Second(round(Int, t * 60))
    w = GeoSurrogates.predict(m.wind_field, (lon, lat, dt))
    u, v = w.u, w.v
    (!isfinite(u) || !isfinite(v)) && return zero(T)
    wind_kmh = sqrt(u^2 + v^2) * 3.6 * m.wind_adjustment

    # Terrain: LANDFIRE slope in degrees → rise/run fraction
    slp_deg = GeoSurrogates.predict(m.slp_wrap, (lon, lat))
    slope = isfinite(slp_deg) ? tan(deg2rad(clamp(slp_deg, 0.0, 80.0))) : 0.0

    # Spread rate [m/min]
    T(Wildfires.Rothermel.rate_of_spread(fuel; moisture=m.moisture, wind=wind_kmh, slope=slope))
end

#-----------------------------------------------------------------------------# simulate_marshall_wf
"""
    simulate_marshall_wf(; duration_min=120, nx=200, ny=200, dx=50.0, Mf=0.03,
                           ignition_radius=200.0, cfl=0.5, reinit_every=10,
                           save_every_min=5.0,
                           output_dir=joinpath(@__DIR__, "..", "report", "images"))

Simulate the Marshall Fire using the Wildfires.jl level set propagation engine with
real HRRR wind, LANDFIRE terrain, and LANDFIRE fuel data.

Saves snapshots every `save_every_min` minutes, then produces:
- A static PNG of the final perimeter overlaid on the observed perimeter.
- An animated GIF showing fire growth over time.

# Returns
A named tuple `(; grid, model, snapshots)`.
"""
function simulate_marshall_wf(;
    duration_min::Real = 120.0,
    nx::Int = 300,
    ny::Int = 200,
    dx::Real = 50.0,
    Mf::Real = 0.03,
    wind_adjustment::Real = 0.36,
    nrel_wind::Union{Nothing, NRELWindData} = nothing,
    ignition_radius::Real = 200.0,
    cfl::Real = 0.5,
    reinit_every::Int = 10,
    save_every_min::Real = 1.0,
    output_dir::String = joinpath(@__DIR__, "..", "report", "images"),
    prefix::String = "wildfires_sim",
    T::Type{<:AbstractFloat} = Float64,
    lb_formula::Symbol = :anderson,
)
    model = MarshallSpreadModel(; Mf, wind_adjustment, nrel_wind, T, lb_formula)

    # Grid offset so ignition (0, 0) is in the western third — fire spreads east
    dx_T = T(dx)
    x0 = -nx * dx_T / 3
    y0 = -ny * dx_T / 2
    grid = Wildfires.LevelSet.LevelSetGrid(nx, ny; dx=dx_T, x0, y0)

    # Ignite at center (radius must exceed dx * √2/2 to capture nearest cell)
    r = max(T(ignition_radius), dx_T * T(1.5))
    Wildfires.LevelSet.ignite!(grid, zero(T), zero(T), r)

    # ----- Precompute burnable mask -----
    # Non-burnable cells (water, urban, etc.) act as firebreaks: φ is clamped ≥ 0
    # so the fire front can never enter them.  Without this, reinitialization can
    # propagate negative φ into water/lake cells even though F = 0 there.
    xs = Wildfires.LevelSet.xcoords(grid)
    ys = Wildfires.LevelSet.ycoords(grid)
    burnable = trues(size(grid))
    for j in eachindex(xs), i in eachindex(ys)
        lon, lat = _to_lonlat(model, xs[j], ys[i])
        if model.lon_range[1] <= lon <= model.lon_range[2] &&
           model.lat_range[1] <= lat <= model.lat_range[2]
            code = round(Int, GeoSurrogates.predict(model.fuel_wrap, (lon, lat)))
            if code in NON_BURNABLE_CODES
                burnable[i, j] = false
            end
        else
            burnable[i, j] = false
        end
    end
    n_unburnable = count(!, burnable)
    n_unburnable > 0 && println("  $n_unburnable / $(length(burnable)) cells marked non-burnable (water, urban, etc.)")

    # ----- Run simulation with snapshot saving -----
    println("Simulating Marshall Fire: $(nx)×$(ny) grid, dx=$(dx)m, $(duration_min) min...")

    F = similar(grid.φ)
    snapshots = [(grid.t, copy(grid.φ))]  # (time_min, φ_copy)
    next_save = save_every_min
    step = 0
    t_start = time()
    last_print = t_start

    while grid.t < duration_min
        Wildfires.SpreadModel.spread_rate_field!(F, model, grid)
        step_dt = Wildfires.LevelSet.cfl_dt(grid, F; cfl=T(cfl))
        isinf(step_dt) && break
        Wildfires.SpreadModel.update!(model, grid, step_dt)
        Wildfires.LevelSet.advance!(grid, F, step_dt)
        step += 1
        if step % reinit_every == 0
            # Reinitialization restores |∇φ| ≈ 1 but its pseudo-timestep can flip
            # the sign of cells barely inside the burned region (φ slightly < 0).
            # Preserve the sign so that once a cell burns it stays burned.
            sign_before = sign.(grid.φ)
            Wildfires.LevelSet.reinitialize!(grid)
            @. grid.φ = ifelse(sign_before < 0 && grid.φ > 0, -eps(T), grid.φ)
        end

        # Enforce non-burnable barrier: clamp φ ≥ 0 for water, urban, etc.
        @. grid.φ = ifelse(!burnable && grid.φ < 0, zero(T), grid.φ)

        if grid.t >= next_save
            push!(snapshots, (grid.t, copy(grid.φ)))
            next_save += save_every_min
        end

        # Progress meter (update every 10 seconds)
        now = time()
        if now - last_print >= 10.0
            pct = grid.t / duration_min * 100
            elapsed = now - t_start
            eta = pct > 0 ? elapsed / pct * (100 - pct) : NaN
            area = Wildfires.LevelSet.burn_area(grid) / 1e6
            print("\r  $(round(pct, digits=1))% | t=$(round(grid.t, digits=1))/$(duration_min) min | " *
                  "burned=$(round(area, digits=2)) km² | " *
                  "elapsed=$(round(elapsed, digits=0))s | ETA=$(round(eta, digits=0))s   ")
            last_print = now
        end
    end
    println()  # newline after progress meter

    # Always save final state
    if isempty(snapshots) || last(snapshots)[1] < grid.t
        push!(snapshots, (grid.t, copy(grid.φ)))
    end

    area_km2 = Wildfires.LevelSet.burn_area(grid) / 1e6
    elapsed_total = round(time() - t_start, digits=1)
    println("Done! $(step) steps, t=$(round(grid.t, digits=1)) min, burned ≈ $(round(area_km2, digits=2)) km² ($(elapsed_total)s wall time)")

    # ----- Produce visualizations -----
    mkpath(output_dir)
    _plot_final_perimeter(grid, model, snapshots, output_dir; prefix)
    _plot_fire_animation(grid, model, snapshots, output_dir; prefix)

    return (; grid, model, snapshots)
end

#-----------------------------------------------------------------------------# simulate_marshall_wf_full
"""
    simulate_marshall_wf_full(; kwargs...)

Simulate the full Marshall Fire period (~44 hours, Dec 30 2021 10:00 AM – Jan 1 2022 6:00 AM Denver).

Convenience wrapper around `simulate_marshall_wf` with:
- `duration_min = 2640` (44 hours)
- `save_every_min = 15.0` (176 frames, ~15s at 12fps)
- Output files prefixed `wildfires_sim_full_`

All other keyword arguments are forwarded to `simulate_marshall_wf`.
"""
function simulate_marshall_wf_full(; kwargs...)
    defaults = (duration_min = 2640.0, save_every_min = 15.0, prefix = "wildfires_sim_full")
    merged = merge(defaults, NamedTuple(kwargs))
    simulate_marshall_wf(; merged...)
end

#-----------------------------------------------------------------------------# Visualization helpers

# Convert grid x/y arrays to lon/lat for plotting
function _grid_lonlats(grid, model::MarshallSpreadModel)
    xs = collect(Wildfires.LevelSet.xcoords(grid))
    ys = collect(Wildfires.LevelSet.ycoords(grid))
    lons = model.lon0 .+ xs ./ M_PER_DEG_LON
    lats = model.lat0 .+ ys ./ M_PER_DEG_LAT
    return lons, lats
end

function _plot_final_perimeter(grid, model, snapshots, output_dir; prefix="wildfires_sim")
    println("  Saving final perimeter image...")
    lons, lats = _grid_lonlats(grid, model)
    _, φ_final = last(snapshots)
    t_final = last(snapshots)[1]

    perim = get_perimeter()

    # Axis limits: encompass observed perimeter with padding
    perim_ext = GI.extent(perim)
    pad = 0.005
    xlo = min(perim_ext.X[1], minimum(lons)) - pad
    xhi = max(perim_ext.X[2], maximum(lons)) + pad
    ylo = min(perim_ext.Y[1], minimum(lats)) - pad
    yhi = max(perim_ext.Y[2], maximum(lats)) + pad

    fig = Figure(size = (900, 700))
    ax = Axis(fig[1, 1];
        title = "Marshall Fire — Wildfires.jl Simulation (t = $(round(t_final, digits=1)) min)",
        xlabel = "Longitude", ylabel = "Latitude", aspect = DataAspect(),
    )
    xlims!(ax, xlo, xhi)
    ylims!(ax, ylo, yhi)

    # Burned/unburned heatmap
    heatmap!(ax, lons, lats, φ_final'; colormap = :RdYlBu, colorrange = (-500, 500))

    # Simulated fire front
    contour!(ax, lons, lats, φ_final'; levels = [0.0], color = :black, linewidth = 2.5)

    # Observed perimeter
    poly!(ax, perim.geometry; color = (:transparent), strokecolor = :red, strokewidth = 2, linestyle = :dash)

    # Ignition point
    scatter!(ax, [model.lon0], [model.lat0]; color = :orange, markersize = 15, marker = :star5)

    # Legend
    elem_sim = LineElement(color = :black, linewidth = 2.5)
    elem_obs = LineElement(color = :red, linewidth = 2, linestyle = :dash)
    Legend(fig[1, 2], [elem_sim, elem_obs], ["Simulated front", "Observed perimeter"])

    path = joinpath(output_dir, "$(prefix)_final.png")
    save(path, fig)
    println("  Saved: $path")
end

function _plot_fire_animation(grid, model, snapshots, output_dir; prefix="wildfires_sim")
    println("  Saving fire growth animation...")
    lons, lats = _grid_lonlats(grid, model)
    perim = get_perimeter()
    buildings = get_building_footprints()

    # Axis limits: focus on observed perimeter extent with padding
    perim_ext = GI.extent(perim)
    pad = 0.005  # ~0.5 km padding
    xlo = perim_ext.X[1] - pad
    xhi = perim_ext.X[2] + pad
    ylo = perim_ext.Y[1] - pad
    yhi = perim_ext.Y[2] + pad

    # Precompute cumulative burned masks
    cumulative_burned = [falses(size(snapshots[1][2])) for _ in snapshots]
    for (idx, (_, φ)) in enumerate(snapshots)
        if idx == 1
            cumulative_burned[idx] .= φ .< 0
        else
            cumulative_burned[idx] .= cumulative_burned[idx-1] .| (φ .< 0)
        end
    end

    # ----- Precompute wind field on a grid matching the marshall extent -----
    wf = model.wind_field
    wind_nx, wind_ny = 20, 15
    wind_lons = range(xlo, xhi, length=wind_nx)
    wind_lats = range(ylo, yhi, length=wind_ny)

    # Precompute wind for each snapshot time
    wind_frames = map(snapshots) do (t_min, _)
        dt = model.start_time + Second(round(Int, t_min * 60))
        u = [GeoSurrogates.predict(wf, (lo, la, dt)).u for lo in wind_lons, la in wind_lats]
        v = [GeoSurrogates.predict(wf, (lo, la, dt)).v for lo in wind_lons, la in wind_lats]
        mag = sqrt.(u .^ 2 .+ v .^ 2)
        (; u, v, mag)
    end
    mag_max = maximum(maximum(f.mag) for f in wind_frames)

    # ----- Load elevation raster for terrain subplot -----
    landfire_data = get_landfire_data()
    function _find_layer(data, patterns...)
        for k in keys(data), p in patterns
            occursin(p, string(k)) && return data[k]
        end
        error("No LANDFIRE layer found for $patterns")
    end
    elev_raster = _find_layer(landfire_data, "ELEV")
    elev_xs = collect(lookup(elev_raster, X))
    elev_ys = collect(lookup(elev_raster, Y))
    elev_data = Float64.(elev_raster.data)

    # ----- Build figure with 3 subplots + timeline -----
    fig = Figure(size = (1800, 800))

    ax1 = Axis(fig[1, 1]; title = "Fire Growth", xlabel = "Lon", ylabel = "Lat", aspect = DataAspect())
    ax2 = Axis(fig[1, 2]; title = "Wind Speed & Direction", xlabel = "Lon", ylabel = "Lat", aspect = DataAspect())
    ax3 = Axis(fig[1, 3]; title = "Terrain & Buildings", xlabel = "Lon", ylabel = "Lat", aspect = DataAspect())

    for ax in (ax1, ax2, ax3)
        xlims!(ax, xlo, xhi)
        ylims!(ax, ylo, yhi)
    end

    # Static terrain + buildings subplot (base layer; fire contour added per frame)
    hm_elev = heatmap!(ax3, elev_xs, elev_ys, elev_data'; colormap = :terrain, colorrange = (Float64(minimum(elev_data)), Float64(maximum(elev_data))))
    Colorbar(fig[2, 3], hm_elev; label = "Elevation (m)", vertical = false, flipaxis = false)

    # Placeholder colorbars for wind (created from first frame, updated implicitly)
    hm_wind_ref = heatmap!(ax2, collect(wind_lons), collect(wind_lats), wind_frames[1].mag';
                           colormap = :viridis, colorrange = (0, mag_max))
    Colorbar(fig[2, 2], hm_wind_ref; label = "Wind Speed (m/s)", vertical = false, flipaxis = false)

    # ----- Timeline axis (row 3, spanning all columns) -----
    total_duration = last(snapshots)[1]
    visible_events = [(; t=e.t, label=e.label) for e in MARSHALL_FIRE_EVENTS if e.t <= total_duration]

    ax_tl = Axis(fig[3, 1:3];
        xlabel = "Denver Time (MST)",
        xticklabelsize = 10, xticklabelrotation = π/6,
    )
    hideydecorations!(ax_tl)
    hidespines!(ax_tl, :l, :r, :t)
    xlims!(ax_tl, 0, total_duration)
    ylims!(ax_tl, -0.2, 1.2)
    rowsize!(fig.layout, 3, Fixed(120))

    # Custom x-ticks showing Denver local time
    denver_base = model.start_time - Hour(7)  # UTC → MST
    tick_interval = total_duration > 1440 ? 360.0 : (total_duration > 360 ? 120.0 : 30.0)
    tick_pos = collect(0.0:tick_interval:total_duration)
    prev_date = nothing
    tl_labels = String[]
    for tp in tick_pos
        dt = denver_base + Minute(round(Int, tp))
        d = Dates.Date(dt)
        h = Dates.hour(dt)
        ampm = h < 12 ? "AM" : "PM"
        h12 = h == 0 ? 12 : (h > 12 ? h - 12 : h)
        time_str = "$(h12) $(ampm)"
        if d != prev_date
            push!(tl_labels, time_str * "\n" * Dates.format(d, "m/d"))
            prev_date = d
        else
            push!(tl_labels, time_str)
        end
    end
    ax_tl.xticks = (tick_pos, tl_labels)

    path = joinpath(output_dir, "$(prefix)_growth.gif")
    n_frames = length(snapshots)
    gif_t_start = time()
    gif_last_print = gif_t_start
    record(fig, path, eachindex(snapshots); framerate = 12) do idx
        t_min, φ = snapshots[idx]

        # Progress meter for GIF rendering
        now = time()
        if now - gif_last_print >= 2.0 || idx == n_frames
            pct = idx / n_frames * 100
            elapsed = now - gif_t_start
            eta = pct > 0 ? elapsed / pct * (100 - pct) : NaN
            print("\r  Rendering GIF: $(idx)/$(n_frames) frames ($(round(pct, digits=1))%) | " *
                  "elapsed=$(round(elapsed, digits=0))s | ETA=$(round(eta, digits=0))s   ")
            gif_last_print = now
        end

        # --- Subplot 1: Fire growth ---
        empty!(ax1)
        burned = Float64.(cumulative_burned[idx])
        heatmap!(ax1, lons, lats, burned'; colormap = [:white, (:red, 0.5)], colorrange = (0, 1))
        contour!(ax1, lons, lats, φ'; levels = [0.0], color = :black, linewidth = 2.5)
        poly!(ax1, perim.geometry; color = :transparent, strokecolor = :red, strokewidth = 2)
        scatter!(ax1, [model.lon0], [model.lat0]; color = :orange, markersize = 12, marker = :star5)
        ax1.title = "Fire Growth — t = $(round(t_min, digits=1)) min"

        # --- Subplot 2: Wind ---
        empty!(ax2)
        wf_frame = wind_frames[idx]
        heatmap!(ax2, collect(wind_lons), collect(wind_lats), wf_frame.mag';
                 colormap = :viridis, colorrange = (0, mag_max))
        arrows2d!(ax2, [lo for lo in wind_lons for _ in wind_lats],
                      [la for _ in wind_lons for la in wind_lats],
                      vec(wf_frame.u), vec(wf_frame.v);
                      lengthscale = 0.0005, color = :white, shaftwidth = 1.5, tipwidth = 5.0, tiplength = 5.0)
        heatmap!(ax2, lons, lats, burned'; colormap = [(:white, 0.0), (:red, 0.5)], colorrange = (0, 1))
        contour!(ax2, lons, lats, φ'; levels = [0.0], color = :black, linewidth = 2.5)
        poly!(ax2, perim.geometry; color = :transparent, strokecolor = :red, strokewidth = 2)
        scatter!(ax2, [model.lon0], [model.lat0]; color = :orange, markersize = 12, marker = :star5)
        dt_str = Dates.format(model.start_time + Second(round(Int, t_min * 60)), "HH:MM") * " UTC"
        ax2.title = "Wind — $dt_str"

        # --- Subplot 3: Terrain & Buildings (redraw each frame for fire contour) ---
        empty!(ax3)
        heatmap!(ax3, elev_xs, elev_ys, elev_data'; colormap = :terrain, colorrange = (Float64(minimum(elev_data)), Float64(maximum(elev_data))))
        poly!(ax3, buildings.geometry; color = (:black, 0.8), strokecolor = :black, strokewidth = 0.3)
        heatmap!(ax3, lons, lats, burned'; colormap = [(:white, 0.0), (:red, 0.5)], colorrange = (0, 1))
        contour!(ax3, lons, lats, φ'; levels = [0.0], color = :black, linewidth = 2.5)
        poly!(ax3, perim.geometry; color = :transparent, strokecolor = :red, strokewidth = 2)
        scatter!(ax3, [model.lon0], [model.lat0]; color = :orange, markersize = 12, marker = :star5)
        ax3.title = "Terrain & Buildings"

        # --- Timeline ---
        empty!(ax_tl)
        # Background bar
        poly!(ax_tl, Rect(0, 0.3, total_duration, 0.4); color = :gray90)
        # Progress fill
        if t_min > 0
            poly!(ax_tl, Rect(0, 0.3, t_min, 0.4); color = (:orangered, 0.7))
        end
        # Current position marker
        vlines!(ax_tl, [t_min]; color = :red, linewidth = 2)
        scatter!(ax_tl, [t_min], [0.75]; marker = :dtriangle, markersize = 10, color = :red)
        # Event markers
        for (i, evt) in enumerate(visible_events)
            past = evt.t <= t_min
            lcolor = past ? (:black, 0.7) : (:gray60, 0.5)
            lstyle = past ? :solid : :dash
            vlines!(ax_tl, [evt.t]; color = lcolor, linewidth = 1, linestyle = lstyle)
            # Alternate labels above/below bar to reduce overlap
            y = isodd(i) ? 0.95 : 0.05
            va = isodd(i) ? :bottom : :top
            tcolor = past ? :black : :gray50
            text!(ax_tl, evt.t, y; text = evt.label, fontsize = 9,
                  align = (:center, va), color = tcolor)
        end
        # Current Denver time label
        denver_now = denver_base + Minute(round(Int, t_min))
        denver_str = Dates.format(denver_now, "HH:MM") * " MST"
        text!(ax_tl, t_min, 0.5; text = denver_str, fontsize = 10, color = :white,
              align = (:right, :center), offset = (-5, 0))
    end
    println()  # newline after progress meter
    println("  Saved: $path")
end

#=
Marshall Fire Level Set Simulation with Realistic Data

This example demonstrates fire spread simulation using the level set method with:
- Real HRRR wind data from RapidRefreshData.jl
- Real terrain (slope, aspect) and fuel data from Landfire.jl
- The Rothermel fire spread model

Simulates the first 2 hours of the Marshall Fire (December 30, 2021).
=#

using MarshallWildfire
using MarshallWildfire: GeoSurrogates
using GLMakie
using Dates

#-----------------------------------------------------------------------------# Load Data

println("Loading data (this may take a while on first run)...")
println("  Loading HRRR wind data...")
wind_data = MarshallWildfire.get_hrrr_data()
println("  Loading Landfire data...")
landfire_data = MarshallWildfire.get_landfire_data()
println("Data loaded!")

# Create the fire model
model = MarshallWildfire.FireModel(; wind_data, landfire_data)

#-----------------------------------------------------------------------------# Parameters

# Simulation parameters
const DURATION_HOURS = 2.0
const DURATION_SECONDS = DURATION_HOURS * 3600.0
const GRID_SIZE = 80  # Grid resolution

# Fuel moisture (Marshall Fire had very dry conditions)
const FUEL_MOISTURE = 0.03  # 3% moisture content

# Get domain extent and ignition point
domain = model.domain
ignition = MarshallWildfire.ignition_point

#-----------------------------------------------------------------------------# Precompute Static Fields

println("\nPrecomputing terrain and fuel data on grid...")

# Create level set grid
ls = MarshallWildfire.LevelSet(domain; nx=GRID_SIZE, ny=GRID_SIZE)
xs, ys = ls.xs, ls.ys
nx, ny = length(xs), length(ys)

# Precompute static fields (terrain and fuel don't change with time)
R0_grid = zeros(nx, ny)      # Base spread rate
ϕw_coeff = zeros(nx, ny)     # Wind factor coefficient (depends on fuel)
ϕs_grid = zeros(nx, ny)      # Slope factor
uphill_x = zeros(nx, ny)     # Uphill direction x-component
uphill_y = zeros(nx, ny)     # Uphill direction y-component

for i in 1:nx, j in 1:ny
    x, y = xs[i], ys[j]

    # Get fuel properties
    fuel = MarshallWildfire.get_fuel(model, x, y)

    # Base spread rate
    R0 = MarshallWildfire.no_wind_no_slope_rate(fuel, FUEL_MOISTURE)
    R0_grid[i, j] = isfinite(R0) ? R0 : 0.0

    # Wind factor coefficient (we'll multiply by actual wind speed later)
    # Store the fuel's wind sensitivity
    ϕw_coeff[i, j] = MarshallWildfire.wind_factor(fuel, 1.0)  # factor per m/s

    # Get terrain
    slope_rad, aspect_rad = MarshallWildfire.get_terrain(model, x, y)
    if !isfinite(slope_rad) || !isfinite(aspect_rad)
        slope_rad, aspect_rad = 0.0, 0.0
    end

    # Slope factor
    ϕs = MarshallWildfire.slope_factor(fuel, slope_rad)
    ϕs_grid[i, j] = isfinite(ϕs) ? ϕs : 0.0

    # Uphill direction (opposite of aspect)
    uphill_x[i, j] = -cos(aspect_rad)
    uphill_y[i, j] = -sin(aspect_rad)
end

println("Precomputation complete!")

#-----------------------------------------------------------------------------# Simulation with Precomputed Fields

function run_simulation!(ls, wind_field, start_time, duration_seconds,
                         R0_grid, ϕw_coeff, ϕs_grid, uphill_x, uphill_y,
                         xs, ys, nx, ny)
    # Time stepping parameters
    dt = 10.0  # 10 second time steps
    save_interval = 600.0  # Save every 10 minutes
    reinit_interval = 300.0  # Reinitialize every 5 minutes

    t = 0.0
    last_save = 0.0
    last_reinit = 0.0
    snapshots = [(0.0, deepcopy(ls))]

    step_count = 0

    while t < duration_seconds
        step_count += 1
        if step_count % 100 == 0
            println("  t = $(round(t/60, digits=1)) min ($(round(100*t/duration_seconds, digits=1))%)")
        end

        # Get current datetime for wind lookup
        current_dt = start_time + Second(round(Int, t))

        # Update level set using precomputed fields
        ψ = ls.ψ
        dx, dy = ls.dx, ls.dy
        ψ_new = copy(ψ)

        for i in 2:nx-1, j in 2:ny-1
            x, y = xs[i], ys[j]

            if !isfinite(ψ[i, j])
                continue
            end

            # Get gradient for normal direction
            grad_i = (ψ[i+1, j] - ψ[i-1, j]) / (2dx)
            grad_j = (ψ[i, j+1] - ψ[i, j-1]) / (2dy)
            grad_mag = sqrt(grad_i^2 + grad_j^2)

            if grad_mag < 1e-10
                continue
            end

            # Normal direction (no swap needed - geographic coordinates match matrix storage)
            # grad_i = ∂ψ/∂x (longitude), grad_j = ∂ψ/∂y (latitude)
            norm_x = grad_i / grad_mag
            norm_y = grad_j / grad_mag

            # Get wind at this location and time
            wind = GeoSurrogates.predict(wind_field, (x, y, current_dt))
            u, v = wind.u, wind.v

            if !isfinite(u) || !isfinite(v)
                u, v = 0.0, 0.0
            end

            wind_speed = sqrt(u^2 + v^2)

            # Wind direction (normalized)
            if wind_speed > 0.1
                wind_dir_x = u / wind_speed
                wind_dir_y = v / wind_speed
            else
                wind_dir_x, wind_dir_y = 0.0, 0.0
            end

            # Alignment factors
            wind_alignment = wind_dir_x * norm_x + wind_dir_y * norm_y
            slope_alignment = uphill_x[i, j] * norm_x + uphill_y[i, j] * norm_y

            # Wind factor (scales with wind speed)
            ϕw = ϕw_coeff[i, j] * wind_speed

            # Compute speed
            R0 = R0_grid[i, j]
            ϕs = ϕs_grid[i, j]

            S = R0 * (1.0 + ϕw * max(0.0, wind_alignment) + ϕs * max(0.0, slope_alignment))
            S = min(S, 50.0)  # Clamp to reasonable max (m/s)

            # Convert from m/s to degrees/s
            # At 40°N: 1° longitude ≈ 85 km, 1° latitude ≈ 111 km
            # Use average: ~100 km per degree
            S = S / 100_000.0  # m/s → degrees/s

            if S <= 0 || !isfinite(S)
                continue
            end

            # Upwind gradient for stability
            Dxm = (ψ[i, j] - ψ[i-1, j]) / dx
            Dxp = (ψ[i+1, j] - ψ[i, j]) / dx
            Dym = (ψ[i, j] - ψ[i, j-1]) / dy
            Dyp = (ψ[i, j+1] - ψ[i, j]) / dy

            Dx = max(max(Dxm, 0.0), -min(Dxp, 0.0))
            Dy = max(max(Dym, 0.0), -min(Dyp, 0.0))
            grad_upwind = sqrt(Dx^2 + Dy^2)

            if !isfinite(grad_upwind)
                continue
            end

            ψ_new[i, j] = ψ[i, j] - dt * S * grad_upwind
        end

        ls.ψ .= ψ_new
        t += dt

        # Reinitialize periodically
        if t - last_reinit >= reinit_interval
            MarshallWildfire.reinitialize!(ls; iterations=5)
            last_reinit = t
        end

        # Save snapshot
        if t - last_save >= save_interval
            push!(snapshots, (t, deepcopy(ls)))
            last_save = t
        end
    end

    # Save final state
    push!(snapshots, (t, deepcopy(ls)))

    return snapshots
end

println("\nRunning simulation for $(DURATION_HOURS) hours...")
println("  Grid size: $(GRID_SIZE)×$(GRID_SIZE)")
println("  Fuel moisture: $(FUEL_MOISTURE * 100)%")

# Initialize level set with small fire at ignition point
# 0.001 degrees ≈ 100m radius at this latitude
MarshallWildfire.init_circle!(ls, ignition.lon, ignition.lat, 0.001)

println("  Ignition point: ($(ignition.lon), $(ignition.lat))")

# Run the simulation
snapshots = run_simulation!(ls, model.wind_field, model.start_time, DURATION_SECONDS,
                            R0_grid, ϕw_coeff, ϕs_grid, uphill_x, uphill_y,
                            xs, ys, nx, ny)

println("Simulation complete! $(length(snapshots)) snapshots saved.")

# Print wind info at ignition point
println("\nWind conditions at ignition point:")
for t_check in [0.0, 3600.0, 7200.0]
    dt_check = model.start_time + Second(round(Int, t_check))
    wind = GeoSurrogates.predict(model.wind_field, (ignition.lon, ignition.lat, dt_check))
    speed = sqrt(wind.u^2 + wind.v^2)
    dir_deg = rad2deg(atan(wind.v, wind.u))
    println("  t = $(t_check/3600) hours: wind = $(round(speed, digits=1)) m/s from $(round(dir_deg, digits=0))°")
end

#-----------------------------------------------------------------------------# Load Perimeter for Comparison

println("\nLoading observed fire perimeter...")
perimeter = MarshallWildfire.get_perimeter()

#-----------------------------------------------------------------------------# Visualization

println("Creating visualization...")

fig = Figure(size = (1200, 500))

# Plot final simulated fire
ax1 = Axis(fig[1, 1],
    title = "Simulated Fire (t = $(DURATION_HOURS) hours)",
    xlabel = "Longitude",
    ylabel = "Latitude",
    aspect = DataAspect()
)

# Get final level set
_, ls_final = snapshots[end]

# Plot the level set field
heatmap!(ax1, ls_final.xs, ls_final.ys, ls_final.ψ';
    colormap = :RdYlBu,
    colorrange = (-0.02, 0.02)
)

# Plot the simulated fire front (don't transpose for contour!)
contour!(ax1, ls_final.xs, ls_final.ys, ls_final.ψ;
    levels = [0.0],
    color = :black,
    linewidth = 2
)

# Plot observed fire perimeter
poly!(ax1, perimeter.geometry;
    color = :transparent,
    strokecolor = :red,
    strokewidth = 2
)

# Plot ignition point
scatter!(ax1, [ignition.lon], [ignition.lat];
    color = :orange,
    markersize = 15,
    marker = :star5
)

# Plot time evolution of fire front
ax2 = Axis(fig[1, 2],
    title = "Fire Spread Over Time",
    xlabel = "Longitude",
    ylabel = "Latitude",
    aspect = DataAspect()
)

# Use a color gradient for different times
n_snaps = length(snapshots)
colors = cgrad(:inferno, n_snaps, categorical=true)

for (i, (t_snap, snap)) in enumerate(snapshots)
    contour!(ax2, snap.xs, snap.ys, snap.ψ;
        levels = [0.0],
        color = colors[i],
        linewidth = 1.5
    )
end

# Plot observed fire perimeter for reference
poly!(ax2, perimeter.geometry;
    color = :transparent,
    strokecolor = :red,
    strokewidth = 2,
    linestyle = :dash
)

# Add ignition point
scatter!(ax2, [ignition.lon], [ignition.lat];
    color = :lime,
    markersize = 15,
    marker = :star5
)

# Add colorbar for time
Colorbar(fig[1, 3], limits = (0, DURATION_HOURS), colormap = :inferno,
    label = "Time (hours)")

# Add overall title
Label(fig[0, :], "Marshall Fire Simulation: First $(DURATION_HOURS) Hours (HRRR Wind + Landfire Terrain)",
    fontsize = 18)

save("examples/marshall_fire_realistic.png", fig)
println("\nFigure saved to examples/marshall_fire_realistic.png")
display(fig)

#-----------------------------------------------------------------------------# Print Statistics

println("\n" * "="^60)
println("Simulation Statistics")
println("="^60)
for (i, (t_snap, snap)) in enumerate(snapshots)
    hours = t_snap / 3600
    area = MarshallWildfire.burned_area(snap)
    # Convert from degrees² to km² (approximate at this latitude)
    # 1 degree longitude ≈ 85 km, 1 degree latitude ≈ 111 km at 40°N
    area_km2 = area * 85 * 111
    println("  t = $(round(hours, digits=2)) hours: burned area ≈ $(round(area_km2, digits=2)) km²")
end

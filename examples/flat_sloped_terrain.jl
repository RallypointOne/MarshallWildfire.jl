#=
Toy Example: Level Set Fire Spread on Sloped Terrain (No Wind)

This example demonstrates the level set method for tracking fire front propagation
using the Rothermel fire spread model with:
- Sloped terrain (4 different slope directions)
- No wind
- Uniform short grass fuel (Anderson Model 1)

The fire starts as a small circle and expands based on the local spread rate.
Fire spreads faster uphill due to preheating of fuels above the flame.
=#

using MarshallWildfire
using Extents: Extent
using GLMakie

#-----------------------------------------------------------------------------# Parameters

# Domain: 1km x 1km
const DOMAIN = Extent(X = (0.0f0, 1000.0f0), Y = (0.0f0, 1000.0f0))

# No wind
const WIND_SPEED = 0.0

# Slope angle (15 degrees)
const SLOPE_ANGLE = deg2rad(15.0)

# Use Anderson Fuel Model 1: Short grass
const FUEL = MarshallWildfire.FUEL_MODELS[1]

# Fuel moisture (5% = 0.05)
const MOISTURE = 0.05

# Simulation duration (seconds)
const DURATION = 120.0

#-----------------------------------------------------------------------------# Simulation Function

"""
    run_simulation(uphill_dir_rad)

Run a fire spread simulation with terrain sloping such that uphill is in the
given direction (radians). Returns the final level set state.

The speed function is direction-dependent: fire spreads faster when moving
uphill due to preheating of fuels above the flame.
"""
function run_simulation(uphill_dir_rad)
    # Create fresh level set
    ls = MarshallWildfire.LevelSet(DOMAIN; nx=100, ny=100)
    MarshallWildfire.init_circle!(ls, 500.0, 500.0, 50.0)

    # Get Rothermel base spread rate (no wind, no slope)
    result = MarshallWildfire.rothermel_spread_rate(
        FUEL,
        0.0,  # no wind
        0.0,  # wind direction (irrelevant)
        0.0,  # no slope for base rate
        0.0;  # aspect (irrelevant)
        Mf = MOISTURE
    )
    R0 = result.R0  # Base spread rate

    # Slope factor from Rothermel
    slope_factor = MarshallWildfire.slope_factor(FUEL, SLOPE_ANGLE)

    # Direction-dependent speed function
    # Speed increases when the fire front normal aligns with uphill direction
    function directional_speed(x, y, t, nx, ny)
        # Uphill direction
        uphill_x = cos(uphill_dir_rad)
        uphill_y = sin(uphill_dir_rad)

        # (nx, ny) is the outward normal of the fire front
        # Fire spreads faster when normal aligns with uphill direction
        # Dot product: how aligned is the spread direction with uphill?
        alignment = uphill_x * nx + uphill_y * ny

        # Speed: base rate + slope contribution (only when spreading uphill)
        return R0 * (1.0 + slope_factor * max(0.0, alignment))
    end

    # Maximum speed (when perfectly aligned with uphill)
    max_speed = R0 * (1.0 + slope_factor)

    # Run simulation with directional speed
    snapshots = MarshallWildfire.simulate_directional(ls, directional_speed, DURATION;
        max_speed = max_speed,
        save_interval = DURATION,  # Only save final
        reinit_interval = 30.0,
        show_progress = false
    )

    return snapshots[end][2]  # Return final level set
end

#-----------------------------------------------------------------------------# Run 4 Simulations

# Uphill directions: East, North, West, South
# (Fire spreads faster in these directions due to slope)
slope_directions = [
    (0.0, "Uphill: East", (1, 0)),
    (π/2, "Uphill: North", (0, 1)),
    (π, "Uphill: West", (-1, 0)),
    (3π/2, "Uphill: South", (0, -1))
]

println("Running 4 simulations with different slope directions...")
println("Slope angle: $(rad2deg(SLOPE_ANGLE))°, No wind")
results = []
for (dir, name, _) in slope_directions
    println("  $name...")
    ls_final = run_simulation(dir)
    push!(results, (name, ls_final))
end
println("Done!")

#-----------------------------------------------------------------------------# Visualization (2x2 grid)

fig = Figure(size = (1000, 1000))

for (i, ((name, ls), (_, _, arrow_dir))) in enumerate(zip(results, slope_directions))
    row = (i - 1) ÷ 2 + 1
    col = (i - 1) % 2 + 1

    ax = Axis(fig[row, col],
        title = name,
        xlabel = "x (m)",
        ylabel = "y (m)",
        aspect = 1
    )

    # Plot level set field and fire front
    heatmap!(ax, ls.xs, ls.ys, ls.ψ'; colormap = :RdYlBu)
    contour!(ax, ls.xs, ls.ys, ls.ψ'; levels = [0.0], color = :black, linewidth = 2)

    # Add ignition point marker
    scatter!(ax, [500.0], [500.0]; color = :orange, markersize = 12, marker = :star5)

    # Add uphill arrow (showing direction fire spreads faster)
    ax_x, ay_y = arrow_dir
    arrows2d!(ax, [150.0], [850.0], [80.0 * ax_x], [80.0 * ay_y];
        color = :forestgreen, shaftwidth = 2)
end

# Add overall title
Label(fig[0, :], "Fire Spread on Sloped Terrain ($(round(Int, rad2deg(SLOPE_ANGLE)))° slope, no wind, t = $(DURATION)s)",
    fontsize = 18)

save("examples/flat_sloped_terrain.png", fig)
println("\nFigure saved to examples/flat_sloped_terrain.png")
display(fig)

using MarshallWildfire, GeoSurrogates, GLMakie, GeoMakie, Rasters, ArchGDAL

#-----------------------------------------------------------------------------# Run fire model on Marshall Fire
# Create the fire model with wind and landfire data
model = MarshallWildfire.FireModel()

# Simulate for 6 hours (the fire spread rapidly due to extreme winds)
duration = 1 * 60 * 60  # 6 hours in seconds

# Run simulation with 60-second time steps on a 150x150 grid
# Fuel moisture ~5% (very dry conditions on Dec 30, 2021)
results = MarshallWildfire.simulate(model, duration; dt=60.0, nx=150, ny=150, Mf=0.05)

@info "Simulation complete: $(length(results)) snapshots"

#-----------------------------------------------------------------------------# Visualize results
# Get fire perimeter for comparison
perim = MarshallWildfire.get_perimeter()

fig = Figure(size=(1000, 800))
ax = GeoAxis(fig[1, 1]; dest="+proj=webmerc", title="Marshall Fire Simulation vs Actual Perimeter")

# Plot actual fire perimeter
# poly!(ax, perim.geometry; color=(:red, 0.2), strokecolor=:red, strokewidth=2, label="Actual perimeter")

# Plot simulated fire spread at different times
subsampled = results[1:10:end]
colors = cgrad(:YlOrRd, length(subsampled))
for (i, (t, ls)) in enumerate(subsampled)
    hours = t / 3600
    contour!(ax, ls.xs, ls.ys, ls.ψ; levels=[0.0], color=colors[i], linewidth=1.5)
end

# Plot ignition point
# scatter!(ax, [MarshallWildfire.ignition_point.lon], [MarshallWildfire.ignition_point.lat];
#          marker=:star5, markersize=20, color=:white, strokecolor=:black, strokewidth=1)

# Final simulated perimeter
final_t, final_ls = results[end]
contour!(ax, final_ls.xs, final_ls.ys, final_ls.ψ; levels=[0.0], color=:blue, linewidth=3, label="Simulated ($(round(final_t/3600, digits=1))h)")

fig

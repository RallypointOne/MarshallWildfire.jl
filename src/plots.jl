#-----------------------------------------------------------------------------# plot_marshall_perimeter
function plot_marshall_perimeter(; output_dir = joinpath(@__DIR__, "..", "output"))
    mkpath(output_dir)
    perim = get_perimeter()
    fig = Figure(size = (800, 800))
    ax = GeoAxis(fig[1, 1]; dest = "+proj=webmerc", title = "Marshall Fire Perimeter")
    hidedecorations!(ax, label=false, ticklabels=false, ticks=false, grid=true)
    m = Tyler.Map(extent2; figure = fig, axis = ax)
    wait(m)  # Wait for tiles to load

    # Plot the perimeter on top of the map (GeoAxis handles the projection)
    poly!(ax, perim.geometry; color = (:red, 0.3), strokecolor = :red, strokewidth = 2)

    # Plot the ignition point as a star (translate in z to render on top)
    p = scatter!(ax, [ignition_point.lon], [ignition_point.lat]; marker = :star5, markersize = 20, color = :black, strokecolor = :black, strokewidth = 0)
    translate!(p, 0, 0, 1)

    # Add label for ignition point
    t = text!(ax, ignition_point.lon, ignition_point.lat; text = "Ignition Point", align = (:left, :bottom), offset = (10, 5), fontsize = 14, color = :black, strokecolor = :black, strokewidth = 0)
    translate!(t, 0, 0, 1)

    # Save the figure
    output_path = joinpath(output_dir, "marshall_fire_perimeter.png")
    save(output_path, fig)
    return
end

#-----------------------------------------------------------------------------# plot_perimeter_with_buildings
function plot_perimeter_with_buildings(; output_dir = joinpath(@__DIR__, "..", "output"))
    mkpath(output_dir)
    perim = get_perimeter()
    buildings = get_building_footprints()
    powerlines = get_power_lines()

    fig = Figure(size = (800, 800))
    ax = GeoAxis(fig[1, 1]; dest = "+proj=webmerc", title = "Marshall Fire: Buildings & Power Lines")
    hidedecorations!(ax, label=false, ticklabels=false, ticks=false, grid=true)
    m = Tyler.Map(extent2; figure = fig, axis = ax)
    wait(m)

    # Plot building footprints in black
    poly!(ax, buildings.geometry; color = :black, strokecolor = :black, strokewidth = 0.5, alpha=0.7)

    # Plot power lines in purple (dotted)
    lines!(ax, powerlines.geometry; color = :purple, linewidth = 2, linestyle = :dot)

    # Plot the perimeter on top
    poly!(ax, perim.geometry; color = (:red, 0.3), strokecolor = :red, strokewidth = 2)

    # Plot the ignition point as a star
    p = scatter!(ax, [ignition_point.lon], [ignition_point.lat]; marker = :star5, markersize = 20, color = :yellow, strokecolor = :black, strokewidth = 1)
    translate!(p, 0, 0, 1)

    # Add legend
    legend_elements = [
        PolyElement(color = :black, strokecolor = :black, strokewidth = 0.5),
        LineElement(color = :purple, linewidth = 2, linestyle = :dot),
        PolyElement(color = (:red, 0.3), strokecolor = :red, strokewidth = 2),
        MarkerElement(marker = :star5, color = :yellow, strokecolor = :black, strokewidth=1, markersize = 15)
    ]
    legend_labels = ["Buildings", "Power Lines", "Fire Perimeter", "Ignition Point"]
    Legend(fig[1, 2], legend_elements, legend_labels; framevisible = true, padding = (10, 10, 10, 10))

    # Save the figure
    output_path = joinpath(output_dir, "marshall_fire_with_buildings.png")
    save(output_path, fig)
    return fig
end

#-----------------------------------------------------------------------------# plot_damage_assessment
"""
    plot_damage_assessment(; output_dir)

Plot fire perimeter with damage assessment points from Boulder County's official data.
- Destroyed: Red
- Major damage: Orange
- Minor damage: Yellow
- Affected: Blue
"""
function plot_damage_assessment(; output_dir = joinpath(@__DIR__, "..", "output"))
    mkpath(output_dir)

    # Load data
    perim = get_perimeter()
    residential = get_residential_damage_assessment()
    commercial = get_commercial_damage_assessment()

    # Combine residential and commercial
    all_points = vcat(collect(residential), collect(commercial))

    # Group by damage category
    destroyed = filter(f -> get(f.properties, :damage, "") == "Destroyed", all_points)
    major = filter(f -> get(f.properties, :damage, "") == "Major", all_points)
    minor = filter(f -> get(f.properties, :damage, "") == "Minor", all_points)
    affected = filter(f -> get(f.properties, :damage, "") == "Affected", all_points)

    @info "Loaded damage assessment: $(length(destroyed)) destroyed, $(length(major)) major, $(length(minor)) minor, $(length(affected)) affected"

    fig = Figure(size = (900, 800))
    ax = GeoAxis(fig[1, 1]; dest = "+proj=webmerc", title = "Marshall Fire Damage Assessment")
    hidedecorations!(ax, label=false, ticklabels=false, ticks=false, grid=true)
    m = Tyler.Map(extent2; figure = fig, axis = ax)
    wait(m)

    # Plot building footprints as background layer
    buildings = get_building_footprints()
    poly!(ax, buildings.geometry; color = (:gray, 0.5), strokecolor = (:gray, 0.7), strokewidth = 0.3)

    # Helper to extract coordinates from point features
    get_coords(features) = begin
        xs = [GI.x(f.geometry) for f in features]
        ys = [GI.y(f.geometry) for f in features]
        (xs, ys)
    end

    # Plot points by damage category (least severe first so severe are on top)
    # Use translate! to ensure points render above building footprints
    if !isempty(affected)
        xs, ys = get_coords(affected)
        p = scatter!(ax, xs, ys; color = :dodgerblue, markersize = 5, strokewidth = 0.5, strokecolor = :black)
        translate!(p, 0, 0, 1)
    end

    if !isempty(minor)
        xs, ys = get_coords(minor)
        p = scatter!(ax, xs, ys; color = :gold, markersize = 5, strokewidth = 0.5, strokecolor = :black)
        translate!(p, 0, 0, 2)
    end

    if !isempty(major)
        xs, ys = get_coords(major)
        p = scatter!(ax, xs, ys; color = :darkorange, markersize = 6, strokewidth = 0.5, strokecolor = :black)
        translate!(p, 0, 0, 3)
    end

    if !isempty(destroyed)
        xs, ys = get_coords(destroyed)
        p = scatter!(ax, xs, ys; color = :red, markersize = 6, strokewidth = 0.5, strokecolor = :black)
        translate!(p, 0, 0, 4)
    end

    # Plot the fire perimeter outline
    poly!(ax, perim.geometry; color = :transparent, strokecolor = :black, strokewidth = 2)

    # Plot the ignition point as a star
    p = scatter!(ax, [ignition_point.lon], [ignition_point.lat];
                 marker = :star5, markersize = 20, color = :white, strokecolor = :black, strokewidth = 1)
    translate!(p, 0, 0, 1)

    # Add legend
    legend_elements = [
        MarkerElement(marker = :circle, color = :red, markersize = 10),
        MarkerElement(marker = :circle, color = :darkorange, markersize = 10),
        MarkerElement(marker = :circle, color = :gold, markersize = 10),
        MarkerElement(marker = :circle, color = :dodgerblue, markersize = 10),
        PolyElement(color = (:gray, 0.5), strokecolor = (:gray, 0.7), strokewidth = 0.3),
        PolyElement(color = :transparent, strokecolor = :black, strokewidth = 2),
        MarkerElement(marker = :star5, color = :white, strokecolor = :black, strokewidth = 1, markersize = 15)
    ]
    legend_labels = [
        "Destroyed ($(length(destroyed)))",
        "Major ($(length(major)))",
        "Minor ($(length(minor)))",
        "Affected ($(length(affected)))",
        "Buildings",
        "Fire Perimeter",
        "Ignition Point"
    ]
    Legend(fig[1, 2], legend_elements, legend_labels; framevisible = true, padding = (10, 10, 10, 10))

    # Save the figure
    output_path = joinpath(output_dir, "damage_assessment.png")
    save(output_path, fig)
    return fig
end

#-----------------------------------------------------------------------------# plot_landfire_layers
function plot_landfire_layers(; output_dir = joinpath(@__DIR__, "..", "output"))
    mkpath(output_dir)
    nt = get_landfire_data()

    for (k, v) in pairs(nt)
        @info "Plotting Landfire layer: $k"
        fig = Figure(size = (800, 800))
        ax = GeoAxis(fig[1, 1]; dest = "+proj=webmerc", title = string(k))
        hidedecorations!(ax, label=false, ticklabels=false, ticks=false, grid=true)
        hm = heatmap!(ax, v)
        Colorbar(fig[1, 2], hm)

        # Save each layer to its own file
        output_path = joinpath(output_dir, "landfire_$k.png")
        save(output_path, fig)
    end
    return
end

#-----------------------------------------------------------------------------# plot_fuel_flammability
"""
    plot_fuel_flammability(; output_dir)

Plot Landfire fuel model data colored by flammability.
Anderson 13 fuel models ranked from most to least flammable based on
typical fire behavior (spread rate and intensity).
"""
function plot_fuel_flammability(; output_dir = joinpath(@__DIR__, "..", "output"))
    mkpath(output_dir)
    lf = get_landfire_data()

    # Find fuel model layer (names include US_250 prefix)
    fuel_key = haskey(lf, :US_250FBFM13) ? :US_250FBFM13 : (haskey(lf, :US_250FBFM40) ? :US_250FBFM40 : nothing)
    if isnothing(fuel_key)
        error("No fuel model layer found in Landfire data. Available layers: $(keys(lf))")
    end
    fuel = lf[fuel_key]

    # Anderson 13 fuel models ranked by flammability (spread rate × intensity)
    # Higher rank = more flammable
    # Based on typical fire behavior characteristics
    flammability_rank = Dict(
        1  => 8,   # Short grass - fast spread, moderate intensity
        2  => 9,   # Timber grass - fast spread, higher intensity
        3  => 10,  # Tall grass - very fast spread, high intensity
        4  => 13,  # Chaparral - extreme fire behavior
        5  => 6,   # Brush - moderate
        6  => 5,   # Dormant brush - moderate
        7  => 4,   # Southern rough - lower (high moisture)
        8  => 2,   # Compact timber litter - slow spread
        9  => 3,   # Hardwood litter - slow/moderate spread
        10 => 7,   # Timber understory - moderate
        11 => 11,  # Light logging slash - high
        12 => 12,  # Medium logging slash - very high
        13 => 14,  # Heavy logging slash - extreme
        # Non-burnable codes (91, 92, 93, 98, 99) -> 0
    )

    # Convert fuel codes to flammability values
    flammability = map(fuel.data) do code
        c = round(Int, code)
        get(flammability_rank, c, 0)
    end
    flammability_raster = Raster(flammability; dims=dims(fuel))

    fig = Figure(size = (900, 800))
    ax = GeoAxis(fig[1, 1]; dest = "+proj=webmerc", title = "Fuel Flammability")
    hidedecorations!(ax, label=false, ticklabels=false, ticks=false, grid=true)

    # Use a fire-themed colormap (yellow -> orange -> red -> dark red)
    hm = heatmap!(ax, flammability_raster; colormap = :YlOrRd, colorrange = (0, 14))

    # Custom colorbar with fuel model labels
    cb = Colorbar(fig[1, 2], hm;
        label = "Flammability",
        ticks = ([1, 4, 7, 10, 13], ["Low", "Moderate", "High", "Very High", "Extreme"])
    )

    # Add fire perimeter overlay
    perim = get_perimeter()
    poly!(ax, perim.geometry; color = :transparent, strokecolor = :black, strokewidth = 2)

    output_path = joinpath(output_dir, "fuel_flammability.png")
    save(output_path, fig)
    return fig
end

#-----------------------------------------------------------------------------# plot_hrrr_wind
function plot_hrrr_wind(; output_dir = joinpath(@__DIR__, "..", "output"))
    mkpath(output_dir)
    wind_data = get_hrrr_data()

    times = lookup(wind_data, Ti)
    xs = lookup(wind_data, X)
    ys = lookup(wind_data, Y)

    # Calculate wind magnitude for colorscale limits
    u_all = wind_data[Band=At(:u)]
    v_all = wind_data[Band=At(:v)]
    mag_all = sqrt.(u_all.^2 .+ v_all.^2)
    clims = (0, maximum(mag_all))

    ext = Extents.extent(wind_data)
    # Nudge limits inward by half a cell to remove whitespace
    dx = abs(step(xs) / 2)
    dy = abs(step(ys) / 2)
    xlims = (ext.X[1] + dx, ext.X[2] - 2dx)
    ylims = (ext.Y[1] + dy, ext.Y[2] - 2dy)
    title_text = Observable(Dates.format(times[1], "yyyy-mm-dd HH:MM") * " UTC")
    fig = Figure(size = (800, 850))

    # Timeline progress bar
    time_nums = datetime2unix.(times)
    current_time = Observable(time_nums[1])
    ax_timeline = Axis(fig[1, 1:2]; height = 30,
                       limits = (extrema(time_nums), (-0.5, 0.5)),
                       title = title_text)
    hideydecorations!(ax_timeline)
    hidexdecorations!(ax_timeline)
    hidespines!(ax_timeline)
    lines!(ax_timeline, time_nums, zeros(length(times)); color = :gray, linewidth = 3)
    scatter!(ax_timeline, current_time, 0; color = :red, markersize = 15)

    ax = GeoAxis(fig[2, 1]; dest = "+proj=webmerc", xlabel = "Longitude", ylabel = "Latitude",
              limits = (xlims, ylims))
    # hidedecorations!(ax, label=false, ticklabels=true, ticks=true, grid=true)

    # Initial data
    u = u_all[Ti=1].data
    v = v_all[Ti=1].data
    mag = sqrt.(u.^2 .+ v.^2)

    hm = heatmap!(ax, collect(xs), collect(ys), mag'; colorrange = clims, colormap = :viridis)
    arrows2d!(ax, collect(xs), collect(ys), u', v';
            lengthscale = 0.003, color = :white)
    cb = Colorbar(fig[2, 2], hm; label = "Wind Speed (m/s)")

    # Make timeline row fixed height, map row fills remaining space
    rowsize!(fig.layout, 1, 50)
    rowsize!(fig.layout, 2, Auto())

    output_path = joinpath(output_dir, "hrrr_wind.gif")
    record(fig, output_path, eachindex(times); framerate = 4) do i
        u = u_all[Ti=i].data
        v = v_all[Ti=i].data
        mag = sqrt.(u.^2 .+ v.^2)

        # Update heatmap
        hm[3] = mag'

        # Update arrows - need to clear and redraw
        empty!(ax)
        heatmap!(ax, collect(xs), collect(ys), mag'; colorrange = clims, colormap = :viridis)
        arrows2d!(ax, collect(xs), collect(ys), u', v';
                lengthscale = 0.003, color = :white)

        title_text[] = Dates.format(times[i], "yyyy-mm-dd HH:MM") * " UTC"
        current_time[] = time_nums[i]
    end

    return output_path
end

#-----------------------------------------------------------------------------# plot_fuel_moisture
"""
    plot_fuel_moisture(; output_dir)

Plot gridMET 100-hour dead fuel moisture on the day of the Marshall Fire (December 30, 2021).
Shows how dry the fuels were when the fire ignited.
"""
function plot_fuel_moisture(; output_dir = joinpath(@__DIR__, "..", "output"))
    mkpath(output_dir)

    # Get fuel moisture data for 2021 - use larger extent to show multiple gridMET pixels (~4km resolution)
    fm_full = get_gridmet_fuel_moisture(; year=2021, variable=:fm100)

    # Extract data for December 30, 2021 (fire day)
    fire_date = Date(2021, 12, 30)
    fm_day = fm_full[Ti=At(fire_date)]

    # Get fire perimeter
    perim = get_perimeter()

    # Get value at ignition point
    fm_at_ignition = fm_day[X=Near(ignition_point.lon), Y=Near(ignition_point.lat)]

    @info "100-hour fuel moisture at ignition point on $(fire_date): $(round(fm_at_ignition, digits=1))%"

    # Use a larger extent to show more gridMET pixels (gridMET is ~4km resolution)
    plot_extent = Extents.grow(extent, 1.0f0)

    fig = Figure(size = (900, 800))
    ax = GeoAxis(fig[1, 1]; dest = "+proj=webmerc",
                 title = "100-Hour Dead Fuel Moisture - December 30, 2021",
                 limits = (plot_extent.X, plot_extent.Y))
    hidedecorations!(ax, label=false, ticklabels=false, ticks=false, grid=true)

    # Plot fuel moisture heatmap (lower = drier = more dangerous)
    # Typical range is 5-30% for dead fuels
    # RdYlBu: red (low/dry) -> yellow -> blue (high/moist)
    hm = heatmap!(ax, fm_day; colormap = :RdYlBu, colorrange = (5, 25))

    # Add fire perimeter overlay
    poly!(ax, perim.geometry; color = :transparent, strokecolor = :black, strokewidth = 2)

    # Plot the ignition point
    p = scatter!(ax, [ignition_point.lon], [ignition_point.lat];
                 marker = :star5, markersize = 20, color = :white, strokecolor = :black, strokewidth = 1)
    translate!(p, 0, 0, 1)

    # Add annotation for ignition point value
    t = text!(ax, ignition_point.lon, ignition_point.lat;
              text = "  $(round(fm_at_ignition, digits=1))%",
              align = (:left, :center), fontsize = 14, color = :black)
    translate!(t, 0, 0, 1)

    # Colorbar
    cb = Colorbar(fig[1, 2], hm;
        label = "Fuel Moisture (%)",
        ticks = ([5, 10, 15, 20, 25], ["5% (Very Dry)", "10%", "15%", "20%", "25% (Moist)"])
    )

    # Save the figure
    output_path = joinpath(output_dir, "fuel_moisture.png")
    save(output_path, fig)
    return fig
end

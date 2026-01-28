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

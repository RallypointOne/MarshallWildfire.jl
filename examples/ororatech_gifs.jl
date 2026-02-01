using GeoJSON
using GLMakie
using GeoMakie
using GeoInterface
using LibGEOS

function create_fire_gif(input_path::String; output_path::String=replace(input_path, ".json" => ".gif"), framerate::Int=4)
    println("Processing: $input_path")

    # Load the GeoJSON data
    data = GeoJSON.read(read(input_path, String))

    # Extract features and convert pps to hours
    features = collect(data)
    pps_values = [f.properties[:pps] for f in features]
    min_pps = minimum(pps_values)

    # Convert pps (seconds) to hours from start
    hours = [floor(Int, (pps - min_pps) / 3600) for pps in pps_values]

    # Create a mapping from hour to features
    hour_to_features = Dict{Int, Vector}()
    for (i, hour) in enumerate(hours)
        if !haskey(hour_to_features, hour)
            hour_to_features[hour] = []
        end
        push!(hour_to_features[hour], features[i])
    end

    # Get sorted unique hours
    unique_hours = sort(collect(keys(hour_to_features)))
    println("  Number of unique hours with data: ", length(unique_hours))
    println("  Hours with data: ", unique_hours)
    println("  Number of features: ", length(features))
    println("  PPS range: ", minimum(pps_values), " - ", maximum(pps_values))

    # Function to get cumulative union up to a given hour
    function get_cumulative_union(up_to_hour)
        geoms = LibGEOS.Geometry[]
        for hour in unique_hours
            hour > up_to_hour && break
            for f in hour_to_features[hour]
                g = GeoInterface.convert(LibGEOS, GeoInterface.geometry(f))
                # Make geometry valid if needed
                if !LibGEOS.isValid(g)
                    g = LibGEOS.makeValid(g)
                end
                push!(geoms, g)
            end
        end

        if isempty(geoms)
            return nothing
        elseif length(geoms) == 1
            return geoms[1]
        else
            result = geoms[1]
            for g in geoms[2:end]
                try
                    result = LibGEOS.union(result, g)
                catch e
                    # If union fails, try making both valid
                    result = LibGEOS.makeValid(result)
                    g = LibGEOS.makeValid(g)
                    result = LibGEOS.union(result, g)
                end
            end
            return result
        end
    end

    # Detect large gaps in hours and start from after the gap if needed
    min_hour = minimum(unique_hours)
    max_hour = maximum(unique_hours)

    # Check for gaps larger than 24 hours
    sorted_hours = sort(unique_hours)
    for i in 1:(length(sorted_hours)-1)
        gap = sorted_hours[i+1] - sorted_hours[i]
        if gap > 24
            min_hour = sorted_hours[i+1]
            println("  Detected gap of $gap hours, starting from hour $min_hour")
        end
    end

    println("  Pre-computing cumulative unions for hours $min_hour to $max_hour...")
    cumulative_unions = Dict{Int, Any}()
    last_geom = nothing
    for hour in min_hour:max_hour
        if hour in unique_hours
            last_geom = get_cumulative_union(hour)
        end
        cumulative_unions[hour] = last_geom
    end
    println("  Done computing unions.")

    # Get bounding box from all features using GeoInterface extent
    ext = GeoInterface.extent(data)
    lon_min, lon_max = ext.X
    lat_min, lat_max = ext.Y

    # Add some padding
    pad = 0.02 * max(lon_max - lon_min, lat_max - lat_min)
    lon_min -= pad
    lon_max += pad
    lat_min -= pad
    lat_max += pad

    # Extract fire name from filename for title
    fire_name = basename(input_path) |> x -> replace(x, ".json" => "")

    # Create the animation
    fig = Figure(size=(800, 800))
    title_obs = Observable{Any}(fire_name)
    ax = GeoAxis(fig[1, 1];
        dest="+proj=longlat",
        limits=(lon_min, lon_max, lat_min, lat_max),
        title=title_obs,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        xticksvisible=false,
        yticksvisible=false
    )

    # Observables for area text
    area_text = Observable("")

    # Conversion factor from deg² to km²
    center_lat = (lat_min + lat_max) / 2
    deg2km² = 111.0^2 * cosd(center_lat)

    # Record the animation - one frame per hour (min_hour to max_hour)
    record(fig, output_path, min_hour:max_hour; framerate=framerate) do hour
        empty!(ax)
        geom = cumulative_unions[hour]
        if geom !== nothing
            # Plot convex hull first (underneath) in yellow
            hull = LibGEOS.convexhull(geom)
            poly!(ax, hull; color=(:yellow, 0.5), strokecolor=:goldenrod, strokewidth=1)
            # Plot fire polygon on top in red
            poly!(ax, geom; color=(:red, 0.7), strokecolor=:darkred, strokewidth=1)

            # Compute areas in km²
            fire_area_km2 = LibGEOS.area(geom) * deg2km²
            hull_area_km2 = LibGEOS.area(hull) * deg2km²

            # Update title with hour and colored area text
            title_obs[] = rich(
                "$fire_name — Hour $hour\n",
                rich("Fire: $(round(fire_area_km2, digits=2)) km²", color=:darkred),
                "  ",
                rich("Hull: $(round(hull_area_km2, digits=2)) km²", color=:goldenrod)
            )
        end
    end

    println("  GIF saved to: $output_path")
    return output_path
end

# Helper function to load and preprocess fire data
function load_fire_data(input_path)
    data = GeoJSON.read(read(input_path, String))
    features = collect(data)
    pps_values = [f.properties[:pps] for f in features]
    min_pps = minimum(pps_values)
    hours = [floor(Int, (pps - min_pps) / 3600) for pps in pps_values]

    hour_to_features = Dict{Int, Vector}()
    for (i, hour) in enumerate(hours)
        if !haskey(hour_to_features, hour)
            hour_to_features[hour] = []
        end
        push!(hour_to_features[hour], features[i])
    end

    unique_hours = sort(collect(keys(hour_to_features)))

    function get_cumulative_union(up_to_hour)
        geoms = LibGEOS.Geometry[]
        for hour in unique_hours
            hour > up_to_hour && break
            for f in hour_to_features[hour]
                g = GeoInterface.convert(LibGEOS, GeoInterface.geometry(f))
                if !LibGEOS.isValid(g)
                    g = LibGEOS.makeValid(g)
                end
                push!(geoms, g)
            end
        end
        isempty(geoms) && return nothing
        length(geoms) == 1 && return geoms[1]
        result = geoms[1]
        for g in geoms[2:end]
            try
                result = LibGEOS.union(result, g)
            catch
                result = LibGEOS.makeValid(result)
                g = LibGEOS.makeValid(g)
                result = LibGEOS.union(result, g)
            end
        end
        return result
    end

    ext = GeoInterface.extent(data)
    fire_name = basename(input_path) |> x -> replace(x, ".json" => "")

    return (; data, unique_hours, get_cumulative_union, ext, fire_name)
end

function create_combined_gif(input_paths::Vector{String}; output_path::String, framerate::Int=4)
    println("Creating combined GIF...")

    # Load all fire data
    fire_data = [load_fire_data(p) for p in input_paths]

    # Use hours 0 to max across all fires
    max_hour = maximum(maximum(fd.unique_hours) for fd in fire_data)
    println("  Max hour across all fires: $max_hour")

    # Pre-compute unions for all fires
    println("  Pre-computing cumulative unions...")
    all_unions = []
    for fd in fire_data
        unions = Dict{Int, Any}()
        last_geom = nothing
        for hour in 0:max_hour
            if hour in fd.unique_hours
                last_geom = fd.get_cumulative_union(hour)
            end
            unions[hour] = last_geom
        end
        push!(all_unions, unions)
    end
    println("  Done.")

    # Create 2x2 figure
    fig = Figure(size=(1800, 1200))
    axes = []
    title_obs = []

    positions = [(1, 1), (1, 2), (2, 1), (2, 2)]
    for (i, fd) in enumerate(fire_data)
        row, col = positions[i]
        lon_min, lon_max = fd.ext.X
        lat_min, lat_max = fd.ext.Y
        pad = 0.02 * max(lon_max - lon_min, lat_max - lat_min)

        t_obs = Observable{Any}(fd.fire_name)
        push!(title_obs, t_obs)

        ax = GeoAxis(fig[row, col];
            dest="+proj=longlat",
            limits=(lon_min - pad, lon_max + pad, lat_min - pad, lat_max + pad),
            title=t_obs,
            titlesize=28,
            xticklabelsvisible=false,
            yticklabelsvisible=false,
            xticksvisible=false,
            yticksvisible=false
        )
        push!(axes, (ax, fd, lon_min - pad, lon_max + pad, lat_min - pad, lat_max + pad))
    end

    # Record animation
    record(fig, output_path, 0:max_hour; framerate=framerate) do hour
        for (i, (ax, fd, lon_min, lon_max, lat_min, lat_max)) in enumerate(axes)
            empty!(ax)
            geom = all_unions[i][hour]
            if geom !== nothing
                hull = LibGEOS.convexhull(geom)
                poly!(ax, hull; color=(:yellow, 0.5), strokecolor=:goldenrod, strokewidth=1)
                poly!(ax, geom; color=(:red, 0.7), strokecolor=:darkred, strokewidth=1)

                center_lat = (lat_min + lat_max) / 2
                deg2km² = 111.0^2 * cosd(center_lat)
                fire_area_km2 = LibGEOS.area(geom) * deg2km²
                hull_area_km2 = LibGEOS.area(hull) * deg2km²

                title_obs[i][] = rich(
                    "$(fd.fire_name) — Hour $hour\n\n",
                    rich("Fire: $(round(fire_area_km2, digits=2)) km²", color=:darkred),
                    "   ",
                    rich("Hull: $(round(hull_area_km2, digits=2)) km²", color=:goldenrod)
                )
            end
        end
    end

    println("  Combined GIF saved to: $output_path")
    return output_path
end

# Process all ororatech datasets
data_dir = joinpath(@__DIR__, "..", "data", "ororatech")
output_dir = joinpath(@__DIR__, "ororatech_gifs")
mkpath(output_dir)

json_files = filter(f -> endswith(f, ".json"), readdir(data_dir))

for json_file in json_files
    input_path = joinpath(data_dir, json_file)
    output_name = replace(json_file, ".json" => ".gif") |> x -> replace(x, " " => "_")
    output_path = joinpath(output_dir, output_name)
    create_fire_gif(input_path; output_path=output_path)
    println()
end

# Create combined 2x2 GIF (excluding Lee Fire)
combined_files = [
    joinpath(data_dir, "Elk Fire 24h.json"),
    joinpath(data_dir, "Red Canyon Fire.json"),
    joinpath(data_dir, "South Rim 24h (1).json"),
    joinpath(data_dir, "Turner Gulch Fire 24h.json")
]
create_combined_gif(combined_files; output_path=joinpath(output_dir, "combined_fires.gif"))

println("\nAll GIFs created in: $output_dir")

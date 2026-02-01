#-----------------------------------------------------------------------------# init_map
function init_map(ext::Extents.Extent = extent2;
        title = "Marshall Wildfire Map Viewer",
        provider = TileProviders.OpenStreetMap(), #TileProviders.Google(:terrain),
        figure = Figure(size=(1200, 1000))
    )
    navbar = Box(figure[1, 1:2], color=:lightgray, height=40, strokevisible=false)
    Label(figure[1, 1:2][1,1], title, fontsize=18, halign=:left)
    Label(figure[1, 1:2][1,2], "(ctrl + click to reset view)", fontsize=12, halign=:right)

    sidebar = Box(figure[2, 1], color=:lightgray, width=200, strokevisible=false)
    axis = GeoAxis(figure[2, 2], dest="+proj=webmerc +datum=WGS84", panbutton=Mouse.left)
    deregister_interaction!(axis, :rectanglezoom)  # Needed to enable panning with left mouse button
    hidedecorations!(axis, label=false, ticklabels=false, ticks=false, grid=true)

    m = Tyler.Map(ext; figure, axis)
    wait(m)

    display(figure)
    (; figure, axis)  # Return figure/axis so further interactions can be added
end

function add_marshall_perimeter!(ax::GeoAxis, color=:red)
    data = get_perimeter()
    poly!(ax, data.geometry, color=(color, .2), strokecolor=color, strokewidth=1)
end

function add_buildings!(ax::GeoAxis, color=:black)
    data = get_building_footprints()
    poly!(ax, data.geometry, color=(color, 0.5), strokecolor=(color, 0.7), strokewidth=0.3)
end

function add_power_lines!(ax::GeoAxis, color=:orange)
    data = get_power_lines()
    lines!(ax, data.geometry, color=color, linewidth=1, linestyle=:dash)
end

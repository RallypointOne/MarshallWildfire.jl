using OSMGeocoder, GeoMakie

boulder = geocode(county="Boulder", state="CO")

plot(boulder.geometry, axis=(;type = GeoAxis))



using Landfire, Rasters, ArchGDAL

f13 = Landfire.products(layer="FBFM13", conus=true);

data = Landfire.Dataset(f13, boulder);

file = get(data);  # Send request to Landfire Service

plot(Raster(file))

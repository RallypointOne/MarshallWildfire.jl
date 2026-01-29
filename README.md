# MarshallWildfire.jl

[![CI](https://github.com/RallypointOne/MarshallWildfire.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/RallypointOne/MarshallWildfire.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/RallypointOne/MarshallWildfire.jl/actions/workflows/deploy-site.yml/badge.svg)](https://github.com/RallypointOne/MarshallWildfire.jl/actions/workflows/deploy-site.yml)

A Julia package for analyzing and simulating the **Marshall Fire** - the most destructive wildfire in Colorado history.

## The Marshall Fire

On December 30, 2021, a wildfire ignited in Boulder County, Colorado. Driven by extreme winds (gusts over 100 mph), the fire spread rapidly through suburban areas, destroying over 1,000 structures and burning approximately 6,000 acres in just hours.

| Attribute | Value |
|-----------|-------|
| **Start Date** | December 30, 2021 |
| **Location** | Boulder County, CO |
| **Ignition Point** | -105.231°, 39.955° |
| **Cause** | Downed power lines (wind) |
| **Structures Destroyed** | 1,084 |

## Features

- **Data Acquisition**: Automated download and processing of:
  - Fire perimeter from Boulder County ArcGIS
  - LANDFIRE fuel, slope, and aspect data
  - HRRR wind field data (hourly)
  - Building footprints and power lines from OpenStreetMap

- **Visualization**: Publication-quality maps using GLMakie, GeoMakie, and Tyler.jl:
  - Fire perimeter with basemap tiles
  - Building and power line overlays
  - LANDFIRE layer visualization
  - Animated wind field GIFs

- **Fire Spread Modeling**:
  - Level set method for front tracking
  - Rothermel fire spread rate model
  - Anderson 13 fuel model properties
  - Spatiotemporal wind interpolation via GeoSurrogates.jl

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/RallypointOne/MarshallWildfire.jl")
```

## Quick Start

```julia
using MarshallWildfire

# Download and view fire perimeter
perim = get_perimeter()

# Get LANDFIRE data (fuel, slope, aspect)
landfire = get_landfire_data()

# Get HRRR wind data for the fire duration
wind = get_hrrr_data()

# Generate visualizations
plot_marshall_perimeter()
plot_perimeter_with_buildings()
plot_fuel_flammability()
plot_hrrr_wind()  # Creates animated GIF
```

## Analysis

Analysis results are available at the [project website](https://rallypointone.github.io/MarshallWildfire.jl/).

## Data Sources

| Dataset | Source | Resolution |
|---------|--------|------------|
| Fire Perimeter | Boulder County ArcGIS | Vector |
| Fuel/Terrain | LANDFIRE (USGS) | 30m |
| Wind | HRRR (NOAA) | 3km, hourly |
| Buildings | OpenStreetMap | Vector |
| Power Lines | OpenStreetMap | Vector |

## References

- Rothermel, R. C. (1972). A mathematical model for predicting fire spread in wildland fuels. *USDA Forest Service Research Paper INT-115*.
- Osher, S., & Sethian, J. A. (1988). Fronts propagating with curvature-dependent speed. *Journal of Computational Physics*, 79(1), 12-49.

## License

MIT License

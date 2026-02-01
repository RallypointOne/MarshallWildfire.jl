#=
Rothermel Fire Spread Model

Implementation of Rothermel's (1972) mathematical model for predicting fire spread
in wildland fuels.

Reference:
    Rothermel, R. C. (1972). A mathematical model for predicting fire spread in
    wildland fuels. USDA Forest Service Research Paper INT-115.
=#

#-----------------------------------------------------------------------------# Physical Constants
const ρp = 513.0      # Ovendry particle density (kg/m³)
const Se = 0.01       # Effective mineral content (fraction)
const St = 0.0555     # Total mineral content (fraction)

#-----------------------------------------------------------------------------# Anderson 13 Fuel Models
# Standard fuel model properties from Anderson (1982)
# w0: total fuel load (kg/m²)
# σ: surface-area-to-volume ratio (1/m)
# Mx: dead fuel moisture of extinction (fraction)
# h: heat content (kJ/kg)
# δ: fuel bed depth (m)

"""
    FuelModel

Properties for a fuel model used in Rothermel's fire spread equations.

# Fields
- `name::String` - Descriptive name
- `w0::Float64` - Total fuel load (kg/m²)
- `σ::Float64` - Surface-area-to-volume ratio (1/m)
- `Mx::Float64` - Dead fuel moisture of extinction (fraction)
- `h::Float64` - Heat content (kJ/kg)
- `δ::Float64` - Fuel bed depth (m)
"""
struct FuelModel
    name::String
    w0::Float64   # Fuel load (kg/m²)
    σ::Float64    # SAV ratio (1/m)
    Mx::Float64   # Moisture of extinction (fraction)
    h::Float64    # Heat content (kJ/kg)
    δ::Float64    # Fuel bed depth (m)
end

# Anderson 13 Fuel Models
const FUEL_MODELS = Dict{Int, FuelModel}(
    # Grass Group
    1  => FuelModel("Short grass (1 ft)",           0.166, 11483.0, 0.12, 18622.0, 0.305),
    2  => FuelModel("Timber grass/understory",      0.896,  9843.0, 0.15, 18622.0, 0.305),
    3  => FuelModel("Tall grass (2.5 ft)",          1.345,  4921.0, 0.25, 18622.0, 0.762),

    # Shrub Group
    4  => FuelModel("Chaparral (6 ft)",             2.242,  6562.0, 0.20, 18622.0, 1.829),
    5  => FuelModel("Brush (2 ft)",                 0.448,  6562.0, 0.20, 18622.0, 0.610),
    6  => FuelModel("Dormant brush/hardwood slash", 0.673,  5741.0, 0.25, 18622.0, 0.762),
    7  => FuelModel("Southern rough",               0.507,  5741.0, 0.40, 18622.0, 0.762),

    # Timber Litter Group
    8  => FuelModel("Closed timber litter",         0.673,  6562.0, 0.30, 18622.0, 0.061),
    9  => FuelModel("Hardwood litter",              1.307,  8203.0, 0.25, 18622.0, 0.061),
    10 => FuelModel("Timber w/ understory",         1.345,  6562.0, 0.25, 18622.0, 0.305),

    # Slash Group
    11 => FuelModel("Light logging slash",          0.673,  4921.0, 0.15, 18622.0, 0.305),
    12 => FuelModel("Medium logging slash",         1.794,  4921.0, 0.20, 18622.0, 0.701),
    13 => FuelModel("Heavy logging slash",          3.140,  4921.0, 0.25, 18622.0, 0.914),
)

# Non-burnable fuel codes (from LANDFIRE)
const NON_BURNABLE_CODES = Set([91, 92, 93, 98, 99])

const DEFAULT_FUEL = FUEL_MODELS[1]

"""
    get_fuel_model(code::Int) -> FuelModel

Get fuel model properties for a given Anderson 13 fuel code.
Returns default (short grass) for unknown codes.
"""
function get_fuel_model(code::Int)
    if code in NON_BURNABLE_CODES
        return nothing
    end
    return get(FUEL_MODELS, code, DEFAULT_FUEL)
end

#-----------------------------------------------------------------------------# Intermediate Calculations

"""
    packing_ratio(fuel::FuelModel) -> NamedTuple

Calculate packing ratios.
- β: actual packing ratio
- βop: optimum packing ratio
- βratio: ratio of actual to optimum
"""
function packing_ratio(fuel::FuelModel)
    β = fuel.w0 / (ρp * fuel.δ)
    βop = 3.348 * fuel.σ^(-0.8189)
    βratio = β / βop
    return (; β, βop, βratio)
end

"""
    reaction_intensity(fuel::FuelModel, Mf::Real) -> IR

Calculate reaction intensity (kW/m²).

# Arguments
- `fuel`: FuelModel properties
- `Mf`: Fuel moisture content (fraction)
"""
function reaction_intensity(fuel::FuelModel, Mf::Real)
    (; β, βop, βratio) = packing_ratio(fuel)

    # Maximum reaction velocity
    Γmax = fuel.σ^1.5 / (495.0 + 0.0594 * fuel.σ^1.5)

    # Optimum reaction velocity
    A = 133.0 * fuel.σ^(-0.7913)

    # Guard against numerical issues with extreme β ratios
    # βratio^A can overflow/underflow for extreme values
    βratio_safe = clamp(βratio, 0.01, 100.0)
    exp_term = A * (1.0 - βratio_safe)
    exp_term_clamped = clamp(exp_term, -20.0, 20.0)

    Γ = Γmax * βratio_safe^A * exp(exp_term_clamped)

    # Moisture damping coefficient
    rm = Mf / fuel.Mx  # Moisture ratio
    ηM = clamp(1.0 - 2.59*rm + 5.11*rm^2 - 3.52*rm^3, 0.0, 1.0)

    # Mineral damping coefficient
    ηs = 0.174 * Se^(-0.19)

    # Reaction intensity
    IR = Γ * fuel.w0 * fuel.h * ηM * ηs

    return IR
end

"""
    propagating_flux_ratio(fuel::FuelModel) -> ξ

Calculate the propagating flux ratio.
"""
function propagating_flux_ratio(fuel::FuelModel)
    (; β) = packing_ratio(fuel)

    # Clamp the exponent to prevent overflow
    # The original formula can produce very large values for high SAV ratios
    exponent = (0.792 + 0.681 * sqrt(fuel.σ)) * (β + 0.1)
    exponent_clamped = min(exponent, 20.0)  # exp(20) ≈ 485 million

    ξ = exp(exponent_clamped) / (192.0 + 0.2595 * fuel.σ)
    return ξ
end

"""
    heat_sink(fuel::FuelModel, Mf::Real) -> Q

Calculate the heat sink term (effective heating number × heat of preignition × bulk density).
"""
function heat_sink(fuel::FuelModel, Mf::Real)
    ε = exp(-138.0 / fuel.σ)           # Effective heating number
    Qig = 250.0 + 1116.0 * Mf          # Heat of preignition (kJ/kg)
    ρb = fuel.w0 / fuel.δ              # Bulk density (kg/m³)
    return ρb * ε * Qig
end

"""
    no_wind_no_slope_rate(fuel::FuelModel, Mf::Real) -> R0

Calculate the no-wind, no-slope spread rate (m/s).
"""
function no_wind_no_slope_rate(fuel::FuelModel, Mf::Real)
    IR = reaction_intensity(fuel, Mf)
    ξ = propagating_flux_ratio(fuel)
    Q = heat_sink(fuel, Mf)

    # Guard against division by zero
    if Q <= 0.0
        return 0.0
    end

    # R0 in m/min, convert to m/s
    # Clamp to reasonable maximum (100 m/s ≈ 360 km/h)
    R0 = min((IR * ξ) / Q / 60.0, 100.0)

    # Return 0 if result is NaN or negative
    return (isfinite(R0) && R0 > 0) ? R0 : 0.0
end

"""
    wind_factor(fuel::FuelModel, wind_speed::Real) -> ϕw

Calculate the wind factor for spread rate.

# Arguments
- `fuel`: FuelModel properties
- `wind_speed`: Wind speed at midflame height (m/s)
"""
function wind_factor(fuel::FuelModel, wind_speed::Real)
    (; βratio) = packing_ratio(fuel)

    C = 7.47 * exp(-0.133 * fuel.σ^0.55)
    B = 0.02526 * fuel.σ^0.54
    E = 0.715 * exp(-3.59e-4 * fuel.σ)

    # Convert wind speed to m/min for original formula
    # Clamp to reasonable range (0 to 50 m/s = 112 mph)
    U = clamp(wind_speed, 0.0, 50.0) * 60.0

    # Wind factor (dimensionless multiplier)
    # Clamp to realistic maximum (~50 is typical upper bound for extreme winds)
    ϕw = min(C * U^B * βratio^(-E), 1000.0)

    return ϕw
end

"""
    slope_factor(fuel::FuelModel, slope::Real) -> ϕs

Calculate the slope factor for spread rate.

# Arguments
- `fuel`: FuelModel properties
- `slope`: Terrain slope (radians)
"""
function slope_factor(fuel::FuelModel, slope::Real)
    (; β) = packing_ratio(fuel)

    # Clamp slope to prevent tan() overflow (max ~80 degrees)
    slope_clamped = clamp(abs(slope), 0.0, 1.4)

    # Slope factor (dimensionless multiplier)
    # Clamp to realistic maximum (~20 is typical upper bound for steep slopes)
    ϕs = min(5.275 * β^(-0.3) * tan(slope_clamped)^2, 1000.0)

    return ϕs
end

#-----------------------------------------------------------------------------# Main Spread Rate Function

"""
    rothermel_spread_rate(fuel, wind_speed, wind_dir, slope, aspect; Mf=0.05)

Calculate fire spread rate and direction using Rothermel's (1972) model.

# Arguments
- `fuel::FuelModel`: Fuel model properties
- `wind_speed::Real`: Wind speed at midflame height (m/s)
- `wind_dir::Real`: Wind direction (radians, direction wind is coming FROM, meteorological convention)
- `slope::Real`: Terrain slope (radians)
- `aspect::Real`: Terrain aspect (radians, downslope direction)
- `Mf::Real=0.05`: Fuel moisture content (fraction, default 5%)

# Returns
Named tuple with:
- `rate`: Fire spread rate (m/s)
- `direction`: Direction of maximum spread (radians)
- `R0`: No-wind, no-slope spread rate (m/s)

# Notes
- Wind direction follows meteorological convention (direction FROM which wind blows)
- Aspect is the downslope direction (direction water would flow)
- Returns zero spread rate for non-burnable fuels or moisture above extinction

# Example
```julia
fuel = FUEL_MODELS[1]  # Short grass
result = rothermel_spread_rate(fuel, 10.0, 0.0, deg2rad(10), 0.0; Mf=0.05)
println("Spread rate: \$(result.rate) m/s")
```
"""
function rothermel_spread_rate(fuel::FuelModel, wind_speed::Real, wind_dir::Real,
                                slope::Real, aspect::Real; Mf::Real=0.05)
    # Check for extinction
    if Mf >= fuel.Mx
        return (rate = 0.0, direction = 0.0, R0 = 0.0)
    end

    # Base spread rate (no wind, no slope)
    R0 = no_wind_no_slope_rate(fuel, Mf)

    # Wind and slope factors
    ϕw = wind_factor(fuel, wind_speed)
    ϕs = slope_factor(fuel, slope)

    # Direction vectors
    # Wind pushes fire in direction it's blowing TO (opposite of wind_dir)
    wind_push_dir = wind_dir + π
    # Slope pushes fire uphill (opposite of aspect which points downhill)
    slope_push_dir = aspect + π

    # Vector combination of base rate + wind effect + slope effect
    Rw = R0 * ϕw
    Rs = R0 * ϕs

    Rx = R0 + Rw * cos(wind_push_dir) + Rs * cos(slope_push_dir)
    Ry = Rw * sin(wind_push_dir) + Rs * sin(slope_push_dir)

    # Resultant spread rate and direction
    R = sqrt(Rx^2 + Ry^2)
    θ = atan(Ry, Rx)

    return (rate = R, direction = θ, R0 = R0)
end

# Convenience method using fuel code instead of FuelModel
function rothermel_spread_rate(fuel_code::Int, wind_speed::Real, wind_dir::Real,
                                slope::Real, aspect::Real; Mf::Real=0.05)
    fuel = get_fuel_model(fuel_code)
    if isnothing(fuel)
        return (rate = 0.0, direction = 0.0, R0 = 0.0)
    end
    return rothermel_spread_rate(fuel, wind_speed, wind_dir, slope, aspect; Mf)
end

#-----------------------------------------------------------------------------# Diagnostic Functions

"""
    diagnose_spread(fuel::FuelModel, wind_speed::Real, slope::Real; Mf=0.05)

Print diagnostic information about fire spread components.
Useful for understanding which factors dominate spread.
"""
function diagnose_spread(fuel::FuelModel, wind_speed::Real, slope::Real; Mf::Real=0.05)
    R0 = no_wind_no_slope_rate(fuel, Mf)
    ϕw = wind_factor(fuel, wind_speed)
    ϕs = slope_factor(fuel, slope)

    IR = reaction_intensity(fuel, Mf)
    (; β, βop, βratio) = packing_ratio(fuel)

    println("Rothermel Spread Rate Diagnostics")
    println("=" ^ 40)
    println("Fuel Model: $(fuel.name)")
    println()
    println("Fuel Properties:")
    println("  Load (w0):     $(fuel.w0) kg/m²")
    println("  SAV ratio (σ): $(fuel.σ) 1/m")
    println("  Depth (δ):     $(fuel.δ) m")
    println("  Mx:            $(fuel.Mx * 100)%")
    println()
    println("Packing Ratios:")
    println("  β:       $(round(β, digits=4))")
    println("  βop:     $(round(βop, digits=4))")
    println("  β/βop:   $(round(βratio, digits=4))")
    println()
    println("Environmental Inputs:")
    println("  Fuel moisture: $(Mf * 100)%")
    println("  Wind speed:    $(wind_speed) m/s ($(round(wind_speed * 2.237, digits=1)) mph)")
    println("  Slope:         $(round(rad2deg(slope), digits=1))°")
    println()
    println("Spread Components:")
    println("  Reaction intensity (IR): $(round(IR, digits=1)) kW/m²")
    println("  Base rate (R0):          $(round(R0 * 60, digits=3)) m/min")
    println("  Wind factor (ϕw):        $(round(ϕw, digits=2))")
    println("  Slope factor (ϕs):       $(round(ϕs, digits=2))")
    println()
    println("Spread Rates:")
    println("  No-wind/no-slope: $(round(R0 * 60, digits=3)) m/min")
    println("  With wind only:   $(round(R0 * (1 + ϕw) * 60, digits=3)) m/min")
    println("  With slope only:  $(round(R0 * (1 + ϕs) * 60, digits=3)) m/min")
    println("  Combined max:     $(round(R0 * (1 + ϕw + ϕs) * 60, digits=3)) m/min")
end

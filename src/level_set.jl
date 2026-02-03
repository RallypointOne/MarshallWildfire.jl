# Level set model based on
# [1] J. Mandel, J. D. Beezley, and A. K. Kochanski, “Coupled atmosphere-wildland fire modeling with WRF 3.3 and SFIRE 2011,” Geosci. Model Dev., vol. 4, no. 3, pp. 591–610, Jul. 2011, doi: 10.5194/gmd-4-591-2011.
# [2] R. C. Rothermel, Ogden, UT: Intermountain Forest and Range Experiment Station, Forest Service, United States Department of Agriculture, 1972.

function fire_spread_rate(x, y, Ψ_grad_x, Ψ_grad_y, t, wind_fn, terrain_slope_fn, fuel_properties_fn)
    # Unpack fuel properties from named tuple
    fuel = fuel_properties_fn(x, y)
    a = fuel.a
    w = fuel.w
    w_l = fuel.w_l
    δ_m = fuel.δ_m
    σ = fuel.σ
    M_x = fuel.M_x
    ρ_P = fuel.ρ_P
    S_T = fuel.S_T
    S_E = fuel.S_E
    h = fuel.h
    M_f = fuel.M_f

    # Get terrain slope at this location
    terrain_grad_x, terrain_grad_y = terrain_slope_fn(x, y)  # Returns (∂z/∂x, ∂z/∂y)

    # Get wind at this location and time
    wind_x, wind_y = wind_fn(x, y, t)

    # Compute normal vector to fire front
    gn = sqrt(Ψ_grad_x^2 + Ψ_grad_y^2)
    n_x = Ψ_grad_x / gn
    n_y = Ψ_grad_y / gn

    # Wind-normal interaction (wind correction factor)
    wind_dot_normal = wind_x * n_x + wind_y * n_y
    U = abs(wind_dot_normal)

    # Slope in direction of fire spread
    tanϕ = terrain_grad_x * n_x + terrain_grad_y * n_y

    ## ROTHERMEL FIRE SPREAD RATE EQUATIONS

    βop = 3.348 * σ^(-0.8189)        # from Rothermel, eq (37)
    w0 = w_l / (1 + M_f)
    ρb = w0 / δ_m                    # bulk density
    β = ρb / ρ_P                     # packing ratio
    ξ = exp((0.792 + 0.618 * σ^0.5) * (β + 0.1)) / (192 + 0.25965 * σ)
    ηs = 0.174 * S_E^(-0.19)         # mineral damping coefficient
    ηM = 1 - 2.59 * M_f / M_x + 5.11 * (M_f / M_x)^2 - 3.52 * (M_f / M_x)^3  # moisture damping
    wn = w0 / (1 + S_T)
    Γmax = σ^1.5 / (495 + 0.594 * σ^1.5)
    A = 1 / (4.77 * σ^0.1 - 7.27)
    Γ = Γmax * (β / βop)^A * exp(A * (1 - β / βop))  # optimum reaction velocity
    ϵ = exp(-138 / σ)                # effective heating number
    Qig = 250 * β + 1116 * M_f       # heat of preignition
    C = 7.47 * exp(-0.133 * σ^0.55)  # eq (48)
    B = 0.02526 * σ^0.54              # eq (49) - wind exponent
    Ua = a * U                        # adjusted wind speed
    E = 0.715 * exp(-0.000359 * σ)   # eq (50)
    IR = Γ * wn * h * ηM * ηs        # reaction intensity
    R0 = IR * ξ / (ρb * ϵ * Qig)     # spread rate without wind/slope
    ϕw = C * (Ua)^B * (β / βop)^(-E)  # eq (47) - wind factor
    ϕS = 5.275 * β^(-0.3) * tanϕ^2              # slope factor

    fuel_scale = 1.0  # scaling factor (reduced from 10000 to prevent gradient explosion)
    S = fuel_scale * R0 * (1 + ϕw + ϕS)  # final fire spread rate

    return S
end

# Register fire_spread_rate at module level
@register_symbolic fire_spread_rate(x, y, Ψ_grad_x, Ψ_grad_y, t, wind_fn, terrain_slope_fn, fuel_properties_fn)

"""
    level_set_equation(wind_fn, terrain_slope_fn, fuel_properties_fn)

Creates the Rothermel level set PDE equation: ∂Ψ/∂t + S * |∇Ψ| = 0

# Arguments
- `wind_fn`: Function with signature (x, y, t) -> (wind_x, wind_y)
- `terrain_slope_fn`: Function with signature (x, y) -> (slope_x, slope_y)
- `fuel_properties_fn`: Function with signature (x, y) -> NamedTuple

# Returns
A named tuple with:
- `equation`: The PDE equation
- `x`, `y`: Spatial parameters
- `Ψ`: The level set function variable
- `Dt`, `Dx`, `Dy`: Differential operators

# Example
```julia
# Define custom functions
my_wind(x, y, t) = (50.0, 50.0)
my_slope(x, y) = (0.0, 0.0)
my_fuel(x, y) = get_fuel_properties(x, y)

# Create the PDE equation
pde = level_set_equation(my_wind, my_slope, my_fuel)
```
"""
function level_set_equation(wind_fn, terrain_slope_fn, fuel_properties_fn)
    # Declare symbolic variables
    @parameters x y
    @variables Ψ(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)

    # Fire spread rate at each point - pass the functions directly
    S = fire_spread_rate(x, y, Dx(Ψ(t, x, y)), Dy(Ψ(t, x, y)), t, wind_fn, terrain_slope_fn, fuel_properties_fn)

    # Gradient magnitude
    gn = sqrt(Dx(Ψ(t, x, y))^2 + Dy(Ψ(t, x, y))^2)

    # Level set equation: ∂Ψ/∂t + S * |∇Ψ| = 0
    equation = Dt(Ψ(t, x, y)) + S * gn ~ 0

    return (
        equation = equation,
        t = t,
        x = x,
        y = y,
        Ψ = Ψ,
        Dt = Dt,
        Dx = Dx,
        Dy = Dy
    )
end


#-----------------------------------------------------------------------------# LevelSet
# Level set method for tracking propagating fronts
# Reference: Osher & Sethian (1988), Sethian (1999)

"""
    LevelSet{T<:AbstractFloat}

Represents a propagating front as the zero level set of a scalar function ψ.

Convention:
- ψ < 0: Inside the front (burned area for fire)
- ψ = 0: The front itself (fire perimeter)
- ψ > 0: Outside the front (unburned area)

# Fields
- `ψ::Matrix{T}` - Level set function values on grid
- `xs::Vector{T}` - X coordinates of grid points
- `ys::Vector{T}` - Y coordinates of grid points
- `dx::T` - Grid spacing in x direction
- `dy::T` - Grid spacing in y direction
"""
mutable struct LevelSet{T<:AbstractFloat}
    ψ::Matrix{T}
    xs::Vector{T}
    ys::Vector{T}
    dx::T
    dy::T
end

"""
    LevelSet(extent::Extent; nx=100, ny=100, T=Float64)

Create a LevelSet on a regular grid covering the given extent.
Initializes ψ to a constant positive value (no burned area).
"""
function LevelSet(extent::Extent; nx::Int=100, ny::Int=100, T::Type{<:AbstractFloat}=Float64)
    xs = collect(T, range(extent.X[1], extent.X[2], length=nx))
    ys = collect(T, range(extent.Y[1], extent.Y[2], length=ny))
    dx = T(xs[2] - xs[1])
    dy = T(ys[2] - ys[1])
    ψ = ones(T, nx, ny)
    LevelSet{T}(ψ, xs, ys, dx, dy)
end

"""
    LevelSet(xs::AbstractVector, ys::AbstractVector)

Create a LevelSet from explicit coordinate vectors.
"""
function LevelSet(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:AbstractFloat
    nx, ny = length(xs), length(ys)
    dx = T(xs[2] - xs[1])
    dy = T(ys[2] - ys[1])
    ψ = ones(T, nx, ny)
    LevelSet{T}(ψ, collect(xs), collect(ys), dx, dy)
end

Base.size(ls::LevelSet) = size(ls.ψ)
Base.eltype(::LevelSet{T}) where T = T

#-----------------------------------------------------------------------------# Initialization

"""
    init_circle!(ls::LevelSet, x0, y0, r)

Initialize the level set as a signed distance function from a circle.
Points inside the circle (distance < r) have ψ < 0.
"""
function init_circle!(ls::LevelSet{T}, x0::Real, y0::Real, r::Real) where T
    for (i, x) in enumerate(ls.xs), (j, y) in enumerate(ls.ys)
        ls.ψ[i, j] = sqrt((x - x0)^2 + (y - y0)^2) - r
    end
    return ls
end

"""
    init_point!(ls::LevelSet, x0, y0; burn_radius=0.001)

Initialize the level set from a point ignition.
Creates a small burned circle of given radius around the point.
"""
function init_point!(ls::LevelSet{T}, x0::Real, y0::Real; burn_radius::Real=0.001) where T
    init_circle!(ls, x0, y0, burn_radius)
end

"""
    init_signed_distance!(ls::LevelSet, x0, y0)

Initialize as signed distance from a point (positive everywhere except at the point).
Useful for tracking expansion from a point source.
"""
function init_signed_distance!(ls::LevelSet{T}, x0::Real, y0::Real) where T
    for (i, x) in enumerate(ls.xs), (j, y) in enumerate(ls.ys)
        ls.ψ[i, j] = sqrt((x - x0)^2 + (y - y0)^2)
    end
    return ls
end

"""
    init_ellipse!(ls::LevelSet, x0, y0, a, b; θ=0.0)

Initialize the level set as a signed distance function from an ellipse.
Points inside the ellipse have ψ < 0.

# Arguments
- `x0`, `y0`: Center of the ellipse
- `a`: Semi-major axis length
- `b`: Semi-minor axis length
- `θ`: Rotation angle in radians (default: 0.0)
"""
function init_ellipse!(ls::LevelSet{T}, x0::Real, y0::Real, a::Real, b::Real; θ::Real=0.0) where T
    cosθ, sinθ = cos(θ), sin(θ)
    for (i, x) in enumerate(ls.xs), (j, y) in enumerate(ls.ys)
        # Translate to ellipse center
        dx = x - x0
        dy = y - y0
        # Rotate to ellipse coordinates
        xr = cosθ * dx + sinθ * dy
        yr = -sinθ * dx + cosθ * dy
        # Approximate signed distance (exact on boundary, approximate elsewhere)
        ls.ψ[i, j] = sqrt((xr/a)^2 + (yr/b)^2) - 1.0
    end
    return ls
end

#-----------------------------------------------------------------------------# Gradient Computation

"""
    gradient_upwind(ψ, i, j, dx, dy)

Compute upwind gradient magnitude using Godunov's method.
This ensures numerical stability for the level set equation.
"""
function gradient_upwind(ψ::Matrix, i::Int, j::Int, dx::Real, dy::Real)
    nx, ny = size(ψ)

    # Boundary check - use one-sided differences at boundaries
    if i == 1 || i == nx || j == 1 || j == ny
        return zero(eltype(ψ))
    end

    # Check for NaN in neighboring cells
    if !isfinite(ψ[i-1, j]) || !isfinite(ψ[i+1, j]) ||
       !isfinite(ψ[i, j-1]) || !isfinite(ψ[i, j+1])
        return zero(eltype(ψ))
    end

    # Backward and forward differences
    Dxm = (ψ[i, j] - ψ[i-1, j]) / dx  # D⁻ₓ
    Dxp = (ψ[i+1, j] - ψ[i, j]) / dx  # D⁺ₓ
    Dym = (ψ[i, j] - ψ[i, j-1]) / dy  # D⁻ᵧ
    Dyp = (ψ[i, j+1] - ψ[i, j]) / dy  # D⁺ᵧ

    # Godunov upwind scheme for expanding front (S > 0)
    # max(D⁻, 0)² + min(D⁺, 0)² for each direction
    Dx = max(max(Dxm, 0.0), -min(Dxp, 0.0))
    Dy = max(max(Dym, 0.0), -min(Dyp, 0.0))

    return sqrt(Dx^2 + Dy^2)
end

"""
    gradient_central(ψ, i, j, dx, dy)

Compute gradient using central differences.
Less stable but useful for visualization.
"""
function gradient_central(ψ::Matrix, i::Int, j::Int, dx::Real, dy::Real)
    nx, ny = size(ψ)

    if i == 1 || i == nx || j == 1 || j == ny
        return (zero(eltype(ψ)), zero(eltype(ψ)))
    end

    dψdx = (ψ[i+1, j] - ψ[i-1, j]) / (2dx)
    dψdy = (ψ[i, j+1] - ψ[i, j-1]) / (2dy)

    return (dψdx, dψdy)
end

#-----------------------------------------------------------------------------# Evolution

"""
    step!(ls::LevelSet, speed::Function, t::Real, dt::Real)

Advance the level set by one time step.

Solves: ∂ψ/∂t + S(x,y,t)|∇ψ| = 0

# Arguments
- `ls`: Level set to update (modified in place)
- `speed`: Function `speed(x, y, t) -> Real` returning local front speed
- `t`: Current time
- `dt`: Time step size
"""
function step!(ls::LevelSet{T}, speed::Function, t::Real, dt::Real) where T
    ψ = ls.ψ
    nx, ny = size(ψ)
    dx, dy = ls.dx, ls.dy

    ψ_new = copy(ψ)

    for i in 2:nx-1, j in 2:ny-1
        x, y = ls.xs[i], ls.ys[j]

        # Get local speed
        S = speed(x, y, t)

        # Skip if speed is zero, negative, NaN, or Inf (no spread or invalid)
        if !(S > 0) || !isfinite(S)
            continue
        end

        # Skip if current ψ value is already NaN (prevent propagation)
        if !isfinite(ψ[i, j])
            continue
        end

        # Compute upwind gradient magnitude
        grad_mag = gradient_upwind(ψ, i, j, dx, dy)

        # Skip if gradient is invalid
        if !isfinite(grad_mag)
            continue
        end

        # Level set equation: ψ_t + S|∇ψ| = 0
        # Forward Euler: ψⁿ⁺¹ = ψⁿ - dt * S * |∇ψ|
        ψ_new[i, j] = ψ[i, j] - dt * S * grad_mag
    end

    ls.ψ .= ψ_new
    return ls
end

"""
    step!(ls::LevelSet, speed::Matrix, dt::Real)

Advance level set using a precomputed speed field (matrix of speeds at each grid point).
"""
function step!(ls::LevelSet{T}, speed::Matrix, dt::Real) where T
    ψ = ls.ψ
    nx, ny = size(ψ)
    dx, dy = ls.dx, ls.dy

    @assert size(speed) == size(ψ) "Speed matrix must match level set dimensions"

    ψ_new = copy(ψ)

    for i in 2:nx-1, j in 2:ny-1
        S = speed[i, j]

        # Skip if speed is zero, negative, NaN, or Inf
        if !(S > 0) || !isfinite(S)
            continue
        end

        # Skip if current ψ value is already NaN
        if !isfinite(ψ[i, j])
            continue
        end

        grad_mag = gradient_upwind(ψ, i, j, dx, dy)

        # Skip if gradient is invalid
        if !isfinite(grad_mag)
            continue
        end

        ψ_new[i, j] = ψ[i, j] - dt * S * grad_mag
    end

    ls.ψ .= ψ_new
    return ls
end

"""
    step_directional!(ls::LevelSet, speed::Function, t::Real, dt::Real)

Advance the level set by one time step with direction-dependent speed.

Solves: ∂ψ/∂t + S(x,y,t,nx,ny)|∇ψ| = 0

# Arguments
- `ls`: Level set to update (modified in place)
- `speed`: Function `speed(x, y, t, nx, ny) -> Real` returning local front speed
          where (nx, ny) is the outward normal direction (unit vector)
- `t`: Current time
- `dt`: Time step size

This method allows the speed to vary based on the direction of propagation,
which is essential for modeling wind-driven fire spread.
"""
function step_directional!(ls::LevelSet{T}, speed::Function, t::Real, dt::Real) where T
    ψ = ls.ψ
    nx_grid, ny_grid = size(ψ)
    dx, dy = ls.dx, ls.dy

    ψ_new = copy(ψ)

    for i in 2:nx_grid-1, j in 2:ny_grid-1
        x, y = ls.xs[i], ls.ys[j]

        # Skip if current ψ value is already NaN (prevent propagation)
        if !isfinite(ψ[i, j])
            continue
        end

        # Compute gradient using central differences for normal direction
        # Note: gradient_central returns (d/di, d/dj) where i is first matrix index
        # In our LevelSet, ψ[i,j] is at (xs[i], ys[j]), so first index = x, second = y
        # But for heatmap plotting with transpose, the visual x corresponds to j and y to i
        # So we swap: dψ/dx uses j-differences, dψ/dy uses i-differences
        grad_i, grad_j = gradient_central(ψ, i, j, dx, dy)
        grad_mag_central = sqrt(grad_i^2 + grad_j^2)

        # Skip if gradient is too small (can't determine direction)
        if grad_mag_central < 1e-10
            continue
        end

        # Outward normal direction (points from burned to unburned)
        # Swap components to match visual coordinate system
        nx = grad_j / grad_mag_central  # x-component from j-derivative
        ny = grad_i / grad_mag_central  # y-component from i-derivative

        # Get direction-dependent speed
        S = speed(x, y, t, nx, ny)

        # Skip if speed is zero, negative, NaN, or Inf (no spread or invalid)
        if !(S > 0) || !isfinite(S)
            continue
        end

        # Compute upwind gradient magnitude for stability
        grad_mag = gradient_upwind(ψ, i, j, dx, dy)

        # Skip if gradient is invalid
        if !isfinite(grad_mag)
            continue
        end

        # Level set equation: ψ_t + S|∇ψ| = 0
        # Forward Euler: ψⁿ⁺¹ = ψⁿ - dt * S * |∇ψ|
        ψ_new[i, j] = ψ[i, j] - dt * S * grad_mag
    end

    ls.ψ .= ψ_new
    return ls
end

#-----------------------------------------------------------------------------# CFL Condition

"""
    cfl_dt(ls::LevelSet, max_speed::Real; cfl=0.5)

Compute stable time step satisfying the CFL condition.

dt ≤ CFL * min(dx, dy) / max_speed
"""
function cfl_dt(ls::LevelSet, max_speed::Real; cfl::Real=0.5)
    return cfl * min(ls.dx, ls.dy) / max_speed
end

#-----------------------------------------------------------------------------# Reinitialization

"""
    reinitialize!(ls::LevelSet; iterations=10, dt=nothing)

Reinitialize the level set to a signed distance function.

Over time, the level set can develop steep gradients or flat regions
that degrade numerical accuracy. Reinitialization restores |∇ψ| ≈ 1.

Solves: ∂ψ/∂τ + sign(ψ₀)(|∇ψ| - 1) = 0
"""
function reinitialize!(ls::LevelSet{T}; iterations::Int=10, dt::Union{Nothing,Real}=nothing) where T
    ψ = ls.ψ
    nx, ny = size(ψ)
    dx, dy = ls.dx, ls.dy

    # Store original sign
    ψ0 = copy(ψ)
    sign_ψ0 = sign.(ψ0)

    # Default dt for reinitialization
    if isnothing(dt)
        dt = 0.5 * min(dx, dy)
    end

    for _ in 1:iterations
        ψ_new = copy(ψ)

        for i in 2:nx-1, j in 2:ny-1
            s = sign_ψ0[i, j]

            # One-sided differences
            Dxm = (ψ[i, j] - ψ[i-1, j]) / dx
            Dxp = (ψ[i+1, j] - ψ[i, j]) / dx
            Dym = (ψ[i, j] - ψ[i, j-1]) / dy
            Dyp = (ψ[i, j+1] - ψ[i, j]) / dy

            # Godunov scheme based on sign
            if s > 0
                Dx = max(max(Dxm, 0.0), -min(Dxp, 0.0))
                Dy = max(max(Dym, 0.0), -min(Dyp, 0.0))
            else
                Dx = max(max(-Dxm, 0.0), min(Dxp, 0.0))
                Dy = max(max(-Dym, 0.0), min(Dyp, 0.0))
            end

            grad_mag = sqrt(Dx^2 + Dy^2)

            # Reinitialization equation
            ψ_new[i, j] = ψ[i, j] - dt * s * (grad_mag - 1.0)
        end

        ψ .= ψ_new
    end

    return ls
end

#-----------------------------------------------------------------------------# Analysis Functions

"""
    burned_area(ls::LevelSet)

Calculate the total burned area (where ψ < 0).
Returns area in coordinate units squared.
"""
function burned_area(ls::LevelSet)
    cell_area = ls.dx * ls.dy
    return count(x -> x < 0, ls.ψ) * cell_area
end

"""
    burned_mask(ls::LevelSet)

Return a boolean matrix indicating burned cells (ψ < 0).
"""
burned_mask(ls::LevelSet) = ls.ψ .< 0

"""
    front_points(ls::LevelSet; threshold=nothing)

Find approximate grid points near the fire front (where |ψ| < threshold).
Returns vector of (x, y) tuples.
"""
function front_points(ls::LevelSet{T}; threshold::Union{Nothing,Real}=nothing) where T
    if isnothing(threshold)
        threshold = max(ls.dx, ls.dy)
    end

    points = Tuple{T, T}[]
    for (i, x) in enumerate(ls.xs), (j, y) in enumerate(ls.ys)
        if abs(ls.ψ[i, j]) < threshold
            push!(points, (x, y))
        end
    end
    return points
end

"""
    perimeter_length(ls::LevelSet)

Estimate the perimeter length of the burned area.
Uses the gradient magnitude along the zero contour.
"""
function perimeter_length(ls::LevelSet)
    ψ = ls.ψ
    nx, ny = size(ψ)
    dx, dy = ls.dx, ls.dy

    length_sum = 0.0

    for i in 2:nx-1, j in 2:ny-1
        # Check if this cell crosses zero
        if ψ[i, j] * ψ[i+1, j] < 0 || ψ[i, j] * ψ[i, j+1] < 0
            # Estimate gradient magnitude
            dψdx, dψdy = gradient_central(ψ, i, j, dx, dy)
            grad_mag = sqrt(dψdx^2 + dψdy^2)

            # Add contribution (cell face length / gradient)
            if grad_mag > 1e-10
                length_sum += sqrt(dx^2 + dy^2) / grad_mag
            end
        end
    end

    return length_sum
end

#-----------------------------------------------------------------------------# Simulation

"""
    simulate(ls::LevelSet, speed::Function, duration::Real;
             dt=nothing, save_interval=nothing, reinit_interval=nothing, max_speed=1.0)

Run a level set simulation.

# Arguments
- `ls`: Initial level set (will be modified)
- `speed`: Speed function `speed(x, y, t) -> Real`
- `duration`: Total simulation time

# Keyword Arguments
- `dt`: Time step (default: computed from CFL with estimated max speed)
- `save_interval`: Time interval between saved snapshots (default: duration/10)
- `reinit_interval`: Time interval for reinitialization (default: no reinitialization)
- `max_speed`: Maximum expected speed for CFL calculation (default: 1.0)

# Returns
Vector of `(time, LevelSet)` tuples representing simulation snapshots.
"""
function simulate(ls::LevelSet, speed::Function, duration::Real;
                  dt::Union{Nothing,Real}=nothing,
                  save_interval::Union{Nothing,Real}=nothing,
                  reinit_interval::Union{Nothing,Real}=nothing,
                  max_speed::Real=1.0,
                  show_progress::Bool=true)

    # Defaults
    if isnothing(dt)
        dt = cfl_dt(ls, max_speed)
    end
    if isnothing(save_interval)
        save_interval = duration / 10
    end

    snapshots = [(0.0, deepcopy(ls))]
    t = 0.0
    last_save = 0.0
    last_reinit = 0.0

    n_steps = ceil(Int, duration / dt)
    prog = Progress(n_steps; desc="Simulating fire spread: ", enabled=show_progress)

    while t < duration
        # Take a step
        step!(ls, speed, t, dt)
        t += dt

        # Reinitialize if needed
        if !isnothing(reinit_interval) && t - last_reinit >= reinit_interval
            reinitialize!(ls)
            last_reinit = t
        end

        # Save snapshot if needed
        if t - last_save >= save_interval
            push!(snapshots, (t, deepcopy(ls)))
            last_save = t
        end

        next!(prog)
    end

    finish!(prog)

    # Always save final state
    if last_save < t
        push!(snapshots, (t, deepcopy(ls)))
    end

    return snapshots
end

"""
    simulate_directional(ls::LevelSet, speed::Function, duration::Real; kwargs...)

Run a level set simulation with direction-dependent speed.

# Arguments
- `ls`: Initial level set (will be modified)
- `speed`: Speed function `speed(x, y, t, nx, ny) -> Real` where (nx, ny) is the
          outward normal direction (unit vector pointing from burned to unburned)
- `duration`: Total simulation time

# Keyword Arguments
- `dt`: Time step (default: computed from CFL with estimated max speed)
- `save_interval`: Time interval between saved snapshots (default: duration/10)
- `reinit_interval`: Time interval for reinitialization (default: no reinitialization)
- `max_speed`: Maximum expected speed for CFL calculation (default: 1.0)

# Returns
Vector of `(time, LevelSet)` tuples representing simulation snapshots.

# Example
```julia
# Wind blowing east at 5 m/s
wind_dir = 0.0  # radians
R0 = 0.1  # base spread rate
wind_factor = 2.0

function directional_speed(x, y, t, nx, ny)
    # Speed increases when spreading in wind direction
    wind_component = cos(wind_dir) * nx + sin(wind_dir) * ny
    return R0 * (1 + wind_factor * max(0, wind_component))
end

snapshots = simulate_directional(ls, directional_speed, 60.0; max_speed=R0*(1+wind_factor))
```
"""
function simulate_directional(ls::LevelSet, speed::Function, duration::Real;
                              dt::Union{Nothing,Real}=nothing,
                              save_interval::Union{Nothing,Real}=nothing,
                              reinit_interval::Union{Nothing,Real}=nothing,
                              max_speed::Real=1.0,
                              show_progress::Bool=true)

    # Defaults
    if isnothing(dt)
        dt = cfl_dt(ls, max_speed)
    end
    if isnothing(save_interval)
        save_interval = duration / 10
    end

    snapshots = [(0.0, deepcopy(ls))]
    t = 0.0
    last_save = 0.0
    last_reinit = 0.0

    n_steps = ceil(Int, duration / dt)
    prog = Progress(n_steps; desc="Simulating fire spread: ", enabled=show_progress)

    while t < duration
        # Take a step using directional speed
        step_directional!(ls, speed, t, dt)
        t += dt

        # Reinitialize if needed
        if !isnothing(reinit_interval) && t - last_reinit >= reinit_interval
            reinitialize!(ls)
            last_reinit = t
        end

        # Save snapshot if needed
        if t - last_save >= save_interval
            push!(snapshots, (t, deepcopy(ls)))
            last_save = t
        end

        next!(prog)
    end

    finish!(prog)

    # Always save final state
    if last_save < t
        push!(snapshots, (t, deepcopy(ls)))
    end

    return snapshots
end

#-----------------------------------------------------------------------------# Makie Recipe

@recipe(LevelSetPlot, levelset) do scene
    Attributes(
        colormap = :RdYlBu,
        frontcolor = :black,
        frontwidth = 2.0,
        showfield = true,
        colorrange = automatic,
    )
end

function Makie.plot!(p::LevelSetPlot)
    ls = p[:levelset][]

    # Plot the level set field as a heatmap
    if p[:showfield][]
        crange = p[:colorrange][]
        if crange === automatic
            maxabs = maximum(abs, ls.ψ)
            crange = (-maxabs, maxabs)
        end
        heatmap!(p, ls.xs, ls.ys, ls.ψ'; colormap=p[:colormap], colorrange=crange)
    end

    # Plot the zero contour (fire front)
    contour!(p, ls.xs, ls.ys, ls.ψ'; levels=[0.0], color=p[:frontcolor], linewidth=p[:frontwidth])

    return p
end

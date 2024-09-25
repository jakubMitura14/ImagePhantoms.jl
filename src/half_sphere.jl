# src/half_sphere.jl
#=
half_sphere.jl
=#
export HalfSphere, half_sphere
include("cuboid.jl")
"""
    HalfSphere <: AbstractShape{3}
"""
struct HalfSphere <: AbstractShape{3}
    axiss::NTuple{3, Real}
end
# constructor
"""
    half_sphere(cx, cy, cz, r, Φ=0, Θ=0, value::Number = 1)
    half_sphere(center::NTuple{3,RealU}, radius::RealU, angle::NTuple{3,RealU}, v)
Construct `Object{HalfSphere}` from parameters.
"""
half_sphere(axiss,args... ; kwargs...) = Object(HalfSphere(axiss), args...; kwargs...)

# methods
volume1(::HalfSphere) = 2/3 * π # volume of unit half-sphere

ℓmax1(::HalfSphere) = 2 # maximum chord through a unit half-sphere

"""
evaluates weather a point is inside the unit half-sphere
first part is identical as in unit sphere second part is added to consider only the upper half-sphere
"""
phantom1(ob::Object3d{HalfSphere}, xyz::NTuple{3,Real}) = (sum(abs2, xyz) ≤ 1) && (xyz[3] >= 0)


# x-ray transform (line integral) of unit half-sphere
# `u,v` should be unitless
function xray1(
    hs::HalfSphere,
    u::Ru,
    v::Rv,
    ϕ::RealU, # irrelevant
    θ::RealU, # irrelevant
) where {Ru <: Real, Rv <: Real}
    T = promote_type(Ru, Rv, Float32)
    r2 = u^2 + v^2

    (sϕ, cϕ) = sincos(ϕ)
    (sθ, cθ) = sincos(θ)

    p1 = u * cϕ + v * sϕ * sθ
    p2 = u * sϕ - v * cϕ * sθ
    p3 = v * cθ

    e1 = -sϕ * cθ # x = p1 + ℓ * e1
    e2 = cϕ * cθ  # y = p2 + ℓ * e2
    e3 = sθ       # z = p3 + ℓ * e3

    ℓxmin, ℓxmax = cube_bounds(p1, e1)
    ℓymin, ℓymax = cube_bounds(p2, e2)
    ℓzmin, ℓzmax = cube_bounds(p3, e3)

    minℓ = max(ℓxmin, ℓymin, ℓzmin)
    maxℓ = min(ℓxmax, ℓymax, ℓzmax)
    ℓ = max(maxℓ - minℓ, zero(T))
    x = p1 + ℓ * e1
    y = p2 + ℓ * e2
    z = p3 + ℓ * e3

    if(hs.axiss==1)
        if(y<0)
            return zero(T)
        end
    elseif(hs.axiss==2)
        if(x<0)
            return zero(T)
        end
    elseif(hs.axiss==3)
        if(z<0)
            return zero(T)
        end
    end    
    if r2 < 1
        return  sqrt(one(T) - r2)
    end

    return zero(T)

end

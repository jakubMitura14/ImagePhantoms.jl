# src/half_sphere.jl
#=
half_sphere.jl
=#
export HalfSphere, half_sphere

"""
    HalfSphere <: AbstractShape{3}
"""
struct HalfSphere <: AbstractShape{3} end

# constructor
"""
    half_sphere(cx, cy, cz, r, Φ=0, Θ=0, value::Number = 1)
    half_sphere(center::NTuple{3,RealU}, radius::RealU, angle::NTuple{3,RealU}, v)
Construct `Object{HalfSphere}` from parameters.
"""
half_sphere(args... ; kwargs...) = Object(HalfSphere(), args...; kwargs...)

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
    ::HalfSphere,
    u::Ru,
    v::Rv,
    ϕ::RealU, # irrelevant
    θ::RealU, # irrelevant
) where {Ru <: Real, Rv <: Real}
    T = promote_type(Ru, Rv, Float32)
    r2 = u^2 + v^2
    if(ϕ<0)
        return zero(T)
    end

    if r2 < 1
        return  sqrt(one(T) - r2)
    end

    return zero(T)

end

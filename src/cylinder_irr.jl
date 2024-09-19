#=
cylinder_irr.jl
=#


export Cylinder, cylinder_irr


"""
    Cylinder_irr <: AbstractShape{3}
"""
struct Cylinder_irr <: AbstractShape{3} end


# constructor


"""
    cylinder_irr(cx, cy, cz, wx, wy, wz, Φ, Θ, value::Number)
    cylinder_irr(center::NTuple{3,RealU}, width::NTuple{3,RealU}, angle::NTuple{3,RealU}, v)
Construct `Object{Cylinder_irr}` from parameters;
here `width` is the *radius* in x,y and the *height* in z.
"""
cylinder_irr(args... ; kwargs...) = Object(Cylinder_irr(), args...; kwargs...)


# methods


volume1(::Cylinder_irr) = π

ℓmax1(::Cylinder_irr) = √5 # max line integral through unit cylinder_irr

ℓmax(ob::Object3d{Cylinder_irr}) = sqrt(sum(abs2, (2maximum(ob.width[1:2]), ob.width[3])))


"""
    phantom1(ob::Object3d{Cylinder_irr}, (x,y,z))
Evaluate unit cylinder_irr at `(x,y,z)`,
for unitless coordinates.
"""
phantom1(ob::Object3d{Cylinder_irr}, xyz::NTuple{3,Real}) =
    (abs(xyz[3]) ≤ 0.5) && (sum(abs2, xyz[1:2]) ≤ 1)


# radon


# line integral through rectangle of width (wx,wy) at (r,ϕ)
function _rect_proj(wx::Real, wy::Real, r::Real, ϕ::Real)
    (sϕ, cϕ) = sincos(ϕ)
    rp = sqrt((wx * cϕ)^2 + (wy * sϕ)^2) # projected radius
    dis = r / rp # scaled distance from center
    abs_cos_ang_pi = wx * abs(cϕ) / rp
    abs_sin_ang_pi = wy * abs(sϕ) / rp

    # break points of the trapezoid
    len = 1 / max(abs_cos_ang_pi, abs_sin_ang_pi)
    dmax = (abs_cos_ang_pi + abs_sin_ang_pi) / 2
    dbreak = abs(abs_cos_ang_pi - abs_sin_ang_pi) / 2
    dmax_break = dmax - dbreak
    scale = wx * wy / rp * len #lmax

    return scale * trapezoid(dis, -dmax, -dbreak, dbreak, dmax)
end


# x-ray transform (line integral) of unit cylinder_irr
# `u,v` should be unitless
function xray1(
    ::Cylinder_irr,
    u::Real,
    v::Real,
    ϕ::RealU, # azim (irrelevant)
    θ::RealU, # polar
)
    T = promote_type(typeof(u), typeof(v), Float32)

    r = abs(u)
    if r > 1
        return zero(T)
    end

    # rectangle in plane of distance `r` from origin
    wz = 1 # from unit-height of cylinder_irr
    wy = 2 * sqrt(1 - r^2)

    return T(_rect_proj(wz, wy, v, θ))+(randn()*0.1)
end


# spectrum

"""
    spectrum(::Object3d{Cylinder_irr}, (kx,ky,kz))
Spectrum of unit cylinder_irr at `(kx,ky,kz)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(::Object3d{Cylinder_irr}, kxyz::NTuple{3,Real})
    return 4 * jinc(2 * sqrt(sum(abs2, kxyz[1:2]))) * sinc(kxyz[3])
end

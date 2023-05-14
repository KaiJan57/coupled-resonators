using GLMakie

Base.@kwdef mutable struct Lorenz
    dt::Float64 = 0.01
    σ::Float64 = 10
    ρ::Float64 = 28
    β::Float64 = 8/3
    x::Float64 = 1
    y::Float64 = 1
    z::Float64 = 1
end

function step!(l::Lorenz)
    dx = l.σ * (l.y - l.x)
    dy = l.x * (l.ρ - l.z) - l.y
    dz = l.x * l.y - l.β * l.z
    l.x += l.dt * dx
    l.y += l.dt * dy
    l.z += l.dt * dz
    Point3f(l.x, l.y, l.z)
end

set_theme!(theme_black())
fig = Figure()
ax = Axis3(fig[1, 1], protrusions = (0, 0, 0, 0),
        viewmode = :fit, limits = (-30, 30, -30, 30, 0, 50))

slider = Slider(fig[1, 2], range = 0:0.01:10, horizontal = false, startvalue = 3)
label = Label(fig[2, 2], string(slider.value))

points = lift(slider.value) do σ
    label.text = string(σ)
    pts = Point3f[]
    attractor = Lorenz()
    attractor.σ = σ
    for i in 1:1000
        push!(pts, step!(attractor))
    end
    pts
end

lines!(points, color = :white, transparency = true)

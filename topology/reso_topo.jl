using DifferentialEquations
using GLMakie
using FFTW

function pt!(dα, α, p, t) # mappings for p: 1-N, 2-J, 3-U, 4-p, 5-δ, 6-κ
    N = Int(p[1])
    dα[1] = -(p[6]/2+im*p[5])*α[1]-im*p[2]*α[2]-2*im*p[3]*α[1]*conj(α[1])*α[1]
    for i in 2:N-1
        dα[i] = im*(p[5]*α[1]-p[2]*α[i+1]-p[2]*α[i-1]-2*p[3]*α[i]*conj(α[i])*α[i])
        if i == (N+1)/2
            dα[i] += -im*p[4]
        end
    end
    dα[N] = -(p[6]/2+im*p[5])*α[N]-im*p[2]*α[N-1]-2*im*p[3]*α[N]*conj(α[N])*α[N]
end


set_theme!(theme_black())
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "time t", ylabel = "field strength |α|", title = "Time evolution of field strength")

N=3
α0 = fill(1.0im, N)
tspan = (0.0, 100.0)
p = [N, 1.0, 1.0, 10.0, 0.0, 1.0]
prob = ODEProblem(pt!, α0, tspan, p);

s=solve(prob, RK4(), dt=0.001)
for i in 1:N
    lines!(s.t, abs.(s[i,:]), label = "Cavity nr. " * string(i), transparency = true)
end
axislegend(ax)

fig

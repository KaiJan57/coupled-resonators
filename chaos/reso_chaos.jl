using DifferentialEquations
using GLMakie
using FFTW

function pt!(dα, α, p, t) # mappings for p: 1-N, 2-J, 3-U, 4-p, 5-δ, 6-κ
    N = Int(p[1])
    dα[1] = -(p[6]/2+im*p[5])*α[1]-im*p[2]*α[2]-2*im*p[3]*α[1]*conj(α[1])*α[1]-im*p[4]
    for i in 2:N-1
        dα[i] = im*(p[5]*α[1]-p[2]*α[i+1]-p[2]*α[i-1]-2*p[3]*α[i]*conj(α[i])*α[i])
    end
    dα[N] = -(p[6]/2+im*p[5])*α[N]-im*p[2]*α[N-1]-2*im*p[3]*α[N]*conj(α[N])*α[N]
end


set_theme!(theme_black())
fig = Figure()
tdax = Axis(fig[1, 1], xlabel = "time t", ylabel = "field strength |α|", title = "Time evolution of field strength")
fftax = Axis(fig[1, 2], xlabel = "frequency f", ylabel = "field strength |α|", title = "FFT of field strength")

N=3
α0 = fill(1.0im, N)
tspan = (0.0, 100.0)
p = [N, 1.0, 1.0, 10.0, 0.0, 1.0]
prob = ODEProblem(pt!, α0, tspan, p);

samplerate = 1000
s=solve(prob, RK4(), dt=1/samplerate)
for i in 1:N
    lines!(tdax, s.t, abs.(s[i,:]), label = "Cavity nr. " * string(i), transparency = true)
end
F=fft(abs.(s[1,:])) |> fftshift
F[round(Int, size(F)[1]/2)]=0 # remove unwanted 0 Hz peak
freqs = fftfreq(length(s.t), samplerate) |> fftshift

lines!(fftax, freqs, abs.(F)/maximum(abs.(F)), label = "Spectrum of α_1", transparency = true)

axislegend(tdax)
axislegend(fftax)

fig

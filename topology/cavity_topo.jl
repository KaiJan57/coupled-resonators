using DifferentialEquations
using GLMakie
using FFTW

function pt!(dα, α, p, t) # mappings for p: 1-N, 2-J1, 3-J2, 4-U, 5-p, 6-δ, 7-κ
    N = Int(p[1])
    dα[1] = -(p[7]/2+im*p[6])*α[1]-im*p[2]*α[2]-2*im*p[4]*α[1]*conj(α[1])*α[1]#-im*p[5]
    for i in 2:N-1
        if iseven(i)
            dα[i] = im*(-p[6]*α[i]-p[3]*α[i+1]-p[2]*α[i-1]-2*p[4]*α[i]*conj(α[i])*α[i])-p[7]/2*α[i]
        else 
            dα[i] = im*(-p[6]*α[i]-p[2]*α[i+1]-p[3]*α[i-1]-2*p[4]*α[i]*conj(α[i])*α[i])-p[7]/2*α[i]
        end
        if i == 5
            dα[i] += -im*p[5]
        end
    end
    dα[N] = -(p[7]/2+im*p[6])*α[N]-im*p[3]*α[N-1]-2*im*p[4]*α[N]*conj(α[N])*α[N]
end


set_theme!(theme_black())
figt = Figure()
figc = Figure()
tax = Axis(figt[1, 1], xlabel = "time t", ylabel = "field strength |α|", title = "Time evolution of field strength")
cax = Axis(figc[1, 1], xlabel = "Cavity", ylabel = "field strength |α|", title = "Evolution of field strength over the cavities")

N=10
α0 = fill(1.0im, N)

#=x = rand()
y = rand()
α0 = [x+y*im]
for i in 2:N
    x1 = rand()
    y1 = rand()
    push!(α0, x1+y1*im)
end=#

tspan = (0.0, 200.0)
p = [N, 1.0, 8.0, 1.0, 10.0, 0.0, 0.25]
dt_save = 0.001
prob = ODEProblem(pt!, α0, tspan, p, saveat = dt_save);

s=solve(prob, RK4())
for i in 1:N
    #val = abs.(s[i,:]).*abs.(s[i,:])
    lines!(tax, s.t, abs.(s[i,:]), label = "Cavity nr. " * string(i), transparency = true)
end
#lines!(s.t, abs.(s[3,:]), label = "Cavity nr. 3", transparency = true)
cavity = range(1,N)
res = [abs.(s[1, size(s[1,:])[1]])]
for i in 2:N
 push!(res, abs.(s[i, size(s[i,:])[1]]))
end

lines!(cax, cavity, res, label = "Limit value", transparency = true)

cavityA = [1]
resA = [res[1]]
cavityB = [2]
resB = [res[2]]
for i in 3:N
    if (isodd(i))
        push!(cavityA, i)
        push!(resA, res[i])
    else
        push!(cavityB, i)
        push!(resB, res[i])
    end
end

scatter!(cax, cavityA, resA, label = "Limit value for cavity A", transparency = true, markersize = 15, color = :red)
scatter!(cax, cavityB, resB, label = "Limit value for cavity B", transparency = true, markersize = 15, color = :blue)

#println(s[2, :])
#println(s[2])
#println(size(s[2,:]))
#println(s[2, size(s[2,:])[1]])

axislegend(tax)
axislegend(cax)

display(figt)

println("Key")
readline()

display(figc)

println("Key")
readline()

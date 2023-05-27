#para_space_3d
using Distributed
using DifferentialEquations
using Plots
using Random, Statistics
using Dates
using JLD2
using Printf
using PyCall
using MATLAB

rng = MersenneTwister(1234)

function pt!(dα, α, para, t) # mappings for p: 1-N, 2-J1, 3-J2, 4-U, 5-p, 6-δ, 7-κ, 8-Np
    N = para[1]
    J1 = para[2]
    J2 = para[3]
    U = para[4]
    p = para[5]
    δ = para[6]
    κ = para[7]
    Np = para[8]

    dα[1] = -(im*δ)*α[1]-im*J1*α[2]-2*im*U*α[1]*conj(α[1])*α[1]-κ/2*α[1]#-im*p
    for i in 2:N-1
        if iseven(i)
            dα[i] = im*(-δ*α[i]-J2*α[i+1]-J1*α[i-1]-2*U*α[i]*conj(α[i])*α[i])-κ/2*α[i]
        else 
            dα[i] = im*(-δ*α[i]-J1*α[i+1]-J2*α[i-1]-2*U*α[i]*conj(α[i])*α[i])-κ/2*α[i]
        end
        if i == Np
            dα[i] = dα[i] -im*p
        end
    end
    dα[N] = -(im*δ)*α[N]-im*J2*α[N-1]-2*im*U*α[N]*conj(α[N])*α[N]-κ/2*α[N]
end

function calc(jj,kk)
    U = y_vec[kk]

    δ = x_vec[jj]

    para = (N, J1, J2, U, p, δ, κ, Np)
    dt_save = 0.001
    prob = ODEProblem(pt!, α0, tspan, para, saveat = dt_save);

    s=solve(prob, RK4())
    amplitude_alpha_N[jj,kk] = abs.(s[N,size(s.u)[1]])/(abs.(s[Np,size(s.u)[1]]))
    #println(amplitude_alpha_N[jj,kk])
end

x_min = -10
x_max = 10
number_pointsx = 100
x_vec = LinRange(x_min, x_max, number_pointsx)
x_vec_length = length(x_vec)

y_min = 10^(-6)
y_max = 10
number_pointsy = 100
y_vec = LinRange(y_min, y_max, number_pointsy)
y_vec_length = length(y_vec)

N = 10

α0 = randn(rng, ComplexF32, (1, N))/10


tspan = (0.0, 200.0)
J1 = 1.0
J2 = 8.0
#U = 1.0
p = 10.0
#δ  = 0.0
κ  = 1.0
Np = 5

amplitude_alpha_N = zeros((x_vec_length, y_vec_length))
#println(amplitude_alpha_N[1,1])

#=for jj = 1:x_vec_length
    for kk = 1:y_vec_length
    
        U = y_vec[kk]

        δ = x_vec[jj]

        para = (N, J1, J2, U, p, δ, κ, Np)
        dt_save = 0.001
        prob = ODEProblem(pt!, α0, tspan, para, saveat = dt_save);

        s=solve(prob, RK4())
        amplitude_alpha_N[jj,kk] = abs.(s[N,size(s.u)[1]])/(abs.(s[Np,size(s.u)[1]]))
        println(amplitude_alpha_N[jj,kk])
    end
end=#


pmap(jj -> pmap(kk -> calc(jj,kk), 1:y_vec_length), 1:x_vec_length)

z_mat = amplitude_alpha_N
#println(amplitude_alpha_N)
#pcolor(x_vec, y_vec, z_mat, cmap="inferno")
plt = plot(x_vec, y_vec, z_mat, yaxis=:log, st=:surface)
savefig(plt, "plt-3d-test")
display(plt)
println("Key")
readline()


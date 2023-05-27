using Distributed
using DifferentialEquations
#using GLMakie
using Plots
using Random, Statistics
using Dates
using JLD2
using Printf

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

function trajectories(traj)
    traj = Int(traj)
    
    #set_theme!(theme_black())
    #fig = Figure()

    #tax = Axis(fig[1, 1], xlabel = "time t", ylabel = "field strength |α|", title = "Time evolution of field strength " *string(traj))
    #cax = Axis(fig[1, 2], xlabel = "Cavity", ylabel = "field strength |α|", title = "Evolution of field strength over the cavities")


    N=10


    α0 = randn(rng, ComplexF32, (1, N))/10

    tspan = (0.0, 200.0)
    J1 = 1.0
    J2 = 8.0
    U = 1.0
    p = 10.0
    δ  = 0.0
    κ  = 1.0
    Np = 5
    para = (N, J1, J2, U, p, δ, κ, Np)
    dt_save = 0.001
    prob = ODEProblem(pt!, α0, tspan, para, saveat = dt_save);

    s=solve(prob, RK4())
    #=for i in 1:N
        lines!(tax, s.t, abs.(s[i,:]), label = "Cavity nr. " * string(i), transparency = true)
    end=#

    cavity = range(1,N)
    res = [abs.(s[1, size(s[1,:])[1]])]
    for i in 2:N
    push!(res, abs.(s[i, size(s[i,:])[1]]))
    end

    #lines!(cax, cavity, res, label = "Limit value", transparency = true)

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

    #scatter!(cax, cavityA, resA, label = "Limit value for cavity A", transparency = true, markersize = 15, color = :red)
    #scatter!(cax, cavityB, resB, label = "Limit value for cavity B", transparency = true, markersize = 15, color = :blue)

    #axislegend(tax)
    #axislegend(cax)

    #display(fig)
    plt = plot(s.t, abs.(s[1,:]), title = "Time evolution of field strength" *string(traj), xlabel = "time t", ylabel = "field strength |α|" *string(1))
    for i in 2:N
        plot!(plt, s.t, abs.(s[i,:]), title = "Time evolution of field strength " *string(traj), xlabel = "time t", ylabel = "field strength |α|" *string(i))
    end
    savefig("plt"*string(traj)*".png")
    println("Key")
    readline()

    
    folder_path = joinpath("Data", "BH_study_","")
    plt2 = plot(s.t, abs.(s[N,:]), title = "Time evolution of field strength" *string(traj), label = "field strength |α| at" *string(N))
    savefig(plt2, folder_path*"plot_"*string(traj))
    #jldsave(folder_path*"data_edge"*string(traj)*".txt", 42)

    dtt = Dates.format(Dates.now(),"yyyy-mm-dd_HH-MM-SS")
    jldsave(joinpath(folder_path*"sol_vec_$(dtt).jld2"); s)

    io = open(folder_path*"trajectory"*string(traj)*".tex", "w")
    #io = open(output_file, "w")
    for time in 1:size(s.u)[1]
        @printf(io, "%.2f %.2f\n", (time-1)*dt_save, abs.(s[N, time]))
    end

    close(io)
end

array = Int[1,2]
pmap(traj -> trajectories(traj), array)




    using Distributed
    using DifferentialEquations
    using Plots
    using Random
    using Statistics
    #using Dates
    #using JLD2
    #using Printf
    #using PyCall
    #using MATLAB

    #addprocs(2)

    rng = MersenneTwister(1234)

    function pt!(dα, α, para, t) # mappings for p: 1-N, 2-J1, 3-J2, 4-U, 5-p, 6-δ, 7-κ, 8-Np

        N, J1, J2, U, p, δ, κ, Np = para

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

    function calc(jj)
        #p = y_vec[kk]

        δ = x_vec[jj]

        para = (N, J1, J2, U, p, δ, κ, Np)
        dt_save = 0.1
        prob = ODEProblem(pt!, α0, tspan, para, saveat = dt_save);

        s=solve(prob, RK4())
        #amplitude_alpha_N[jj,kk] = abs.(s[N,size(s.u)[1]])/(abs.(s[Np,size(s.u)[1]]))
        #println(jj, kk)
        amplitude_alpha_N_cut[jj] = abs.(s[N,size(s.u)[1]])/(abs.(s[Np,size(s.u)[1]]))
        println(jj)
        #println(amplitude_alpha_N[jj,kk])
    end


    x_min = -10
    x_max = 10
    number_pointsx = 100
    x_vec = LinRange(x_min, x_max, number_pointsx)
    x_vec_length = length(x_vec)

    y_vec =  collect(10 .^ (range(start = -4, length = 101, stop = 2)))
    #=y_min = 0.1
    y_max = 10
    number_pointsy = 101
    y_vec = LinRange(y_min, y_max, number_pointsy)=#
    y_vec_length = length(y_vec)

    N = 30 #cavities

    #α0 = randn(rng, ComplexF32, (1, N))/10
    α0 = zeros(ComplexF32, 1,N)

    tspan = (0.0, 500.0)
    J1 = 1.0
    J2 = 8.0
    U = 10^(-1)
    p = 10^(1)
    #δ  = 0.0
    κ  = 1.0
    Np = 5

    
    amplitude_alpha_N = zeros(x_vec_length, y_vec_length)
    amplitude_alpha_N_cut = zeros(x_vec_length)
    #println(amplitude_alpha_N[1,1])


    #=for jj = 1:x_vec_length
        for kk = 1:y_vec_length
            calc(jj, kk)
        end
    end=#

    for jj = 1:x_vec_length
        calc(jj)
    end

    #pmap(jj -> pmap(kk -> calc(jj,kk), 1:y_vec_length), 1:x_vec_length)
    #=@sync @distributed for jj in 1:x_vec_length
        @sync @distributed for kk in 1:y_vec_length
            calc(jj,kk)
        end
    end=#
    println(x_vec)
    println(log10.(amplitude_alpha_N_cut))

    #println(amplitude_alpha_N)
    #pcolor(x_vec, y_vec, z_mat, cmap="inferno")
    #plt = heatmap(x_vec, y_vec, log10.(amplitude_alpha_N)', title = "α_{N}/α_{Np}, N="*string(N)*", δ="*string(δ)*", U="*string(U)*", p="*string(p)*", κ="*string(κ)*", Np="*string(Np)*", T="*string(tspan[2])*"", titlefont = font(11,"Computer Modern"), ylabel = "J2", xlabel = "J1", zlabel = "α_{N}/α_{Np}")#, title="|α1|\n$(Parameters1)")
    #plt = heatmap(x_vec, y_vec, log10.(amplitude_alpha_N)', yscale=:log10, title = "α_{N}/α_{Np}, N="*string(N)*", J1="*string(J1)*", J2="*string(J2)*", p="*string(p)*", κ="*string(κ)*", Np="*string(Np)*", T="*string(tspan[2])*"", titlefont = font(11,"Computer Modern"), ylabel = "U", xlabel = "δ", zlabel = "α_{N}/α_{Np}")#, title="|α1|\n$(Parameters1)")
    #plt = heatmap(x_vec, y_vec, log10.(amplitude_alpha_N)', yscale=:log10, title = "α_{N}/α_{Np}, N="*string(N)*", J1="*string(J1)*", J2="*string(J2)*", U="*string(U)*", κ="*string(κ)*", Np="*string(Np)*", T="*string(tspan[2])*"", titlefont = font(11,"Computer Modern"), ylabel = "p", xlabel = "δ", zlabel = "α_{N}/α_{Np}")#, title="|α1|\n$(Parameters1)")
    plt = plot(x_vec, log10.(amplitude_alpha_N_cut), title = "α_{N}/α_{Np}, N="*string(N)*#=", δ ="*string(δ)*=#", J1="*string(J1)*", J2="*string(J2)*", U="*string(U)*", p="*string(p)*", κ="*string(κ)*", Np="*string(Np)*", T="*string(tspan[2])*"", titlefont = font(10,"Computer Modern"),  xlabel = "δ", ylabel = "α_{N}/α_{Np}")


    #xlims!(1e-6,1e+1)
    savefig(plt, "versus-plots/cuts/d-U10(-1)-cut-N30.png")


    display(plt)
    println("Key")
    readline()

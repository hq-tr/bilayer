include("/home/trung/_qhe-julia/FQH_bilayer.jl")
using .BilayerFQH
using BenchmarkTools
using ArgMacros
using SparseArrays
using Arpack # For sparse matrix ED
using LinearAlgebra

function ground_state(basis_line::Vector{BitVector}, z_c::ComplexF64,ϵ::Float64)
    H =  10I - 100*bilayer_density_matrix(basis_line,ϵ, z_c)
    E,V = eigs(H, nev=8, which=:SM)
    return V[:,1]
end

function berry_phase(R::Float64, basis_line::Vector{BitVector},ϵ::Float64)
    θ_step  = minimum((0.0005/R,0.01)) # make number of points scale with R
    θ_list  = 2π * (0:θ_step:1)
    gs_list = Vector{ComplexF64}[]

    for θ in θ_list
        push!(gs_list, ground_state(basis_line,R*exp(im*θ),ϵ))
    end

    push!(gs_list,gs_list[1]) # "Periodic boundary"
    co = prod(i->gs_list[i]⋅gs_list[i+1], 1:length(θ_list))
    println((abs(co), angle(co)))
    return abs(co), angle(co), gs_list[1]
end


function main()
    @inlinearguments begin
        @argumentrequired String fname "-f" "--filename"
        @argumentrequired Float64 epsilon "-e" "--epsilon"
    end

    #epsilon = 0.3

    if !isfile(fname)
        println("$(fname): File does not exist. Terminating.")
        return
    end
    #fname = "one_particle_basis"
    println("Reading basis from file [$(fname)]")

    basis = bilayerreadwf(fname).basis
    dim = length(basis[1])
    basis_line = [vcat(basis[1][i], basis[2][i]) for i in 1:dim]
    #display(basis_line[1])

    if !isdir("ground_states")
        mkdir("ground_states")
    end

    if !isdir("berry_phase")
        mkdir("berry_phase")
    end

    if !isdir("eigenvalues")
        mkdir("eigenvalues")
    end


    N_o    = length(basis[1][1]) 
    R_max  = √(2N_o)
    R_step = 0.15


    R_list  = 0.1:R_step:R_max
    D_list  = Float64[]
    γ_list  = Float64[]

    Lz_list = Float64[]

    for R in R_list
        println("====== R = $R")
        @time D, γ ,gs_coef= berry_phase(R, basis_line, epsilon)
        push!(D_list,D)
        push!(γ_list,γ)
        gs_wf = bilayer_state(basis, gs_coef)

        push!(Lz_list,findLz(gs_wf,epsilon))
        bilayerprintwf(gs_wf;fname="ground_states/$(fname)_R_$(R)_ϵ_$(epsilon)")
    end

    open("berry_phase/$(fname)_one_pin_ϵ_$(epsilon).txt","w+") do f
        for i in 1:length(R_list)
            write(f,"$(R_list[i])\t$(D_list[i])\t$(γ_list[i])\n")
        end
    end

    open("berry_phase/$(fname)_one_pin_Lz_ϵ_$(epsilon).txt","w+") do f
        for i in 1:length(R_list)
            write(f,"$(R_list[i])\t$(Lz_list[i])\n")
        end
    end
    return 
end



@time main()
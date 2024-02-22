include("/home/trung/_qhe-julia/FQH_bilayer.jl")
using .BilayerFQH
using BenchmarkTools
using ArgMacros
using SparseArrays
using Arpack # For sparse matrix ED
using LinearAlgebra

function main()
    @inlinearguments begin
        @argumentrequired String fname "-f" "--filename"
    end

    if !isfile(fname)
        println("$(fname): File does not exist. Terminating.")
        return
    end
    #fname = "one_particle_basis"
    println("Reading basis from file [$(fname)]")

    println("Apply the same pins to both layers (1) or separate pins for separate layers (2)?")
    popt = parse(Int, readline())

    if popt == 2
        println("How many pins in total? ")
        npins = parse(Int, readline())

        x_list = Float64[]
        y_list = Float64[]
        λ_list = Int[]

        for i in 1:npins
            println("Input x-position of potential pin: ")
            x = parse(Float64,readline())
            println("Input y-position of potential pin: ")
            y = parse(Float64,readline())
            println("Apply this pin to first layer (1), second layer (2) or both(3)?")
            λ = parse(Int, readline())
            push!(x_list, x)
            push!(y_list, y)
            push!(λ_list, λ)
        end

        posstring = join(["($(x_list[i]),$(y_list[i]),$(λ_list[i]))" for i in 1:length(x_list)],"_")
        outdir = "$(fname)_output_$(posstring)"
        if !isdir(outdir)
            mkdir(outdir)
        end

        z_c = 0.5 .* (x_list + im * y_list)

        basis = bilayerreadwf(fname).basis
        dim = length(basis[1])
        basis_line = [vcat(basis[1][i], basis[2][i]) for i in 1:dim]
        #display(basis_line[1])

        println("-------")
        println("Constructing matrix")

        dim = length(basis_line)
        mat = spzeros(ComplexF64, (dim,dim))

        @time begin
            for i in 1:length(z_c)
                mat += bilayer_density_matrix(basis_line, 0.5,z_c[i];layer_choice = λ_list[i])
            end

            H =  10I - 10*mat
        end

        if !ishermitian(H) 
            println("WARNING: Hamiltonian is not Hermitian! Check for bugs!") 
        end

        println("-------")
        E,V = eigs(H, nev=8, which=:SM)
        println("Eigenvalues:")
        display(real.(E))

        open("$(outdir)/eigen.txt","w") do f
            for ε in E
                write(f,"$ε\n")
            end
        end

        for i in 1:length(E)
            gs = bilayer_state(basis, V[:,i])
            bilayerprintwf(gs;fname="$(outdir)/g_$(i-1)")
        end
    elseif popt == 1
        println("How many pins? ")
        npins = parse(Int, readline())

        x_list = Float64[]
        y_list = Float64[]

        for i in 1:npins
            println("Input x-position of potential pin: ")
            x = parse(Float64,readline())
            println("Input y-position of potential pin: ")
            y = parse(Float64,readline())
            push!(x_list, x)
            push!(y_list, y)
        end

        posstring = join(["($(x_list[i]),$(y_list[i]))" for i in 1:length(x_list)],"_")
        outdir = "$(fname)_output_$(posstring)"
        if !isdir(outdir)
            mkdir(outdir)
        end

        z_c = 0.5 .* (x_list + im * y_list)

        basis = bilayerreadwf(fname).basis
        dim = length(basis[1])
        basis_line = [vcat(basis[1][i], basis[2][i]) for i in 1:dim]
        #display(basis_line[1])

        println("-------")
        println("Constructing matrix")

        dim = length(basis_line)
        mat = spzeros(ComplexF64, (dim,dim))

        @time begin
            for z in z_c
                mat += bilayer_density_matrix(basis_line, 0.5,z)
            end

            H =  10I - 10*mat
        end

        if !ishermitian(H) 
            println("WARNING: Hamiltonian is not Hermitian! Check for bugs!") 
        end

        println("-------")
        E,V = eigs(H, nev=8, which=:SM)
        println("Eigenvalues:")
        display(real.(E))

        open("$(outdir)/eigen.txt","w") do f
            for ε in E
                write(f,"$ε\n")
            end
        end

        for i in 1:length(E)
            gs = bilayer_state(basis, V[:,i])
            bilayerprintwf(gs;fname="$(outdir)/g_$(i-1)")
        end
    else
        return 
    end
end



@time main()
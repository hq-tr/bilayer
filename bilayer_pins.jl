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
println("Input x-position of potential pin: ")
x = parse(Float64,readline())
println("Input y-position of potential pin: ")
y = parse(Float64,readline())

z_c = x+im*y

basis = bilayerreadwf(fname).basis
dim = length(basis[1])
basis_line = [vcat(basis[1][i], basis[2][i]) for i in 1:dim]
#display(basis_line[1])

println("-------")
println("Constructing matrix")

@time H =  10I - 100*bilayer_density_matrix(basis_line, 0.5, z_c)
if !ishermitian(H) 
    println("WARNING: Hamiltonian is not Hermitian! Check for bugs!") 
end

println("-------")
E,V = eigs(H, nev=8, which=:SM)
println("Eigenvalues:")
display(real.(E))
display(V[:,1])

for i in 1:8
    gs = bilayer_state(basis, V[:,i])
    bilayerprintwf(gs;fname="g_$(i-1)")
end

return 
end



@time main()
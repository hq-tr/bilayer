include("/home/hqtr/QHE_Julia/FQH_bilayer.jl")
using .BilayerFQH
using BenchmarkTools

println("Plot the density of single-particle coherent state centered at position (x,y)")

print("Input number of orbital: "); No = parse(Int64,readline())
print("Input x: "); x = parse(Float64, readline())
print("Input y: "); y = parse(Float64, readline())
print("Half-flux shift? (0 for no, 1 for yes) "); S = parse(Bool, readline())

@time begin
# Construct the states
basis_single = map(i -> BitVector(0:(No-1) .== i), 0:(No-1))
zero_vector  = BitVector(zeros(No)); zero_vector[1]=1
if S
	basis = map(v ->[zero_vector, v], basis_single)
else
	basis = map(v ->[v,zero_vector], basis_single)
end
a = (x + im*y)/sqrt(2)
coef = map(i -> a^i/factorial(big(i)),0:(No-1))

mystate = bilayer_state(basis, coef)
#printwf(mystate)


# Get the density
disk_density(mystate,"coherent_$(No)o_$(x)_$(y)_$(S)")
end

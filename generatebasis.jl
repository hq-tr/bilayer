include("/home/hqtr/QHE_Julia/HilbertSpace.jl")
include("/home/hqtr/QHE_Julia/FQH_bilayer.jl")
using .BilayerFQH
using .HilbertSpaceGenerator

basis = bilayerhilbertspace(0,1,30)
coef  = zeros(length(basis[1]))

state = bilayer_state(basis,coef)

bilayerprintwf(state;fname="basis_layer2")

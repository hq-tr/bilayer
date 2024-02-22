include("/home/hqtr/QHE_Julia/FQH_state_v2.jl")
include("/home/hqtr/QHE_Julia/FQH_bilayer.jl")


state = Main.BilayerFQH.bilayerreadwf("g_0")
basis = state.basis[1]

state2 = Main.FQH_states.FQH_state(basis, state.coef)
Main.FQH_states.printwf(state2;fname="g_0_lay1")

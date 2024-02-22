include("/home/hqtr/QHE_Julia/FQH_bilayer.jl")
include("/home/hqtr/QHE_Julia/FQH_state_v2.jl")
using .FQH_states
using .BilayerFQH

using BenchmarkTools
using ArgMacros
using Plots

ENV["GKSwstype"] = "100" # Disconnect from Display

function main()
    @inlinearguments begin
        @argumentrequired String fname "-f" "--filename"
    end

ϵ = 0.5

# NOTE: The file is read as a single-layer state anyway
# It is easier to feed into function
state = readwf(fname)
No = length(state.basis[1]) ÷ 2
R_max = sqrt(2*No)+0.2

println("Each layer has $(No) orbitals.")
println("Check norm = $(wfnorm(state))")

x = -R_max:0.05:R_max
y = -R_max:0.05:R_max

N = length(x)

Y = repeat(y, inner=(1,N))
X = collect(transpose(Y))

Z = 0.5(X+Y*im)
ϕ = single_particle_state_disk # Alias. Use ϕ(z, m ,Δ)
#single_particle_function = map(m->ϕ(Z, m-No*(m>No), ϵ*(m>No)),1:2*No)
single_particle_function= [ϕ.(Z, m-No*(m>=No), ϵ*(m>=No)) for m in 0:length(state.basis[1])]

println("---- test")
println(single_particle_function[3][8,8])
println("---------")

D = dim(state)

den = zeros(size(Z))
@time begin
for i in 1:D
    print("\r$i\t")
    for j in i:D
        print("\r$i\t$j\t\t\t")
        coef = 2^(i!=j) * conj(state.coef[i]) * state.coef[j]
        bilayer_density_element!(No, den, coef, state.basis[i], state.basis[j], single_particle_function)
    end
end
end

println("-----")
#den = @time get_density_disk(state, X, Y)

open("$(fname)_density.dat", "w") do f
    for i in 1:N, j in 1:N
        write(f, "$(X[i,j])\t$(Y[i,j])\t$(den[i,j])\n")
    end
end


#p = plot(heatmap(z=density, aspect_ratio=:equal)) 

p = heatmap(x, y, den, aspect_ratio=:equal)
savefig(p, "$(fname)_density.svg")
end

main()
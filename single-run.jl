include("main.jl")

# conf = []
# for i = 1:5
#     cf = main(true,false,"new-params-$i")
#     push!(conf,cf)
# end
# plot(1:T,conf,title="Longer Fibre, More Filaments, More Motors",xlabel="Time",ylabel="Contractile Force",legend=false)
# savefig("new-params-cf.pdf")

cf = main(true,false,"plot-test")
plot(1:T, cf)

# I = [10, 20, 30, 40]
# J = [4.0,6.0,8.0,10.0]
# for j in J
#     println("B = $j")
#     global B = j
#     for i in I
#         println("M = $i")
#         conf = []
#         global M = i
#         for k = 1:3
#             cf = main(false,false,"M-$i")
#             push!(conf,cf)
#         end
#         # plot(1:T,conf,title="η = $i and α = $j",xlabel="Time",ylabel="Contractile Force",labels=string.(I'),legend=:outerright)
#         plot(1:T,conf,title="M = $i and B = $j",xlabel="Time",ylabel="Contractile Force",legend=false)
#         savefig("M-$i-B-$j-cf.pdf")
#     end
# end


### Computational study plan
# Reduce friction to increase evolution speed
# Implement filament turnover, randomising filament movement and removing previously attached motors
# Longer stress fibre 10 μm (get reference)
# Work out how to check if stress fibre is broken (prevent fibre breaking with filament turnover)
# Do test runs for reference parameters, finding a way to evaluate the steady contractile force

# Take avg cf for last 50 timesteps after long time step, 40-50 trials

# Eqn for stress: Take derivative wrt length of stress fibre and so F=-dE/dL

# Loess for υ, α/β, η 40-50 trials for mean
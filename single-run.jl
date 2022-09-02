include("main.jl")

# conf = []
# for i = 1:5
#     cf = main(true,false,"filament-turnover-rate-$α-less-drag-$i")
#     push!(conf,cf)
# end
# plot(1:T,conf,title="Longer Stress Fibre",xlabel="Time",ylabel="Contractile Force",legend=false)
# savefig("filament-turnover-rate-$α-less-drag-cf.pdf")

conf = []
I = 10:10:50
for i in I
    global M = i
    cf = main(true,false,"M-$i")
    push!(conf,cf)
end
plot(1:T,conf,title="Number of Motors",xlabel="Time",ylabel="Contractile Force",labels=string.(I'),legend=:outerright)
savefig("M-cf.pdf")

### Computational study plan
# Reduce friction to increase evolution speed
# Implement filament turnover, randomising filament movement and removing previously attached motors
# Longer stress fibre 10 μm (get reference)
# Work out how to check if stress fibre is broken (prevent fibre breaking with filament turnover)
# Do test runs for reference parameters, finding a way to evaluate the steady contractile force

# Take avg cf for last 50 timesteps after long time step, 40-50 trials

# Eqn for stress: Take derivative wrt length of stress fibre and so F=-dE/dL
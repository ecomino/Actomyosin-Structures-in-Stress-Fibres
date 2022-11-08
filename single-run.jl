include("src/main.jl")

# Visualise simulation with reference parameters
main(Parameters(),true,false,"");

# Plot contractile force of stress fibres with early motor drop-offs and motor turnover
conf = []
for i in 1:3
    cf = main(Parameters(δ=1, β=0.01),false,false,"")
    push!(conf,cf)
end
plot(1:250,conf,xlabel="Time",ylabel="Contractile Force",legend=false)
savefig("example-figures/ref-cf.pdf")
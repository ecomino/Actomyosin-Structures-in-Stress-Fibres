include("main.jl")

num_trials = 5 # Number of trials per parameter value
param_vec = 10 .^ (-6:0.1:1) # Parameter test values to test
mean_cf_vec = Vector{Float64}(undef,length(param_vec)) # Mean steady contractile force for each parameter test value
filename = "test" # Label parameter test
param_str = "eta" # Parameter label

for (i,param) in enumerate(param_vec)
    println("---Parameter = $(round(param,digits=6))---")
    global η = param # Update parameter
    contractile_force_trials = Vector{Vector{Float64}}(undef, num_trials)
    
    for j in 1:num_trials
        println("--Trial $j--")
        global anim = Animation()
        global contractile_force_trials[j] = main(false,false,"$(filename)-$i-$j")
    end

    plot(1:T,contractile_force_trials,title="Stress Fibre Contraction",xlabel="Time",ylabel="Contractile Force",legend=false)
    savefig("$(filename)-$(param_str)-$(round(param,digits=6))-trials.pdf")

    steady_contractile_forces = maximum.([abs.(k[.!isnan.(k)]) for k in contractile_force_trials])
    mean_cf_vec[i] = mean(steady_contractile_forces)

end

# Plot LOESS Regression with confidence intervals
loess_plot(param_vec,mean_cf_vec,"Drag due to Crosslinkers","loess-$(filename)-$(param_str)-$(num_trials)-trials-$(length(param_vec))-points.pdf")
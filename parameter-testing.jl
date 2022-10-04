include("main.jl")

# Eta
num_trials = 25 # Number of trials per parameter value
param_vec = 10 .^ (-4:0.1:1) # Parameter test values to test
# param_vec = collect(0.0:0.025:L/2)
mean_cf_vec = Vector{Float64}(undef,length(param_vec)) # Mean steady contractile force for each parameter test value
filename = "loess-eta/v1" # Label parameter test
param_str = "eta" # Parameter label

for (i,param) in enumerate(param_vec)
    println("---Parameter = $(round(param,digits=6))---")
    global η = param # Update parameter
    contractile_force_trials = Vector{Vector{Float64}}(undef, num_trials)
    
    for j in 1:num_trials
        print(" Trial $j...")
        global contractile_force_trials[j] = main(false,false,"$(filename)-$i-$j")
        println("finished")
    end

    plot(1:T,contractile_force_trials,title="Stress Fibre Contraction",xlabel="Time",ylabel="Contractile Force",legend=false)
    savefig("$(filename)-$(param_str)-$(round(param,digits=6))-trials.pdf")

    steady_contractile_forces = mean.(contractile_force_trials)
    mean_cf_vec[i] = mean(steady_contractile_forces)
end


# Plot LOESS Regression with confidence intervals
loess_plot_log(param_vec,mean_cf_vec,"Drag Effective due to Crosslinking","$(filename)-loess-$(param_str)-$(num_trials)-trials-$(length(param_vec))-points.pdf")
println("Loess Complete. Success!")

# Alpha
# num_trials = 25 # Number of trials per parameter value
# param_vec = 10 .^ (-6:0.1:1) # Parameter test values to test
# # param_vec = collect(0.0:0.025:L/2)
# mean_cf_vec = Vector{Float64}(undef,length(param_vec)) # Mean steady contractile force for each parameter test value
# filename = "loess-alpha/v1" # Label parameter test
# param_str = "alpha" # Parameter label

# for (i,param) in enumerate(param_vec)
#     println("---Parameter = $(round(param,digits=6))---")
#     global η = param # Update parameter
#     contractile_force_trials = Vector{Vector{Float64}}(undef, num_trials)
    
#     for j in 1:num_trials
#         print(" Trial $j...")
#         global contractile_force_trials[j] = main(false,false,"$(filename)-$i-$j")
#         println("finished")
#     end

#     plot(1:T,contractile_force_trials,title="Stress Fibre Contraction",xlabel="Time",ylabel="Contractile Force",legend=false)
#     savefig("$(filename)-$(param_str)-$(round(param,digits=6))-trials.pdf")

#     steady_contractile_forces = mean.(contractile_force_trials)
#     mean_cf_vec[i] = mean(steady_contractile_forces)
# end

# # Plot LOESS Regression with confidence intervals
# loess_plot_log(param_vec,mean_cf_vec,"Filament Turnover Rate","$(filename)-loess-$(param_str)-$(num_trials)-trials-$(length(param_vec))-points.pdf")
# println("Loess Complete. Success!")

# Beta
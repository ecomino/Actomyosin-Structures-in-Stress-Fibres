include("src/main.jl")

### TESTING EARLY MOTOR DROP-OFFS ###
num_trials = 10 # Number of trials per parameter value
param_vec = collect(0.0:0.2:2) # Parameter test values to test
mean_cf_vec = Vector{Float64}(undef,length(param_vec)) # Mean steady contractile force for each parameter test value
param_str = "delta" # Parameter label
filename = "loess-v1-$(param_str)-$(num_trials)-trials-$(length(param_vec))-points"

for (i,param) in enumerate(param_vec)
    print("Parameter = $(round(param,digits=6)). Trials...")
    contractile_force_trials = Vector{Vector{Float64}}(undef, num_trials)
    
    for j in 1:num_trials
        print(" $j")
        contractile_force_trials[j] = main(Parameters(δ=param),false,false,"")
    end
    println("...finished")

    steady_contractile_forces = mean.(contractile_force_trials)
    mean_cf_vec[i] = mean(steady_contractile_forces)
end

# Plot LOESS Regression with confidence intervals
loess_plot(param_vec,mean_cf_vec,"Motor Drop-off Distance","example-figures/$(filename).pdf")
println("Loess Complete. Success!")

## COMPARING MOTOR DROP-OFFS AND FILAMENT TURNOVER ###
num_trials = 10
I = collect(0.0:0.2:2)
J = 10 .^ (-6:0.5:-1)
heat = Matrix{Float64}(undef,length(I),length(J))
for (ind_j, j) in enumerate(J)
    print("α = $j... δ = ")
    for (ind_i,i) in enumerate(I)
        print("$i,")
        contractile_force_trials = Vector{Vector{Float64}}(undef, num_trials)
        for k = 1:num_trials
            contractile_force_trials[k] = main(Parameters(α=j,δ=i),false,false,"")
        end
        steady_contractile_forces = mean.(contractile_force_trials)
        heat[ind_i,ind_j] = mean(steady_contractile_forces)
    end
    println("done.")
end
gr()
heatmap(I, J, heat', c=cgrad([:red, :white, :blue]), clims=(-0.9, 0.9).*maximum(abs, heat), xlabel="Motor Drop-Off Distance", ylabel="Filament Turnover",title="Combined Effect of α and δ on Contractile Force", yaxis=:log)
savefig("example-figures/alpha-delta-heat.pdf")
using Plots, Optim, DataFrames, CSV
using Distributions: Poisson, sample
using StatsBase: mean
include("model-parameters.jl")
include("model.jl")
include("plot-simulation.jl")
include("write-simulation.jl")

function main(PLOTSIM,WRITESIM,filename)
    # Set up simulation
    X, Y, contractile_force = reset_model()
    WRITESIM && (df = initialise_dataframe())

    # Evolve simulation
    for t in 1:T
        # State progress
        # (t % 50 == 0) && println("Iteration $t/$T")
        # Store current simulation status
        contractile_force[t] =  k * (X[t][focal_adhesions[2]] - X[t][focal_adhesions[1]] - (B - A)) # Calculate current contractile force between focal adhesions
        (WRITESIM && t > 1) && push!(df, sim_status(t,X[t],Y,contractile_force[t])) 
        PLOTSIM && plot_sim(X[t],Y,t)

        # Perform minimisation to advance simulation to next time step
        od = OnceDifferentiable(x -> E(x, X[t], O(X[t][filaments]), Oᵩ(X[t][filaments], X[t][focal_adhesions]), Y), X[t]; autodiff=:forward)
        X[t+1] = Optim.minimizer(optimize(od, X[t], LBFGS()))

        # Perform random filament turnover
        X[t+1] = turnover(X[t+1])
        # Update attached filaments to motors
        Y = update_af(X[t+1],Y)
    end

    # Output results
    # PLOTSIM && gif(anim, "$(filename).gif", fps=5)
    WRITESIM && CSV.write("$(filename).csv", df)
    return contractile_force[1:end-1]
end;
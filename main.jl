using Plots, Optim, DataFrames, CSV
using Distributions: Poisson, sample
using StatsBase: mean
include("model-parameters.jl")
include("model.jl")
include("plot-simulation.jl")
include("write-simulation.jl")

function main(params::Parameters,PLOTSIM::Bool=false,WRITESIM::Bool=false,filename::String="")
    # Set up simulation
    X, Y, contractile_force, params = reset_model(params)
    WRITESIM && (df = initialise_dataframe())

    # Evolve simulation
    for t in 1:params.T
        # Store current simulation status
        contractile_force[t] =  params.k * (X[t][params.focal_adhesions[2]] - X[t][params.focal_adhesions[1]] - (params.B - params.A)) # Calculate current contractile force between focal adhesions
        (WRITESIM && t > 1) && push!(df, sim_status(t,X[t],Y,contractile_force[t],params.filaments,params.motors,params.focal_adhesions)) 
        PLOTSIM && plot_sim(X[t],Y,t,params)
        PLOTSIM && (params = save_sim(params))

        # Perform minimisation to advance simulation to next time step
        od = OnceDifferentiable(x -> E(x, X[t], O(X[t][params.filaments],params.L), Oᵩ(X[t][params.filaments],X[t][params.focal_adhesions],params.L), Y, params), X[t]; autodiff=:forward)
        X[t+1] = Optim.minimizer(optimize(od, X[t], LBFGS()))

        # Perform random filament turnover
        X[t+1] = turnover(X[t+1],params)
        # Update attached filaments to motors
        Y = update_af(X[t+1],Y,params)
    end

    # Output results
    PLOTSIM && gif(params.anim, "$(filename).gif", fps=5)
    WRITESIM && CSV.write("$(filename).csv", df) 
    return contractile_force[1:end-1]
end;
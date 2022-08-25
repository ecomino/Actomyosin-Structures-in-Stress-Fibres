using Plots, Optim, DataFrames, CSV
using Distributions: Poisson, sample
using StatsBase: mean
include("model-parameters.jl")
include("plot-simulation.jl")
include("model.jl")

function main(PLOTSIM,WRITESIM,filename)
    X = [vcat(B .* rand(N), B .* rand(M), A, B)] # Centre point of N filaments, M motors and focal adhesions centred at end points [A,B]
    # X = [[0.25, 0.75, 0.5, A, B]]
    Y = [[] for m in 1:M] # Y[m]...List of fibres attached to motor m
    # Generate filament-motor connections
    for m in sample(1:M,λ,replace=false)
        Y[m] = gen_af(m, X[end])
    end
    # Y = [[1,2]]
    contractile_force = [] # Contractile force between focal adhesions at each time step
    WRITESIM && (df = DataFrame(Time=Int64[], FilamentPos=Vector{Vector{Float64}}(), MotorPos=Vector{Vector{Float64}}(), FocalTesionPos=Vector{Vector{Float64}}(), MotorConnections=String[], ContractileForce=Float64[], Other=String[]))
    # Evolve simulation
    for t in 1:T
        (t % 10 == 0) && println("Iteration $t/$T")
        cf = k * (X[end][focal_adhesions[2]] - X[end][focal_adhesions[1]] - (B - A)) # Calculate current contractile force between focal adhesions
        push!(contractile_force, cf)
        (WRITESIM && t!=1) && push!(df, (t, round.(X[end][filaments],digits=6), round.(X[end][motors],digits=6), round.(X[end][focal_adhesions],digits=6), serialize_motor_connections(Y), round(cf,sigdigits=6)))
        PLOTSIM && plot_sim(X[end],Y,t)
        od = OnceDifferentiable(x -> E(x, X[end], O(X[end][filaments]), Oᵩ(X[end][filaments], X[end][focal_adhesions]), Y), X[end]; autodiff=:forward)
        next_x = Optim.minimizer(optimize(od, X[end], LBFGS()))
        push!(X, next_x)
        Y = update_af(X[end],Y)
    end
    # Output results
    PLOTSIM && gif(anim, "$(filename).gif", fps=5)
    WRITESIM && CSV.write("$(filename).csv", df)
    return contractile_force
end;

anim = Animation()

contractile_force_trials = Vector{Vector{Float64}}()
# global P = [1,-1]
U = collect(0:0.25:L/2)
for u in U
    println("---Run $u---")
    global anim = Animation()
    global υ = u
    # global P = rand((-1,1),N)
    con = main(true,false,"test")
    push!(contractile_force_trials,con)
end
plot(1:T,contractile_force_trials,title="Many Motor Varied Dropoff",xlabel="Time",ylabel="Contractile Force",legend=false)
# savefig("sim-ex.png")

# Threading
# U = collect(0:0.1:L/2)
# conf = Vector{Vector{Float64}}(undef,length(U))
# Threads.@threads for i in eachindex(U)
#     println("---Run $i---")
#     global anim = Animation()
#     global υ = U[i]
#     # global P = rand((-1,1),N)
#     conf[i] = main(true,false,"threads-$(U[i])")
# end
# Turnover, random filament movement (remove attached motors)
# Longer stress fibre 10 μm

# Computational study
# Mean and variance of contractile force, with reference parameters
# Then reduce/increase parameters (e.g. η vs cf)
# Loess (regression) fit, 99% ci, R package spatialeco 
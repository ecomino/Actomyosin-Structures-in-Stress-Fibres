using LinearAlgebra: eigvals

# Compute adjacency matrix of overlap matrix
adj(O) = (O).!==(0.0)
# Compute degree matrix of overlap matrix
deg(O) = [i == j ? Float64(sum(O[i,:].!==(0.0))) : 0.0 for i in 1:size(O)[1], j in 1:size(O)[2]]
# Compute Laplacian matrix for list of filaments
function laplace(x)
    o = O(x)
    L = deg(o) - adj(o)
    return L
end
# Determine whether filaments are disconnected (and hence if fibre is broken), return true if structure is disconnected
function is_broken(x)
    # Simple check if one fibre is disconnected
    0.0 ∈ sum(O(x),dims=1) && return true
    # If the second smallest eigenvalue of Laplacian matrix (a.k.a. Fiedler value) is greater than 0, the graph is connected
    fiedler = eigvals(laplace(x))[2]
    status = fiedler > 0.0 ? false : true
    return status
end

function main(PLOTSIM,WRITESIM,filename)
    df = DataFrame(Time=Int64[], FilamentPos=Vector{Vector{Float64}}(), MotorPos=Vector{Vector{Float64}}(), FocalTesionPos=Vector{Vector{Float64}}(), MotorConnections=Vector{Vector{Vector{Int64}}}(), ContractileForce=Float64[], Other=String[])
    X = [vcat(B .*rand(N), B .*rand(M),A,B)] # Centre point of N filaments, M motors and focal tesions centred at end points [A,B]
    # Ensure no initial breakage
    count = 0
    while is_broken(X[end][structure])
        count > 10 && error("Could not find valid initial structure. Try adding more filaments, or shortening the stress fibre.")
        count += 1
        X = [vcat(B .*rand(N), B .*rand(M),A,B)]
    end
    println("Accepted initial structure")
    display(adj(O(X[end][structure])))
    display(deg(O(X[end][structure])))
    display(laplace(X[end][structure]))
    println("Has fiedler value: ", eigvals(laplace(X[end][structure]))[2])
    Y = [Int64[] for _ in 1:M] # Y[m]...List of fibres attached to motor m
    # Generate filament-motor connections
    for m in sample(1:M,λ,replace=false)
        Y[m] = gen_af(m, X[end])
    end
    contractile_force = [] # Contractile force between focal tesions
    # Evolve simulation
    for t in 1:T
        println("---Iteration $t/$T---")
        if is_broken(X[end][structure])
            println("Fibroblast is broken.")
            break
        end
        cf = k * (X[end][focal_tesions[2]] - X[end][focal_tesions[1]] - (B - A)) # Calculate current contractile force between focal tesions
        WRITESIM && push!(df, (t, round.(X[end][filaments],digits=6), round.(X[end][motors],digits=6), round.(X[end][focal_tesions],digits=6), Y, round(cf,digits=6), isempty(Y[1]) ? "Motor dettached" : "Motor attached"))
        push!(contractile_force, cf)
        PLOTSIM && plot_sim(X[end],Y,t)
        od = OnceDifferentiable(x -> E(x, X[end], O(X[end][filaments]), Oᵩ(X[end][filaments], X[end][focal_tesions]), Y), X[end]; autodiff=:forward)
        next_x = Optim.minimizer(optimize(od, X[end], LBFGS()))
        push!(X, next_x)
        Y = update_af(X[end],Y)
    end
    # Output results
    PLOTSIM && gif(anim, "$(filename).gif", fps=5)
    WRITESIM && CSV.write("$(filename).csv", df)
    return contractile_force
end;

# anim = Animation()
# con = main(true,true,"write-csv");
# plot(1:T,con,legend=false,title="")
# savefig("system-testing-8iii-motor.png")

conf = []
for i in 1:5
    println("Run $i")
    global anim = Animation()
    con = main(true,true,"single-motor-run-a-$(i)")
    push!(conf,con)
end
plot(1:T,conf,title="Single motor contractile forces")
savefig("single-motor-contractile-force-a.png")
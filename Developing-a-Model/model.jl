using Plots, Optim, DataFrames, CSV
using Distributions: Poisson, sample
using StatsBase: mean
using LinearAlgebra: eigvals

# Simulation parameters
const A = 0.0 # Left end point
const B = 1.1 # Right end point
const Δt = 0.1 # Step size
const T = 1000 # Number of time steps

# Filament parameters
const N = 2 # Number of actin filaments, use odd number to ensure they do not visually overlap with focal tesions
const filaments = 1:N # Indexes of filaments
const L = 1.0 # Length (μm) of filaments
# P = rand((-1,1),N) # Polarities of filaments, left to right 1 represents -ve to +ve
P = [-1,1]

# Motor parameters
const M = 1 # Maximum number of myosin motors
const motors = N+1:N+M # Indexes of motors
const λ = round(Int,0.75*M) # Average number of motors currently attached
const β = 0.05 # Probability motor detaches per second

# Focal tesion parameters
const k = 1.0 # Stiffness between focal tesions
const focal_tesions = N+M+1:N+M+2 # Indexes of focal tesions
const structure = [filaments; focal_tesions] # Indexes of all filaments and focal tesions

# Compute overlap between two fibres x1 and x2
overlap(x1,x2) = max(L - abs(x1-x2), 0.0)
# Symmetric matrix where O[i,j] is half the overlap between filament i and j and diagonal = 0 (filament can't overlap with self)
O(x) = [x1 == x2 ? 0.0 : overlap(x1,x2) for x1 in x, x2 in x]
# An Nx2 matrix when O[i,j] gives the overlap between filament i and focal tesion j
Oᵩ(x,f) = [overlap(xi,fj) for xi in x, fj in f]
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

function plot_sim(x,y,t)
    padding = 2*L
    X_span = A-padding:B+padding # Range for filament movement
    # Set up domain
    plot(xlims=(first(X_span), last(X_span)), ylims=(0,N+1), title="Time = $t")
    # Display filaments
    a = x[filaments] .- L/2 # Start of filaments
    b = x[filaments] .+ L/2 # End of filaments
    line_colour = reshape([p > 0 ? :blue : :red for p in P], 1, length(P)) # Polarity of filaments
    plot!([a b]', [filaments filaments]', legend=false, lc=line_colour)
    # Display motors and their attachments
    for m = motors
        af = y[m-N]
        if isempty(af)
            scatter!([x[m]], [0], markercolor=:black, markerstroke=0, markersize=3)
        else
            scatter!([x[m]], [mean(af)], markercolor=:black, markerstroke=0, markersize=3)
            plot!([x[m];x[m]], [af[1];af[end]], legend=false, lc=:black)
        end
    end
    # Display focal tesions
    N % 2 == 0 ? height = N/2+L/2 : height = N/2 # Ensure focal tesions are not overlapping a filament
    plot!([x[focal_tesions[1]]-L/2 x[focal_tesions[2]]-L/2;x[focal_tesions[1]]+L/2 x[focal_tesions[2]]+L/2], [height height; height height], legend=false, lc=:black, linewidth=2)
    plot!(show=true)
    frame(anim)
end

# Energy functional where 
    # x...the next iterate (minimiser)
    # xn...the previous iterate
    # O_mat...the matrix storing the overlap between filaments
    # Oᵩ_mat... the matrix storing the overlap between filaments and focal tesions
    # y...y[m] stores the list of filaments motor m is attached to
function E(x, xn, O_mat, Oᵩ_mat, y)
    
    # Parameters
    F = 0.0#000001 # Pulling factor
    ξ = 0.0000001 # Drag on filaments due to moving through the cytoplasm (think about units!)
    η = 15.0 # Effective viscous drag (pNs/μm^2) due to cross-linker proteins
    Fs = 5.0 # Motor stall force (pN)
    Vm = 0.5 # Maximum (load-free) motor velocity (μm/s)
    ρ = 1.0 # Drag on filaments due to climbing through extracellular matrix (think about units!)
    
    res = -x[1] * F 
    
    # Filament movement
    for i in filaments
        res += ξ * (x[i]-xn[i])^2/(2*Δt)
        for j in filaments
            res += η * O_mat[i,j] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*Δt) # Cross linkers
        end
        for j in focal_tesions
            res += ρ * Oᵩ_mat[i,j-N-M] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*Δt) # Focal tesions
        end
    end
    
    # Motor movement
    for i in motors
        af = y[i-N] # Indices of attached filaments to current motor
        for j in 1:length(af)
            res += Fs/Vm * (x[af[j]] - x[i] - (xn[af[j]] - xn[i]))^2/(2*Δt) + Fs * (x[af[j]] - x[i]) * P[af[j]]
        end
    end

    # Focal tesion spring attachment
    res += k/2 * (x[focal_tesions[2]] - x[focal_tesions[1]] - (B - A))^2

    return res
end

mf_overlap(xm,x1) = x1 - L/2 < xm < x1 + L/2 # Check if motor is still touching a filament
attached(xm,x1,x2) = x1 - L/2 < xm < x1 + L/2 && x2 - L/2 < xm < x2 + L/2 # Check two filaments are attached to a motor

function gen_af(m,x)
    all_af = [n for n in filaments if mf_overlap(x[N+m],x[n])]
    length(all_af) ≥ 2 ? (return sort!(sample(all_af,2,replace=false))) : (return [])
end

# Update the filaments attached to each fibre, x stores position of N filaments and M motors, y stores the connections
function update_af(x,y)
    # Remove motor-filament attachments if they are no longer overlapping or randomly
    for m in 1:M
        if !isempty(y[m]) && (!attached(x[N+m],x[y[m][1]],x[y[m][end]]))# || rand() < β*Δt)
            y[m] = []
        end
    end

    m_status = (!isempty).(y) # m_status[m]... 1 if motor m is currently attached to filaments, 0 otherwise
    attached_count = sum(m_status)
    unattached_count = M - attached_count
    
    # Add motor-filament attachment randomly
    desired_attached_count = min(M,rand(Poisson(λ)))
    to_attach_count = min(max(0,desired_attached_count-attached_count),unattached_count)
    motors_to_attach = sample(findall(==(0),m_status), to_attach_count)
    for m in motors_to_attach
        y[m] = gen_af(m,x)
    end

    return y # Updated list of motor-filament connections
end

function main(PLOTSIM,WRITESIM,filename)
    df = DataFrame(Time=Int64[], FilamentPos=Vector{Vector{Float64}}(), MotorPos=Vector{Vector{Float64}}(), FocalTesionPos=Vector{Vector{Float64}}(), MotorConnections=Vector{Vector{Vector{Int64}}}(), ContractileForce=Float64[], Other=String[])
    X = [vcat(B .*rand(N), B .*rand(M),A,B)] # Centre point of N filaments, M motors and focal tesions centred at end points [A,B]
    # X = [[0.49520812591209784, 0.9129525143013779, 0.5, 0.0, 1.1]]
    # Ensure no initial breakage
    # count = 0
    # while is_broken(X[end][structure])
    #     count > 10 && error("Could not find valid initial structure. Try adding more filaments, or shrinking the interval.")
    #     count += 1
    #     X = [vcat(B .*rand(N), B .*rand(M),A,B)]
    # end
    # println("Accepted initial structure")
    # display(adj(O(X[end][structure])))
    # display(deg(O(X[end][structure])))
    # display(laplace(X[end][structure]))
    # println("Has fiedler value: ", eigvals(laplace(X[end][structure]))[2])
    Y = [Int64[] for m in 1:M] # Y[m]...List of fibres attached to motor m
    # Generate filament-motor connections
    for m in sample(1:M,λ,replace=false)
        Y[m] = gen_af(m, X[end])
    end
    contractile_force = [] # Contractile force between focal tesions
    # Evolve simulation
    for t in 1:T
        println("---Iteration $t/$T---")
        # if is_broken(X[end][structure])
        #     println("Fibroblast is broken.")
        #     break
        # end
        cf = k * (X[end][focal_tesions[2]] - X[end][focal_tesions[1]] - (B - A)) # Calculate current contractile force between focal tesions
        WRITESIM && push!(df, (t, round.(X[end][filaments],digits=6), round.(X[end][motors],digits=6), round.(X[end][focal_tesions],digits=6), Y, round(cf,digits=6), isempty(Y[1]) ? "Motor dettached" : "Motor attached"))
        push!(contractile_force, cf)
        PLOTSIM && plot_sim(X[end],Y,t)
        od = OnceDifferentiable(x -> E(x, X[end], O(X[end][filaments]), Oᵩ(X[end][filaments],X[end][focal_tesions]), Y), X[end]; autodiff=:forward)
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

## TODO
# Variation in parameters affecting contractile force, what causes the system to break
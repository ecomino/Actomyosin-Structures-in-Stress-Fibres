using Plots, Optim
using Distributions: Poisson, sample
using StatsBase: mean
# using ForwardDiff

# Simulation parameters
const A = 0 # Left end point
const B = 5 # Right end point
const Δt = 0.1 # Step size
const T = 20 # Number of time steps

# Filament parameters
const N = 30 # Number of actin filaments
const L = 1 # Length (μm) of filaments 

# Motor parameters
const M = 15 # Maximum number of myosin motors
const λ = 0.75*M # Average number of motors currently attached
const β = 0.05 # Probability motor detaches per second
P = rand((-1,1),N) # Polarities of filaments, left to right 1 represents -ve to +ve

anim = Animation()

# Compute overlap between two fibres x1 and x2
overlap(x1,x2) = max(L - abs(x1-x2), 0)
# Symmetric matrix where O[i,j] is half the overlap between filament i and j and diagonal = 0 (filament can't overlap with self)
O(x) = [x1 == x2 ? 0 : overlap(x1,x2) for x1 in x, x2 in x]
# An Nx2 matrix when O[i,j] gives the overlap between filament i and focal tesion j
Oᵩ(x,f) = [overlap(xi,fj) for xi in x, fj in f]

function plot_sim(x,y,t)
    X_span = A-L/2:B+L/2 # Range for filament movement
    # Set up domain
    plot(xlims=(first(X_span), last(X_span)), ylims=(0,N), title="Time = $t")
    # Display filaments
    n = 1:N
    a = x[n] .- L/2 # Start of filaments
    b = x[n] .+ L/2 # End of filaments
    line_colour = reshape([p > 0 ? :blue : :red for p in P], 1, length(P)) # Polarity of filaments
    plot!([a b]', [n n]', legend=false, lc=line_colour)
    # Display motors and their attachments
    for m = N+1:N+M
        cf = y[m-N]
        if !isempty(cf)
            scatter!([x[m]], [mean(cf)], markercolor=:black, markerstroke=0, markersize=2.5)
            plot!([x[m];x[m]], [cf[1];cf[end]], legend=false, lc=:black)
        end
    end
    # Display focal tesions
    plot!([x[end-1]-L/2;x[end-1]+L/2], [N/2;N/2], legend=false, lc=:black, linewidth=2)
    plot!([x[end]-L/2;x[end]+L/2], [N/2;N/2], legend=false, lc=:black, linewidth=2)
    plot!(show=true)
    frame(anim)
end

# Energy function where 
    # x...the next iterate (minimiser)
    # xn...the previous iterate
    # O_mat...the matrix storing the overlap between filaments
    # y...y[m] stores the list of filaments motor m is connected to
function E(x, xn, O_mat, Oᵩ_mat, y)
    
    # Parameters
    F = 0.2 # Pulling factor
    ξ = 1 # Drag on filaments due to moving through the extracellular matrix (think about units!)
    η = 15 # Effective viscous drag (pNs/μm^2) due to cross-linker proteins
    Fs = 5 # Motor stall force (pN)
    Vm = 0.5 # Maximum (load-free) motor velocity (μm/s)
    ρ = 0.001 # Drag on filaments due to climbing through extracellular matrix (think about units!)
    
    res = ξ * sum((x-xn).^2/(2*Δt)) - x[1] * F 
    
    # Filament movement
    for i in 1:N
        for j in 1:N
            res += η * O_mat[i,j] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*Δt) # Cross linkers
        end
        for j in N+M+1:N+M+2
            res += ρ * Oᵩ_mat[i,j-N-M] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*Δt)
        end
    end
    
    # Motor movement
    for i in N+1:N+M
        cf = y[i-N] # Indices of connected filaments to current motor
        for j in 1:length(cf)
            res += Fs/Vm * (x[cf[j]] - x[i] - (xn[cf[j]] - xn[i]))^2/(2*Δt) + Fs * (x[cf[j]] - x[i]) * P[cf[j]]
        end
    end
    return res
end

mf_overlap(xm,x1) = x1 - L/2 < xm < x1 + L/2
connected(xm,x1,x2) = x1 - L/2 < xm < x1 + L/2 && x2 - L/2 < xm < x2 + L/2

function gen_cf(m,x)
    all_cf = [n for n in 1:N if mf_overlap(x[N+m],x[n])]
    length(all_cf) ≥ 2 ? (return sort!(sample(all_cf,2,replace=false))) : (return [])
end

# Update the filaments connected to each fibre, x stores position of N filaments and M motors, y stores the connections
function update_cf(x,y)
    # Remove motor-filament attachments if they are no longer overlapping or randomly
    for m in 1:M
        if !isempty(y[m]) && (!connected(x[N+m],x[y[m][1]],x[y[m][end]]) || rand() < β*Δt)
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
        y[m] = gen_cf(m,x)
    end

    return y # Updated list of motor-filament connections
end

function main()
    X = [vcat(B .*rand(N), B .*rand(M),A,B)] # Centre point of N filaments, M motors and focal tesions centred at end points [A,B]
    Y = [[] for m in 1:M] # Y[m]...List of fibres connected to motor m
    if O(X[end][1:N])
    ft_pos = [] # Will store position of focal tesions
    # Generate filament-motor connections
    for m in sample(1:M,round(Int,λ),replace=false)
        Y[m] = gen_cf(m, X[end])
    end
    # Evolve simulation
    for t in 1:T
        println("Iteration $t/$T")
        plot_sim(X[end],Y,t)
        next_x = optimize(x -> E(x, X[end], O(X[end][1:N]), Oᵩ(X[end][1:N],X[end][end-1:end]), Y), X[end]) # , LBFGS(); autodiff = :forward
        push!(X, Optim.minimizer(next_x))
        Y = update_cf(X[end],Y)
        push!(ft_pos, X[end][end-1:end])
    end
    # Output results
    gif(anim, "bottleneck-testing-2.gif", fps=5)
end

main()

# Tam supplementary material, Oelz supp for crosslinker drag (look at units!)
# Add spring and calculate force at each time step
# Variation in parameters affecting contractile force, what causes the system to break
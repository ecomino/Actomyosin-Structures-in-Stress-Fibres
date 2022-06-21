using Plots, Optim
using Distributions: Poisson, sample
using StatsBase: mean

# Simulation parameters
const A = 0 # Left end point
const B = 5 # Right end point
const Δt = 0.1 # Step size
const T = 500 # Number of time steps

# Filament parameters
const N = 30 # Number of actin filaments
const L = 1 # Length (μm) of filaments
const k = 1 # Stiffness between focal tesions

# Motor parameters
const M = 20 # Maximum number of myosin motors
const λ = round(Int,0.75*M) # Average number of motors currently attached
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
        af = y[m-N]
        if !isempty(af)
            scatter!([x[m]], [mean(af)], markercolor=:black, markerstroke=0, markersize=2.5)
            plot!([x[m];x[m]], [af[1];af[end]], legend=false, lc=:black)
        end
    end
    # Display focal tesions
    plot!([x[end-1]-L/2 x[end]-L/2;x[end-1]+L/2 x[end]+L/2], [N/2 N/2; N/2 N/2], legend=false, lc=:black, linewidth=2)
    plot!(show=true)
    frame(anim)
end

# Energy functional where 
    # x...the next iterate (minimiser)
    # xn...the previous iterate
    # O_mat...the matrix storing the overlap between filaments
    # y...y[m] stores the list of filaments motor m is attached to
function E(x, xn, O_mat, Oᵩ_mat, y)
    
    # Parameters
    F = 0 # Pulling factor
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
            res += ρ * Oᵩ_mat[i,j-N-M] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*Δt) # Focal tesions
        end
    end
    
    # Motor movement
    for i in N+1:N+M
        af = y[i-N] # Indices of attached filaments to current motor
        for j in 1:length(af)
            res += Fs/Vm * (x[af[j]] - x[i] - (xn[af[j]] - xn[i]))^2/(2*Δt) + Fs * (x[af[j]] - x[i]) * P[af[j]]
        end
    end

    # Focal tesion spring attachment
    res += k/2 * (x[end] - x[end-1] - (B - A))^2

    return res
end

mf_overlap(xm,x1) = x1 - L/2 < xm < x1 + L/2
attached(xm,x1,x2) = x1 - L/2 < xm < x1 + L/2 && x2 - L/2 < xm < x2 + L/2

function gen_af(m,x)
    all_af = [n for n in 1:N if mf_overlap(x[N+m],x[n])]
    length(all_af) ≥ 2 ? (return sort!(sample(all_af,2,replace=false))) : (return [])
end

# Update the filaments attached to each fibre, x stores position of N filaments and M motors, y stores the connections
function update_af(x,y)
    # Remove motor-filament attachments if they are no longer overlapping or randomly
    for m in 1:M
        if !isempty(y[m]) && (!attached(x[N+m],x[y[m][1]],x[y[m][end]]) || rand() < β*Δt)
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

function main(PLOTSIM)
    X = [vcat(B .*rand(N), B .*rand(M),A,B)] # Centre point of N filaments, M motors and focal tesions centred at end points [A,B]
    X = [[1.869132361544753, 4.850869302868191, 1.4410882371816887, 2.637691108994159, 1.9373120508746218, 1.2844397986073726, 3.6098891422490835, 0.9324901695196425, 4.400263116057822, 1.2680141471185173, 4.4999564851316585, 4.997028360785933, 3.4835766462290008, 4.831418313854082, 0.0029281889183413456, 2.7498062620647135, 1.3026631775391684, 3.552459897773908, 0.04245513213551433, 3.1048005516440873, 2.275424155589181, 3.5017656090150995, 1.929672332565064, 2.664508515378249, 1.9427542606844321, 4.261940895754634, 4.14533574754058, 0.23608459303215756, 2.736539059089296, 2.2175589771404285, 2.355319374169276, 4.94425640709104, 3.19318908374502, 2.902125370606216, 0.2036651440595877, 4.870194704825631, 3.3476663531110424, 4.433407715457577, 0.37570545809965106, 2.2019235508726256, 2.4669365687838813, 0.6586563429464387, 2.616100124626551, 4.101664174622075, 0.8286650023859315, 1.14617907313739, 1.3420598719032322, 1.9580717359505955, 3.5408066436186663, 2.2397129219679073, 0.0, 5.0]]
    # Ensure no initial breakage
    # while 0 ∈ sum(O(X[end]),dims=1)
    #     X = [vcat(B .*rand(N), B .*rand(M),A,B)]
    # end
    Y = [[] for m in 1:M] # Y[m]...List of fibres attached to motor m
    # Generate filament-motor connections
    for m in sample(1:M,λ,replace=false)
        Y[m] = gen_af(m, X[end])
    end
    Y = [[], [2, 11], [13, 20], [4, 24], [15, 19], [2, 14], [13, 20], [2, 11], [19, 28], [25, 30], [4, 29], [8, 28], [], [], [6, 10], [3, 10], [3, 17], [5, 30], [], []]
    contractile_force = [] # Contractile force between focal tesions
    # Evolve simulation
    for t in 1:T
        println("Iteration $t/$T")
        cf = k * (X[end][end] - X[end][end-1] - (B - A)) # Calculate current contractile force between focal tesions
        push!(contractile_force, cf)
        PLOTSIM && plot_sim(X[end],Y,t)
        # next_x = Optim.minimizer(optimize(x -> E(x, X[end], O(X[end][1:N]), Oᵩ(X[end][1:N],X[end][end-1:end]), Y), X[end]))
        od = OnceDifferentiable(x -> E(x, X[end], O(X[end][1:N]), Oᵩ(X[end][1:N],X[end][end-1:end]), Y), X[end]; autodiff = :forward)
        next_x = Optim.minimizer(optimize(od, X[end], LBFGS()))
        push!(X, next_x)
        Y = update_af(X[end],Y)
    end
    # Output results
    PLOTSIM && gif(anim, "same-ic-autodiff.gif", fps=5)
    return contractile_force
end
conf = []
for i in 1:3
    println("Run $i")
    con = main(false)
    push!(conf,con)
end
plot(1:T,conf,legend=false,title="Autodiff LBFGS")
savefig("cf-500-autodiff.png")

## TODO
# Variation in parameters affecting contractile force, what causes the system to break
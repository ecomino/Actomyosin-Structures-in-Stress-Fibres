using Plots, Optim, DataFrames, CSV
using Distributions: Poisson, sample
using StatsBase: mean

# Simulation parameters
const A = 0.0 # Left end point
const B = 1.1 # Right end point
const Δt = 0.1 # Step size
const T = 20 # Number of time steps

# Filament parameters
const N = 2 # Number of actin filaments, use odd number to ensure they do not visually overlap with focal adhesions
const filaments = 1:N # Indexes of filaments
const L = 1.0 # Length (μm) of filaments
P = rand((-1,1),N) # Polarities of filaments, left to right 1 represents -ve to +ve

# Motor parameters
const M = 1 # Maximum number of myosin motors
const motors = N+1:N+M # Indexes of motors
const λ = floor(Int,0.75*M) # Average number of motors currently attached
const β = 0.05 # Probability motor detaches per second

# focal adhesion parameters
const focal_adhesions = N+M+1:N+M+2 # Indexes of focal adhesions
const structure = [filaments; focal_adhesions] # Indexes of all filaments and focal adhesions

# Energy functional parameters
const F = 0.0 #000001 # Pulling factor
const ξ = 0.0000001 # Drag on filaments due to moving through the cytoplasm (think about units!)
const η = 15.0 # Effective viscous drag (pNs/μm^2) due to cross-linker proteins
const Fs = 5.0 # Motor stall force (pN)
const Vm = 0.5 # Maximum (load-free) motor velocity (μm/s)
const ρ = 10.0 # Drag on filaments due to climbing through extracellular matrix (think about units!)
const k = 100.0 # Stiffness between focal adhesions

# Compute overlap between two fibres x1 and x2
overlap(x1,x2) = max(L - abs(x1-x2), 0.0)
# Symmetric matrix where O[i,j] is half the overlap between filament i and j and diagonal = 0 (filament can't overlap with self)
O(x) = [x1 == x2 ? 0.0 : overlap(x1,x2) for x1 in x, x2 in x]
# An Nx2 matrix when O[i,j] gives the overlap between filament i and focal adhesion j
Oᵩ(x,f) = [overlap(xi,fj) for xi in x, fj in f]

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
    # Display focal adhesions
    N % 2 == 0 ? height = N/2+L/2 : height = N/2 # Ensure focal adhesions are not overlapping a filament
    plot!([x[focal_adhesions[1]]-L/2 x[focal_adhesions[2]]-L/2;x[focal_adhesions[1]]+L/2 x[focal_adhesions[2]]+L/2], [height height; height height], legend=false, lc=:black, linewidth=2)
    plot!(show=true)
    frame(anim)
end

# Energy functional where 
    # x...the next iterate (minimiser)
    # xn...the previous iterate
    # O_mat...the matrix storing the overlap between filaments
    # Oᵩ_mat... the matrix storing the overlap between filaments and focal adhesions
    # y...y[m] stores the list of filaments motor m is attached to
function E(x, xn, O_mat, Oᵩ_mat, y)    
    
    res = -x[1] * F 
    
    # Filament movement
    for i in filaments
        res += ξ * (x[i]-xn[i])^2/(2*Δt) # Filament drag through cytoplasm
        for j in filaments
            res += η * O_mat[i,j] / 2 * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*Δt) # Drag effect from cross linker proteins
        end
        for j in focal_adhesions
            res += ρ * Oᵩ_mat[i,j-N-M] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*Δt) # Drag due to "fixed" focal adhesions
        end
    end
    
    # Motor movement
    for i in motors
        af = y[i-N] # Indices of attached filaments to current motor
        for j in 1:length(af)
            res += Fs * (x[af[j]] - x[i]) * P[af[j]] - Fs/Vm * (x[af[j]] - x[i] - (xn[af[j]] - xn[i]))^2/(2*Δt)
        end
    end

    # Focal adhesion spring attachment
    res += k/2 * (x[focal_adhesions[2]] - x[focal_adhesions[1]] - (B - A))^2

    return res
end

function calculate_motor_force(af,m,x,xn)
    res = 0.0
    for j in 1:length(af)
        res += Fs/Vm * (x[af[j]] - x[N+m] - (xn[af[j]] - xn[N+m]))/Δt - Fs * P[af[j]] 
    end
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

function serialize_motor_connections(Y)
    output = "[ "
    for y in Y
        if isempty(y)
            output *= "[] "
        else
            output *= "[$(y[1]),$(y[2])] "
        end
    end
    return output * "]"
end

function main(PLOTSIM,WRITESIM,filename)
    df = DataFrame(Time=Int64[], FilamentPos=Vector{Vector{Float64}}(), MotorPos=Vector{Vector{Float64}}(), FocalTesionPos=Vector{Vector{Float64}}(), MotorConnections=String[], ContractileForce=Float64[], Other=String[])
    X = [vcat(B .* rand(N), B .* rand(M), A, B)] # Centre point of N filaments, M motors and focal adhesions centred at end points [A,B]
    X = [[0.8, 0.4, 0.6, A, B]]
    Y = [[] for m in 1:M] # Y[m]...List of fibres attached to motor m
    # Generate filament-motor connections
    for m in sample(1:M,λ,replace=false)
        Y[m] = gen_af(m, X[end])
    end
    Y = [[1,2]]
    contractile_force = [] # Contractile force between focal adhesions
    mf = 0
    # Evolve simulation
    for t in 1:T
        println("Iteration $t/$T")
        cf = k * (X[end][focal_adhesions[2]] - X[end][focal_adhesions[1]] - (B - A)) # Calculate current contractile force between focal adhesions
        if t > 1
            @show Fs/Vm * ((X[end][2]+X[end][1]-(X[end-1][2]+X[end-1][1])-2*(X[end][3]-X[end-1][3]))/Δt) 
            #mf = ((X[end][2]-X[end][1]-(X[end-1][2]-X[end-1][1]))/Δt) * (Fs/Vm + 2*η*O(X[end][filaments])[1,2]) + 2*Fs # Calculate motor force of single motor
        end
        (WRITESIM && t!=1) && push!(df, (t, round.(X[end][filaments],digits=6), round.(X[end][motors],digits=6), round.(X[end][focal_adhesions],digits=6), serialize_motor_connections(Y), round(cf,sigdigits=6), isempty(Y[1]) ? "Motor dettached / Motor force $(mf)" : "Motor attached / Motor force $(mf)"))
        push!(contractile_force, cf)
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
# con = main(true,false,"pos");
# plot(1:T,con,legend=false,title="")
# savefig("params-rho-1-cf.png")

conf = []
for i in 1:1
    println("---Run $i---")
    global anim = Animation()
    global P = [1,-1] #rand((-1,1),N)
    con = main(true,true,"motor-force-1")
    push!(conf,con)
end
plot(1:T,conf,title="Contractile Forces",xlabel="Time",ylabel="Contractile Force")
# savefig("updated-motor-force-cf.png")
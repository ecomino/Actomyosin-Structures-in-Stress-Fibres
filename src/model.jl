using Distributions: Bernoulli
include("model-parameters.jl")

# Energy functional where 
    # x...the next iterate (minimiser)
    # xn...the previous iterate
    # O_mat...the matrix storing the overlap between filaments
    # Oᵩ_mat... the matrix storing the overlap between filaments and focal adhesions
    # y...y[m] stores the list of filaments motor m is attached to
function E(x, xn, O_mat, Oᵩ_mat, y, p)    
    
    res = 0.0
    
    # Filament movement
    for i in p.filaments
        res += p.ξ * (x[i]-xn[i])^2/(2*p.Δt) # Filament drag through cytoplasm
        for j in p.filaments
            res += p.η * O_mat[i,j] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*p.Δt) # Drag effect from cross linker proteins
        end
        for j in p.focal_adhesions
            res += p.ρ * Oᵩ_mat[i,j-p.N-p.M] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*p.Δt) # Drag due to focal adhesions
        end
    end
    
    # Motor movement
    for i in p.motors
        af = y[i-p.N] # Indices of attached filaments to current motor
        for j in af
            res += - p.Fs * (x[j] - x[i]) * p.P[j] + p.Fs/p.Vm * (x[j] - x[i] - (xn[j] - xn[i]))^2/(2*p.Δt)
        end
    end

    # Focal adhesion spring attachment
    res += p.k/2 * (x[p.focal_adhesions[2]] - x[p.focal_adhesions[1]] - (p.B - p.A))^2

    return res
end

function reset_model(p)
    X = Vector{Vector{Float64}}(undef, p.T+1)
    X[1] = vcat(p.B .* rand(p.N), p.B .* rand(p.M), p.A, p.B) # Initial centre point of N filaments, M motors and focal adhesions centred at end points [A,B]
    p.filaments = 1:p.N
    p.motors = p.N+1:p.N+p.M
    p.focal_adhesions = p.N+p.M+1:p.N+p.M+2
    p.structure = [p.filaments; p.focal_adhesions]
    p.λ = floor(Int,0.75*p.M)
    @assert p.δ ≤ p.L && p.δ ≥ 0.0
    Y = [[] for _ in 1:p.M] # Y[m]...List of fibres attached to motor m
    # Generate filament-motor connections
    for m in sample(1:p.M,p.λ,replace=false)
        Y[m] = gen_af(m, X[1], p.N, p.filaments, p.P, p.L)
    end
    cf = Vector{Float64}(undef,p.T+1) # Contractile force between focal adhesions at each time step
    p.P = rand((-1,1),p.N)
    return X, Y, cf, p
end

# Symmetric matrix where O[i,j] is the overlap between filament i and j and diagonal = 0 (filament can't overlap with self)
O(x,L) = [x1 == x2 ? 0.0 : max(L - abs(x1-x2), 0.0) for x1 in x, x2 in x]
# An Nx2 matrix when O[i,j] gives the overlap between filament i and focal adhesion j
Oᵩ(x,f,l) = [max(l - abs(xi-fj), 0.0) for xi in x, fj in f]

mf_overlap(xm,x1,L) = x1 - L/2 < xm < x1 + L/2 # Check if motor is touching a filament
ma_interval(xm,x1,p,L,δ) = p*(xm-x1) .+ [L/2 - δ, - L/2] # Distance motor is from drop off point and filament's opposite end

# Generate attached filaments
function gen_af(m,x,N,filaments,P,L)
    all_af = [n for n in filaments if mf_overlap(x[N+m],x[n],L)]
    sort!(all_af)
    if length(all_af) ≥ 2
        af = sample(all_af,2,replace=false)
        if P[af[1]] != P[af[2]]
            return af
        end
    end
    return []
end

# Update the filaments attached to each fibre, x stores position of N filaments and M motors, y stores the connections
function update_af(x,y,params)
    for m in 1:params.M
        if !isempty(y[m]) # Check motor is attached to filaments
            for i in y[m] # Indexes of filaments connected to motor m
                K, R = ma_interval(x[params.N+m],x[i],params.P[i],params.L,params.δ) # Interval that motor should be in to stay attached to filament
                # Remove motor-filament attachments if the motor has passed the dropoff point, detached, or randomly
                if (K < 0) || (R > 0)
                    y[m] = []
                end
            end
        end
    end

    m_status = (!isempty).(y) # m_status[m]... 1 if motor m is currently attached to filaments, 0 otherwise
    attached_count = sum(m_status)
    unattached_count = params.M - attached_count
    
    # Add motor-filament attachment randomly
    desired_attached_count = min(params.M,rand(Poisson(params.λ)))
    to_attach_count = min(max(0,desired_attached_count-attached_count),unattached_count)
    motors_to_attach = sample(findall(==(0),m_status), to_attach_count)
    for m in motors_to_attach
        y[m] = gen_af(m,x,params.N,params.filaments,params.P,params.L)
    end

    return y # Updated list of motor-filament connections
end

function turnover(x,params)
    for i in params.filaments
        rand(Bernoulli(params.α)) && (x[i] = params.B * rand())
    end
    for j in params.motors
        rand(Bernoulli(params.β)) && (x[j] = params.B * rand())
    end
    return x
end;


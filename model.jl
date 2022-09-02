using Distributions: Bernoulli
include("model-parameters.jl")

# Energy functional where 
    # x...the next iterate (minimiser)
    # xn...the previous iterate
    # O_mat...the matrix storing the overlap between filaments
    # Oᵩ_mat... the matrix storing the overlap between filaments and focal adhesions
    # y...y[m] stores the list of filaments motor m is attached to
function E(x, xn, O_mat, Oᵩ_mat, y)    
    
    res = 0.0
    
    # Filament movement
    for i in filaments
        res += ξ * (x[i]-xn[i])^2/(2*Δt) # Filament drag through cytoplasm
        for j in filaments
            res += η * O_mat[i,j] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*Δt) # Drag effect from cross linker proteins
        end
        for j in focal_adhesions
            res += ρ * Oᵩ_mat[i,j-N-M] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*Δt) # Drag due to "fixed" focal adhesions
        end
    end
    
    # Motor movement
    for i in motors
        af = y[i-N] # Indices of attached filaments to current motor
        for j in af
            res += - Fs * (x[j] - x[i]) * P[j] + Fs/Vm * (x[j] - x[i] - (xn[j] - xn[i]))^2/(2*Δt)
        end
    end

    # Focal adhesion spring attachment
    res += k/2 * (x[focal_adhesions[2]] - x[focal_adhesions[1]] - (B - A))^2

    return res
end

function reset_model()
    X = Vector{Vector{Float64}}(undef, T+1)
    X[1] = vcat(B .* rand(N), B .* rand(M), A, B) # Initial centre point of N filaments, M motors and focal adhesions centred at end points [A,B]
    global filaments = 1:N
    global motors = N+1:N+M
    global focal_adhesions = N+M+1:N+M+2
    global structure = [filaments; focal_adhesions]
    global λ = floor(Int,0.75*M)
    # X = [[0.25, 0.75, 0.5, A, B]]
    Y = [[] for m in 1:M] # Y[m]...List of fibres attached to motor m
    # Generate filament-motor connections
    for m in sample(1:M,λ,replace=false)
        Y[m] = gen_af(m, X[1])
    end
    # Y = [[1,2]]
    cf = Vector{Float64}(undef,T+1) # Contractile force between focal adhesions at each time step
    global P = rand((-1,1),N)
    return X, Y, cf
end

# Symmetric matrix where O[i,j] is the overlap between filament i and j and diagonal = 0 (filament can't overlap with self)
O(x) = [x1 == x2 ? 0.0 : max(L - abs(x1-x2), 0.0) for x1 in x, x2 in x]
# An Nx2 matrix when O[i,j] gives the overlap between filament i and focal adhesion j
Oᵩ(x,f) = [max(l - abs(xi-fj), 0.0) for xi in x, fj in f]

mf_overlap(xm,x1) = x1 - L/2 < xm < x1 + L/2 # Check if motor is touching a filament
ma_interval(xm,x1,p) = p*(xm-x1) .+ [L/2 - υ, -L/2] # Distance motor is from drop off point and filament's opposite end

# Generate attached filaments
function gen_af(m,x)
    all_af = [n for n in filaments if mf_overlap(x[N+m],x[n])]
    length(all_af) ≥ 2 ? (return sort!(sample(all_af,2,replace=false))) : (return [])
end

# Update the filaments attached to each fibre, x stores position of N filaments and M motors, y stores the connections
function update_af(x,y)
    for m in 1:M
        if !isempty(y[m]) # Check motor is attached to filaments
            for i in y[m] # Indexes of filaments connected to motor m
                K, R = ma_interval(x[N+m],x[i],P[i]) # Interval that motor should be in to stay attached to filament
                # Remove motor-filament attachments if the motor has passed the dropoff point, detached, or randomly
                if (K < 0) || (R > 0) || rand(Bernoulli(β))
                    y[m] = []
                end
            end
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

function filament_turnover(x)
    for i in 1:N
        rand(Bernoulli(α)) && (x[i] = B * rand())
    end
    return x
end


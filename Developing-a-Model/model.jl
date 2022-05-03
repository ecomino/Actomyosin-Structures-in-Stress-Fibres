using Plots, Optim, LinearAlgebra, StatsBase

const N = 2 # Number of actin filaments
const M = 1 # Number of motors
const L = 1 # Length of filaments
const Δt = 0.1 # Step size
#const P = rand((-1,1),N) # Polarities of filaments, left to right 1 represents -ve to +ve
const P = [-1, 1]
anim = Animation()

#=struct Parameters
    F::Float64 # Pulling factor
    ξ::Float64 # Coefficient of drag friction
    η::Float64 # Coefficient for cross-linker proteins drag
end=#

#=
struct Motor
    f1::Int64 # Index of first connected fibre
    f2::Int64 # Index of other connected fibre
end

# Given a motor, return tuple of indexes of connected filaments 
function get_filaments(m::Motor)
    return (m.f1, m.f2)
end
=#

# Compute overlap between two fibres x1 and x2
overlap(x1,x2) = max(L - abs(x1-x2), 0)
# Symmetric Matrix where A[i,j] is half the overlap between filament i and j and diag = 0
A(x) = [x1 == x2 ? 0 : overlap(x1,x2) for x1 in x, x2 in x] 

function plot_sim(x,y,t)
    X_span = -2:10 # Range for filament movement
    plot(xlims=(first(X_span), last(X_span)), ylims=(0,N), title="Time = $t", show=false)
    for n = 1:N
        a = x[n]-L/2 # Start of filament
        b = x[n]+L/2 # End of filament
        line_colour(n) = P[n] > 0 ? :blue : :red
        plot!([a;b], [n;n], legend=false, lc=line_colour(n), xlims=(first(X_span), last(X_span)), ylims=(0,N), show=false)
    end
    for m = N+1:N+M
        cf = y[m-M-1]
        motor_pos(m) = mean(cf)
        scatter!([x[m]], [mean(cf)], markercolor=:black, markerstroke=0, markersize=2.5, show=false)
        plot!([x[m];x[m]], [cf[1];cf[end]], legend=false, lc=:black, show=false)
    end
    plot!(show=true)
    frame(anim)
end

# Energy function where 
    # x...the next iterate (minimiser)
    # xn...is the previous
    # A_mat...is the matrix storing the overlap between filaments
function E(x, xn, A_mat, y)
    
    # Parameters
    F = 0.2 # Pulling factor
    ξ = 1 # Coefficient of drag friction
    η = 0.2 # Coefficient for cross-linker proteins drag
    Fs = 0.2 # Motor stall force
    Vm = 1 # Maximum motor velocity
    
    res = ξ * sum((x-xn).^2/2*Δt) - x[1] * F
    
    # Filament movement
    for i in 1:N
        for j in 1:N
            res += η * A_mat[i,j] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*Δt) 
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

function gen_cf(m,x)
    cf = [n for n in 1:N if x[n]-L/2 < x[N+m] < x[n]+L/2]
    return sort!(sample(cf,2,replace=false))
end

function main()
    T = 10 # Number of time steps
    # X = [-1 .* vcat(rand(N),rand(M))] # Centre point of N filaments and M motors
    X = [[-0.5, 0, -0.25]]
    Y = [gen_cf(m, X[end]) for m in 1:M] # Y[i]...List of fibres connected to motor i
    for t in 1:T-1
        plot_sim(vcat(X[end],Y[end]),Y,t)
        next_x = optimize(x -> E(x, X[end], A(X[end][1:N]), Y), X[end])
        push!(X, Optim.minimizer(next_x))
    end
    gif(anim, "fibres-motors-2.gif", fps=5)
end

main()
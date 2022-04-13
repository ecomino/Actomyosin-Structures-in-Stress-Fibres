using Plots, Optim, LinearAlgebra

const N = 10 # Number of actin filaments
const M = 5 # Number of motors
const L = 1 # Length of filaments
const Δt = 0.1 # Step size
const P = rand((-1,1),N) # Polarities of filaments, left to right 1 represents -ve to +ve
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

function plot_sim(x)
    X_span = -2:10 # Range for filament movement
    plot(xlims=(first(X_span), last(X_span)), ylims=(0,N), show=true)
    for n = 1:N
        a = x[n]-L/2 # Start of filament
        b = x[n]+L/2 # End of filament
        line_colour(n) = P[n] > 0 ? :blue : :red
        plot!([a;b], [n;n], legend=false, lc=line_colour(n), show=false, xlims=(first(X_span), last(X_span)), ylims=(0,N))
    end
    for m = N+1:N+M
        scatter!([x[m]], [m-N+0.5], markercolor=:black, markerstroke=0, markersize=2.5, show=false, xlims=(first(X_span), last(X_span)), ylims=(0,N))
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
    F = 1 # Pulling factor
    ξ = 1 # Coefficient of drag friction
    η = 0.2 # Coefficient for cross-linker proteins drag
    Fs = 1 # Motor stall force
    Vm = 1 # Maximum motor velocity
    
    res = ξ * sum((x-xn).^2/2*Δt) - x[1] * F
    
    # Sum over filaments
    for i in 1:N
        for j in 1:N
            res += η * A_mat[i,j] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*Δt)
        end
    end
    
    # # Sum over motors
    # for i in N+1:N+M
    #     cf = y[i-N] # Indices of filaments connected to current motor
    #     for j in 1:length(cf)
    #         res += -Fs * (x[cf[j]] - x[i]) * P[cf[j]]
    #     end
    # end

    return res
end

function main()
    T = 10 # Number of time steps
    X = [-1.5 .* vcat(rand(N),rand(M))] # Centre point of N filaments and M motors
    Y = [[rand(1:N),rand(1:N)] for i in 1:M] # Y[i]...List of fibres connected to motor i
    θ = [] # θ[i]... tuple of fibres motor i is connected to
    for _ in 1:T-1
        next_x = optimize(x -> E(x, X[end], A(X[end]), Y), X[end])
        push!(X, Optim.minimizer(next_x))
        plot_sim(vcat(X[end],Y[end]))
    end

    gif(anim, "stress-fibres-2.gif", fps=5)
end


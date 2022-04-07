using Plots, Optim, LinearAlgebra

const N = 10 # Number of actin filaments
const L = 0.5 # Length of filaments
const Δt = 0.1 # Step size
anim = Animation()

# Compute overlap between two fibres x1 and x2
overlap(x1,x2) = max(L - abs(x1-x2), 0)
# Symmetric Matrix where A[i,j] is half the overlap between filament i and j and diag = 0
A(x) = [x1 == x2 ? 0 : overlap(x1,x2) for x1 in x, x2 in x] 

function plot_fibres(x)
    X_span = -10:10 # Range for filament movement
    plot(xlims=(first(X_span), last(X_span)), ylims=(0,N), show=true)
    for n = 1:N
        a = x[n]-L/2 # Start of filament
        b = x[n]+L/2 # End of filament
        plot!([a;b], [n;n], legend=false, lc=:blue, show=false, xlims=(first(X_span), last(X_span)), ylims=(0,N))
    end
    plot!(show=true)
    frame(anim)
end

# Energy function where 
    # x...the next iterate (minimiser)
    # xn...is the previous
    # A_mat...is the overlap matrix and 
    # F...the pulling factor
function E(x, xn, A_mat, F)
    res = ξ * sum((x-xn).^2/2*Δt) - x[1] * F
    for i in 1:N
        for j in 1:N
            res += η * A_mat[i,j] * (x[i]-x[j]-(xn[i]-xn[j]))^2/(2*Δt)
        end
    end    
    return res
end

function main()
    X = [-4 .* rand(N)] # Centre point of filaments
    F = 1 # Pulling factor
    T = 10 # Number of time steps
    for n in 1:T-1
        next_x = optimize(x->E(x, X[end], A(X[end]), F), X[end])
        push!(X, Optim.minimizer(next_x))
        plot_fibres(X[end])
    end

    gif(anim, "stress-fibres-4.gif", fps=5)
end


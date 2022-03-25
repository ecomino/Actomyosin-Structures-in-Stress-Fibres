using Plots, Optim, LinearAlgebra
N = 10 # Number of actin filaments
X = [rand(N)] # Centre point of filaments
L = 1 # Length of filaments
anim = Animation()
dist(x1, x2) = abs(x1-x2) 
overlap(x1,x2) = max(L - dist(x1,x2), 0)
A = [overlap(x1,x2) for x1 in X[1], x2 in X[1]] # Symmetric Matrix where A[i,j] is the overlap between filament i and j

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

F = 0.1 #* Matrix(I,3,1) # Pulling factor
Δt = 0.1 # Step size
T = 10 # Number of time steps
for n in 1:T-1
    E(x) = sum((x-X[end]).^2/2*Δt) - x[1] * F
    next_x = optimize(E, last(X))
    push!(X, Optim.minimizer(next_x))
    plot_fibres(X[end])
end

gif(anim, "test.gif", fps=5)


# Assume cross linker proteins cause drag friction -> Overlap between filaments causes drag


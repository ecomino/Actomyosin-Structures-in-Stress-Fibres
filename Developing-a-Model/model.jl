using Plots, Optim, LinearAlgebra
const N = 10 # Number of actin filaments
const L = 1 # Length of filaments
const Δt = 0.1 # Step size
anim = Animation()

dist(x1, x2) = abs(x1-x2) 
overlap(x1,x2) = max(L - dist(x1,x2), 0)

A(Xvec) = [overlap(x1,x2) for x1 in Xvec, x2 in Xvec] # Symmetric Matrix where A[i,j] is the overlap between filament i and j

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

function E(Xtest, Xold, Amat, F)
    res=sum((Xtest-Xold).^2/2*Δt) - Xtest[1] * F
    for i in 1:N
        for j in 1:N
            res += Amat[i,j]*(Xtest[i]-Xtest[j]-(Xold[i]-Xold[j]))^2/(2*Δt)
        end
    end    
    return res
end

function main()
    X = [rand(N)] # Centre point of filaments
    F = 0.1 #* Matrix(I,3,1) # Pulling factor
    T = 10 # Number of time steps
    for n in 1:T-1
        
        next_x = optimize(x->E(x, X, A(last(X)), F), last(X))
        push!(X, Optim.minimizer(next_x))
        plot_fibres(X[end])
    end

    gif(anim, "test.gif", fps=5)
end

# Assume cross linker proteins cause drag friction -> Overlap between filaments causes drag


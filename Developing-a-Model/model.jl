using Plots
N = 10 # Number of actin filaments
X = rand(N) # Centre point of filaments
L = 1 # Length of filaments
X_span = -3:4 # Range for filament movement
p = plot(xlims=(first(X_span), last(X_span)), ylims=(0,N))

for n = 1:N
    a = X[n]-L/2 # Start of filament
    b = X[n]+L/2 # End of filament
    plot!([a;b], [n;n], legend=false, lc=:blue)
end
p

# Assume cross linker proteins cause drag friction -> Overlap between filaments causes drag


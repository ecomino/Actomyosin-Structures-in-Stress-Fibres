using Optim, Plots, LaTeXStrings

Φ(u) = u^2/2 # function to minimise
dΦ(u) = -u # its derivative
u0 = 10 # initial guess

# -----Euler's method----- #
Δt = 0.1 # time step
N = 100 # max number of iterations
u = Array{Float64, 1}(undef, N)
u[1] = u0 

for n in 1:N-1
    u[n+1] = u[n] / (1 + Δt)
end

# Plot Φ(u) against u
x = collect(-u0:0.01:u0)
plot(x, Φ.(x), label=L"\Phi(u)=u^2/2", legend=:bottomright, color=:black, lw=0.25)
scatter!(u, Φ.(u), label="", markersize=3, markerstrokewidth=0, color=:blue)
title!("Finding a Local Minimum")
xlabel!("u")
ylabel!("Φ(u)")
savefig("eulers-method-1.png")

# Plot u against t
scatter(collect(1:N).*Δt, u, legend=false, markersize=3, markerstrokewidth=0, color=:blue)
title!("Finding a Local Minimum")
xlabel!("t")
ylabel!("u")
savefig("eulers-method-2.png")
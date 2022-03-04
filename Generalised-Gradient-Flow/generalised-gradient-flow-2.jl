using Optim, Plots

# -----Generalised gradient flow for a more complicated function----- #
Φ(u) = u^4 + 2*u^3
Δt = 0.1
N = 100
U = Array{Float64, 1}(undef, N)
U[1] = 1

for n in 1:N-1
    f(u) = (u-U[n])^2/(2*Δt) + Φ(u)
    next_u = optimize(f, -2.1, 1)
    U[n+1] = Optim.minimizer(next_u)
end

x = collect(-2.1:0.01:1)
plot(x, Φ.(x), label=L"\Phi(u)=u^4 + 2u^3", legend=:bottomright, color=:black, lw=0.25)
scatter!(U, Φ.(U), label="", markersize=3, markerstrokewidth=0, color=:blue)
title!("Finding a Local Minimum")
xlabel!("u")
ylabel!("Φ(u)")
savefig("generalised-gradient-flow-2-1.png")
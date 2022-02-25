using Optim, Plots

# -----Generalised gradient flow----- #
Φ(u) = u^2/2
u_1 = 10
Δt = 0.1
N = 100
U = Array{Float64, 1}(undef, N)
U[1] = u_1 

for n in 1:N-1
    u_n = U[n]
    @show f(u) = (u-u_n)^2/(2*Δt) + Φ(u)
    U[n+1] = optimize(f, U[n], U[n]+Δt)
end

x = collect(-u0:0.01:u0)
plot(x, Φ.(x), label="Φ(u)=u^2/2", legend=:bottomright, color=:black, lw=0.25)
scatter!(u, Φ.(u), label="", markersize=3, markerstrokewidth=0, color=:blue)
title!("Finding a Local Minimum")
xlabel!("u")
ylabel!("Φ(u)")
savefig("generalised-gradient-flow-1.png")
using Optim, Plots

# -----Generalised gradient flow----- #
Φ(u) = u^2/2
u_1 = 10
Δt = 0.1
N = 100
U = Array{Float64, 1}(undef, N)
U[1] = u_1 

for n in 1:N-1
    f(u) = (u-U[n])^2/(2*Δt) + Φ(u)
    next_u = optimize(f, -u_1, u_1) # assume minimum ∈ [-u_1, u_1]
    U[n+1] = Optim.minimizer(next_u)
end

x = collect(-u_1:0.01:u_1)
plot(x, Φ.(x), label="Φ(u)=u^2/2", legend=:bottomright, color=:black, lw=0.25)
scatter!(U, Φ.(U), label="", markersize=3, markerstrokewidth=0, color=:blue)
title!("Finding a Local Minimum")
xlabel!("u")
ylabel!("Φ(u)")
savefig("generalised-gradient-flow-1.png")

scatter(collect(1:N).*Δt, U, legend=false, markersize=3, markerstrokewidth=0, color=:blue)
title!("Finding a Local Minimum")
xlabel!("t")
ylabel!("u")
savefig("generalised-gradient-flow-2.png")
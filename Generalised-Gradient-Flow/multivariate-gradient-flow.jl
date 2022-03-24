using Optim, Plots, LaTeXStrings

# -----Multivariate gradient flow----- #
Φ(u) = u[1]^2 + u[2]^2
Φ(u,v) = u^2 + v^2
Δt = 0.1
N = 100
U = Vector{Vector{Float64}}(undef, N)
U[1] = [10, 16]

for n in 1:N-1
    f(u) = (u[1]-U[n][1])^2/(2*Δt) + (u[2]-U[n][2])^2/(2*Δt) + Φ(u)
    next_u = optimize(f, U[n])
    U[n+1] = Optim.minimizer(next_u)
end

x, y = collect(-20:0.1:20), collect(-20:0.1:20)
plot(x, y, Φ, st=:surface, camera=[30,20], xlabel="u", ylabel="v", zlabel="Φ(u,v)", colorbar=false, color=cgrad([:blue,:cyan]))
u_path, v_path = map(u -> u[1], U), map(u -> u[2], U)
scatter!(u_path, v_path, Φ.(u_path, v_path), label=L"\Phi(u,v)=u^2/2+v^2/2", legend=:outerbottom, markersize=2, markerstrokewidth=0, color=:black)
title!("Finding a Local Minimum")
savefig("multivariate-gradient-flow-1.png")
using Optim, Plots, LaTeXStrings

# -----Multivariate gradient flow with box constraints----- #
# Objective to minimise
Φ(u) = u[1]^2 + u[2]^2
Φ(u,v) = u^2 + v^2
# Subject to box constraints
lower = [0, 0.1]
upper = [Inf, Inf]
# Define g = ∇Φ(u)
function g!(storage, u)
    storage[1] = 2*u[1]
    storage[2] = 2*u[2]
end
Δt = 0.1
N = 100
U = Vector{Vector{Float64}}(undef, N)
U[1] = [10, 16]

for n in 1:N-1
    f(u) = (u[1]-U[n][1])^2/(2*Δt) + (u[2]-U[n][2])^2/(2*Δt) + Φ(u)
    od = OnceDifferentiable(f, g!, U[n])
    results = optimize(od, U[n], lower, upper, Fminbox{GradientDescent}())
    U[n+1] = Optim.minimizer(results)
end

x, y = collect(-20:0.1:20), collect(-20:0.1:20)
plot(x, y, Φ, st=:surface, xlabel="u", ylabel="v", zlabel="Φ(u,v)", colorbar=false, color=cgrad(:blues, rev=true))
u_path, v_path = map(u -> u[1], U), map(u -> u[2], U)
scatter!(u_path, v_path, Φ.(u_path, v_path), label=L"\Phi(u,v)=u^2/2+v^2/2", legend=:outerbottom, markersize=2, markerstrokewidth=0, color=:black)
title!("Finding a Local Minimum")
savefig("box-constrained-optimisation.png")
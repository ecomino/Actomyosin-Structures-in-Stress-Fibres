include("main.jl")

cf = main(Parameters(T=1000, δ=1.5),true,false,"");

D = collect(0:0.1:0.6)
conf = []
for d in D
    cf = main(Parameters(N = 2, M = 1, L = 1, B = 1.25, T=110, δ=d),false,false,"")
    push!(conf,cf)
end
labs = "δ = " .* string.(D')
plot(1:110,conf,title="Varying Motor Drop-Off Distance δ",xlabel="Time",ylabel="Contractile Force",labels=labs)

# t = 1:1000
# cf, pos = main(Parameters(T=1000,α=0.001), false, true, "");
# plot(pos,t,xlims=(0,8),lc=[:blue :blue :red :red],legend=false)
# xlabel!("Position of Filament on Stress Fibre")
# ylabel!("Time")
# title!("Filament Movement (with turnover)")
# num_trials = 10
# I = collect(8:2:30)
# means = Vector{Float64}(undef,length(I))
# for (ind_i,i) in enumerate(I)
#     print("$i...")
#     contractile_force_trials = Vector{Vector{Float64}}(undef, num_trials)
#     for k = 1:num_trials
#         contractile_force_trials[k] = main(Parameters(M=i),false,false,"")
#     end
#     steady_contractile_forces = mean.(contractile_force_trials)
#     means[ind_i] = mean(steady_contractile_forces)
# end
# scatter(I,means,xlabel="Number of Motors (M)",ylabel="Contractile Force",title="Effect of Number of Motors on Contractile Force",legend=false,ms=3,mc=:red,msw=0.01)

# num_trials = 10
# I = collect(10:10:40)
# J = collect(4:2:10)
# heat = Matrix{Float64}(undef,length(I),length(J))
# for (ind_j, j) in enumerate(J)
#     print("L = $j... N = ")
#     for (ind_i,i) in enumerate(I)
#         print("$i,")
#         contractile_force_trials = Vector{Vector{Float64}}(undef, num_trials)
#         for k = 1:num_trials
#             contractile_force_trials[k] = main(Parameters(L=1,B=j,N=i,M=10,T=200),false,false,"")
#         end
#         steady_contractile_forces = mean.(contractile_force_trials)
#         heat[ind_i,ind_j] = mean(steady_contractile_forces)
#     end
#     println("done.")
# end
# gr()
# heatmap(I, J, heat, c=cgrad([:red,:white]), xlabel="Number of Filaments (N)", ylabel="Length of Stress Fibre (L)",title="Combined Effect of N and L on Contractile Force")
# savefig("n-b-heat.pdf")
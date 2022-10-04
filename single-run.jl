include("main.jl")

cf = main(Parameters(),false,true,"no-more-globals!")

# conf = []
# for i = 1:5
#     cf = main(true,false,"new-params-$i")
#     push!(conf,cf)
# end
# plot(1:T,conf,title="Longer Fibre, More Filaments, More Motors",xlabel="Time",ylabel="Contractile Force",legend=false)
# savefig("new-params-cf.pdf")

# conf = []
# for i in 1:3
#     cf = main(false,false,"")
#     push!(conf,cf)
# end
# plot(1:T,conf,title="N = 20 and B = 8.0",xlabel="Time",ylabel="Contractile Force",legend=false)
# savefig("b-n-2.pdf")


# I = [10, 20, 30, 40]
# J = [4.0,6.0,8.0,10.0]
# for j in J
#     println("B = $j")
#     global B = j
#     for i in I
#         println("M = $i")
#         conf = []
#         global M = i
#         for k = 1:3
#             cf = main(false,false,"M-$i")
#             push!(conf,cf)
#         end
#         # plot(1:T,conf,title="η = $i and α = $j",xlabel="Time",ylabel="Contractile Force",labels=string.(I'),legend=:outerright)
#         plot(1:T,conf,title="M = $i and B = $j",xlabel="Time",ylabel="Contractile Force",legend=false)
#         savefig("M-$i-B-$j-cf.pdf")
#     end
# end

num_trials = 10
I = collect(0.0:0.2:1.8)
J = 10 .^ (-4:0.5:0)
heat = Matrix{Float64}(undef,length(I),length(J))
for (ind_j, j) in enumerate(J)
    print("β = $j... υ = ")
    global β = j
    for (ind_i,i) in enumerate(I)
        print("$i,")
        global υ = i
        contractile_force_trials = Vector{Vector{Float64}}(undef, num_trials)
        for k = 1:num_trials
            contractile_force_trials[k] = main(false,false,"υ-$i")
        end
        steady_contractile_forces = mean.(contractile_force_trials)
        heat[ind_i,ind_j] = mean(steady_contractile_forces)
    end
    println("done.")
end
gr()
heatmap(I,J, heat, c=cgrad([:red,:white]), xlabel="Motor Dropoff (δ)", ylabel="Motor Turnover (β)",title="Combined Effect of β and δ on Contractile Force",yaxis=:log)
savefig("heat-beta-delta.pdf")
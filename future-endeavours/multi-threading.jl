include("main.jl")

function start_sim!(mean_contractile_forces,params,l=ReentrantLock())
    println("Thread spawned: α = $(params.α), δ = $(params.δ)...")
    cf = main(params)
    lock(l)
    try
        push!(mean_contractile_forces,mean(cf))
    finally
        unlock(l)
        println("Done: α = $(params.α), δ = $(params.δ).")
    end
end

num_trials = 1
D = [0.0, 0.1]
A = [0.01, 0.1]
# D = collect(0.0:0.2:1.8)
# A = 10 .^ (-4:0.5:0)
heat = Matrix{Float64}(undef,length(D),length(A))
for (j, a) in enumerate(A)
    for (i,d) in enumerate(D)
        contractile_forces = Float64[]
        Threads.@threads for k = 1:num_trials
            contractile_forces[k] = Threads.@spawn start_sim!(contractile_forces, Parameters(α = j, δ = i))
        end
        heat[i,j] = mean(contractile_forces)
    end
end

Threads.@threads for i = 1:10            
        a[i] = Threads.threadid()        
    end

# @show heat
# gr()
# heatmap(D, A, heat, c=cgrad([:red,:white]), xlabel="Motor Dropoff (δ)", ylabel="Filament Turnover (α)",title="Combined Effect of α and δ on Contractile Force",yaxis=:log)
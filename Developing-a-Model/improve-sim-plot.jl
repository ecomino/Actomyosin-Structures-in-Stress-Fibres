# Old plot sim function which takes 2secs per iteration to run (N=30, M=15)
function plot_sim_old(x,y,t)
    X_span = A-L/2:B+L/2 # Range for filament movement
    # Set up domain
    plot(xlims=(first(X_span), last(X_span)), ylims=(0,N), title="Time = $t")
    # Display filaments
    for n = 1:N
        a = x[n]-L/2 # Start of filament
        b = x[n]+L/2 # End of filament
        line_colour(n) = P[n] > 0 ? :blue : :red
        plot!([a,b], [n,n], legend=false, lc=line_colour(n))
    end
    # Display motors and their attachments
    for m = N+1:N+M
        cf = y[m-N]
        if !isempty(cf)
            scatter!([x[m]], [mean(cf)], markercolor=:black, markerstroke=0, markersize=2.5)
            plot!([x[m];x[m]], [cf[1];cf[end]], legend=false, lc=:black)
        end
    end
    # Display focal tesions
    plot!([x[end-1]-L/2;x[end-1]+L/2], [N/2;N/2], legend=false, lc=:black, linewidth=2)
    plot!([x[end]-L/2;x[end]+L/2], [N/2;N/2], legend=false, lc=:black, linewidth=2)
    plot!(show=true)
    frame(anim)
end

const A = 0 # Left end point
const B = 5 # Right end point
const N = 20 # Number of actin filaments
const L = 1 # Length of filaments
const M = 10 # Maximum number of motors
x = vcat(B .*rand(N), B .*rand(M),A,B)
y = [[1, 5], [10, 17], [17, 19], [6, 9], [7, 15], [], [2, 9], [13, 15], [], [7, 15]]
t = 1

X_span = A-L/2:B+L/2 # Range for filament movement
# Set up domain
plot(xlims=(first(X_span), last(X_span)), ylims=(0,N), title="Time = $t")
# Display filaments
i = 1:N
a = x[i] .- L/2 # Start of filaments
b = x[i] .+ L/2 # End of filaments
line_colour = reshape([p > 0 ? :blue : :red for p in P], 1, length(P)) # Polarity of filaments
plot!([a b]', [i i]', legend=false, lc=line_colour)
# Display motors and their attachments
i = N+1:N+M




function plot_sim(x,y,t)
    X_span = A-L/2:B+L/2 # Range for filament movement
    # Set up domain
    plot(xlims=(first(X_span), last(X_span)), ylims=(0,N), title="Time = $t")
    # Display filaments
    n = 1:N
    a = x[n] .- L/2 # Start of filaments
    b = x[n] .+ L/2 # End of filaments
    line_colour = reshape([p > 0 ? :blue : :red for p in P], 1, length(P)) # Polarity of filaments
    plot!([a b]', [n n]', legend=false, lc=line_colour)
    # Display motors and their attachments
    for m = N+1:N+M
        cf = y[m-N]
        if !isempty(cf)
            scatter!([x[m]], [mean(cf)], markercolor=:black, markerstroke=0, markersize=2.5)
            plot!([x[m];x[m]], [cf[1];cf[end]], legend=false, lc=:black)
        end
    end
    # Display focal tesions
    plot!([x[end-1]-L/2;x[end-1]+L/2], [N/2;N/2], legend=false, lc=:black, linewidth=2)
    plot!([x[end]-L/2;x[end]+L/2], [N/2;N/2], legend=false, lc=:black, linewidth=2)
    plot!(show=true)
end

plot_sim(x,y,t)
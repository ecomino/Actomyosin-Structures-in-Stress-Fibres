using Plots, RCall
include("model-parameters.jl")

function plot_sim(x,y,t)
    padding = L
    X_span = (A-padding,B+padding) # Range for filament movement
    # Set up domain
    plot(xlims=X_span, ylims=(0,N+1), title="Time = $t")
    # Display filaments
    a = x[filaments] .- L/2 # Start of filaments
    b = x[filaments] .+ L/2 # End of filaments
    line_colour = reshape([p > 0 ? :blue : :red for p in P], 1, length(P)) # Polarity of filaments
    plot!([a b]', [filaments filaments]', legend=false, lc=line_colour)
    # Display motors and their attachments
    for m = motors
        af = y[m-N]
        if isempty(af)
            scatter!([x[m]], [0], markercolor=:black, markerstroke=0, markersize=3)
        else
            scatter!([x[m]], [mean(af)], markercolor=:black, markerstroke=0, markersize=3)
            plot!([x[m];x[m]], [af[1];af[end]], legend=false, lc=:black)
        end
    end
    # Display focal adhesions
    N % 2 == 0 ? height = (N+1)/2 : height = N/2 # Ensure focal adhesions are not overlapping a filament
    plot!([x[focal_adhesions[1]]-l/2 x[focal_adhesions[2]]-l/2;x[focal_adhesions[1]]+l/2 x[focal_adhesions[2]]+l/2], [height height; height height], legend=false, lc=:black, linewidth=2)
    plot!(show=true)
    frame(anim)
end

"Constructing LOESS Regression plot using RCall package, and LOESS Regression Package in R"
function loess_plot(statistics::Vector{Vector{Float64}}, parN::Numerical_Parameters)
    gr(); plot(); # Load GR plotting back-end and clear previous plots
    mean_stats = mean(statistics); # Take time-averaged value across various trials
    time::Vector{Float64} = (0:(parN.nT-1))*parN.dt;
    @rput mean_stats time
    R"""
    library(spatialEco)
    loess_model = loess.ci(mean_stats, time, p=0.99, plot = FALSE, span=0.6)
    plot(time, mean_stats, xlab="Time, t (s)", ylab="Contractile Force", main="LOESS regression")
    lines(time, loess_model$loess, col="red")
    lines(time, loess_model$uci, col="red", lty=2)
    lines(time, loess_model$lci, col="red", lty=2)
    polygon(x = c(time, rev(time)),
            y = c(loess_model$lci,
                rev(loess_model$uci)),
            col =  adjustcolor("red", alpha.f = 0.10), border = NA)
    """
end;
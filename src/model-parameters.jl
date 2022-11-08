Base.@kwdef mutable struct Parameters
    # Simulation parameters
    A::Float64 = 0.0 # Left end point
    B::Float64 = 8.0 # Right end point
    Δt::Float64 = 0.1 # Step size
    T::UInt16 = 250 # Number of time steps

    # Filament parameters
    N::UInt16 = 20 # Number of actin filaments, use odd number to ensure they do not visually overlap with focal adhesions
    filaments::UnitRange{UInt16} = 1:N # Indexes of filaments
    L::Float64 = 2.0 # Length (μm) of filaments
    P::Vector{Int64} = rand((-1,1),N) # Polarities of filaments, 1 represents (left to right) minus end to plus end
    α::Float64 = 0.0 # Probability filament turns over

    # Motor parameters
    M::UInt16 = 10 # Maximum number of myosin motors
    motors::UnitRange{UInt16} = N+1:N+M # Indexes of motors
    λ::UInt16 = ceil(Int,0.75*M) # Average number of motors currently attached
    β::Float64 = 0.0 # Probability motor turns over
    δ::Float64 = 0.0 # Motor drop off point, how far the motor has to be from a connected filament's plus end to detach

    # Focal adhesion parameters
    l::Float64 = L/2 # Length (μm) of focal adhesions
    focal_adhesions::UnitRange{UInt16} = N+M+1:N+M+2 # Indexes of focal adhesions
    structure::Vector{UInt16} = [filaments; focal_adhesions] # Indexes of all filaments and focal adhesions

    # Energy functional parameters
    ξ::Float64 = 0.0000001 # Drag on filaments due to moving through the cytoplasm (pNs/μm^2)
    η::Float64 = 15.0 # Effective viscous drag (pNs/μm^2) due to cross-linker proteins
    Fs::Float64 = 5.0 # Motor stall force (pN)
    Vm::Float64 = 0.5 # Maximum (load-free) motor velocity (μm/s)
    ρ::Float64 = 10.0 # Drag on filaments due to climbing through extracellular matrix (pNs/μm^2)
    k::Float64 = 10.0 # Stiffness between focal adhesions

    anim::Animation = Animation()
end;
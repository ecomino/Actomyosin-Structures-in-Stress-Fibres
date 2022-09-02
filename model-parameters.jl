# Simulation parameters
const A = 0.0 # Left end point
B = 10.0 # Right end point
const Δt = 0.1 # Step size
const T = 200 # Number of time steps

# Filament parameters
N = 30 # Number of actin filaments, use odd number to ensure they do not visually overlap with focal adhesions
filaments = 1:N # Indexes of filaments
L = 2.0 # Length (μm) of filaments
P = rand((-1,1),N) # Polarities of filaments, left to right 1 represents -ve to +ve
α = 0.01 # Probability filament turns over

# Motor parameters
M = 10 # Maximum number of myosin motors
motors = N+1:N+M # Indexes of motors
λ = floor(Int,0.75*M) # Average number of motors currently attached
const β = 0.05 # Probability motor detaches per second
υ = 0.0 # Motor drop off point, how far the motor has to be from a connected filament's plus end to detach
@assert υ < L/2 && υ ≥ 0.0

# Focal adhesion parameters
const l = L/2 # Length (μm) of focal adhesions
focal_adhesions = N+M+1:N+M+2 # Indexes of focal adhesions
structure = [filaments; focal_adhesions] # Indexes of all filaments and focal adhesions

# Energy functional parameters
const ξ = 0.0000001 # Drag on filaments due to moving through the cytoplasm (think about units!)
η = 10.0 #15.0 # Effective viscous drag (pNs/μm^2) due to cross-linker proteins
const Fs = 5.0 # Motor stall force (pN)
const Vm = 0.5 # Maximum (load-free) motor velocity (μm/s)
ρ = 10.0 # Drag on filaments due to climbing through extracellular matrix (think about units!)
k = 10.0 # Stiffness between focal adhesions

# Output gif
anim = Animation();
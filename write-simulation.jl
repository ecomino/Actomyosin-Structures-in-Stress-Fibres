initialise_dataframe() = DataFrame(
    Time=Int64[], 
    FilamentPos=Vector{Vector{Float64}}(), 
    MotorPos=Vector{Vector{Float64}}(), 
    FocalTesionPos=Vector{Vector{Float64}}(), 
    MotorConnections=String[], 
    ContractileForce=Float64[], 
    Other=String[]
)

function serialize_motor_connections(Y)
    output = "[ "
    for y in Y
        if isempty(y)
            output *= "[] "
        else
            output *= "[$(y[1]),$(y[2])] "
        end
    end
    return output * "]"
end;

sim_status(t,x,y,cf,f,m,a) = (t, round.(x[f],digits=6), round.(x[m],digits=6), round.(x[a],digits=6), serialize_motor_connections(y), round(cf,sigdigits=6), "")
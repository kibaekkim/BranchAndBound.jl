module BranchAndBound

using JuMP
using MathOptInterface

const MOI = MathOptInterface

include("node.jl")
include("tree.jl")

# specific implementation
include("impl-JuMP.jl")

end # module

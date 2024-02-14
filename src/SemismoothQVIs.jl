module SemismoothQVIs

using Gridap, LineSearches, LinearAlgebra

export GeneralizedThermoformingQVI,
    fixed_point, newtonss, h1,
    FEM_1D_model, FEM_2D_model


include("nls.jl")
include("thermoforming.jl")
include("fem.jl")


end # module SemismoothQVIs

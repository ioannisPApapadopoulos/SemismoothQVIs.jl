module SemismoothQVIs

using Gridap, LineSearches, LinearAlgebra

# include("models.jl")

export GeneralizedThermoformingQVI,
    fixed_point, semismoothnewton, moreau_yosida_newton, h1, EOC,
    fem_model,
    latex_table

include("fem.jl")
include("nls.jl")
include("thermoforming.jl")
include("misc.jl")




end # module SemismoothQVIs

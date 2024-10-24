module SemismoothQVIs

using Gridap, LineSearches, LinearAlgebra

# include("models.jl")

export GeneralizedThermoformingQVI, NonlinearVI,
    fixed_point, semismoothnewton, moreau_yosida_newton, visolver,
    h1, h10, EOC,
    fem_model,
    latex_table, test2_latex_table, test4_latex_table

include("fem.jl")
include("nls.jl")
include("thermoforming.jl")
include("nonlinearvi.jl")
include("misc.jl")
include("linesearch.jl")




end # module SemismoothQVIs

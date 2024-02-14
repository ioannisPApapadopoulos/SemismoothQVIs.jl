using Gridap

function FEM_1D_model(n::Int)
    model = CartesianDiscreteModel((0,1),(n,))
    labels = get_face_labeling(model)
    p = 1
    reffe_u = ReferenceFE(lagrangian,Float64,p)

    Vu = TestFESpace(model,reffe_u,labels=labels,dirichlet_tags="boundary",conformity=:H1)
    VT = TestFESpace(model,reffe_u,labels=labels,conformity=:H1)
    Uu = TrialFESpace(Vu, 0.0)
    UT = TrialFESpace(VT)
    Ω = Triangulation(model)
    dΩ = Measure(Ω,p+20)
    return (dΩ, Uu, Vu, UT, VT)
end

function FEM_2D_model(n::Int)
    domain = (0,1,0,1)
    partition = (n,n)
    model = CartesianDiscreteModel(domain,partition)
    labels = get_face_labeling(model)
    p = 1
    reffe_u = ReferenceFE(lagrangian,Float64, p)

    Vu = TestFESpace(model,reffe_u,labels=labels,dirichlet_tags="boundary",conformity=:H1)
    VT = TestFESpace(model,reffe_u,labels=labels,conformity=:H1)
    Uu = TrialFESpace(Vu, 0.0)
    UT = TrialFESpace(VT)
    Ω = Triangulation(model)
    dΩ = Measure(Ω,p+10)
    return (dΩ, Uu, Vu, UT, VT)
end
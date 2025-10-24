module SOpt

    # Carrega as dependencias
    using LFrame
    using LinearSolve
    using Distributions
    using LinearAlgebra
    #using Test
    using OrderedCollections
    using Lgmsh
    using YAML
    using DelimitedFiles
    using LASS
    using ForwardDiff
    using WallE
    using Plots


    # Carrega os modulos
    include("validacoes.jl")
    include("main_otim.jl")
    include("driver_otm.jl")
    include("forcas.jl")
    include("tensoeS.jl")
    include("volume.jl")
    include("LASS/lass.jl")    
    include("auxiliar.jl")
    #include("tensao.jl")

    export main_otim


end 
module SOpt

    # Carrega as dependencias
    using LFrame
    using LinearSolve
    using Distributions
    using LinearAlgebra
    using Test
    using OrderedCollections
    using Lgmsh
    using YAML
    using DelimitedFiles
    using LASS
    using Enzyme
    using ForwardDiff


    # Carrega os modulos
    #include("main.jl")
    include("main_otim.jl")
    include("tensao2.jl")
    include("forcas.jl")
    include("driver_analise_tensao.jl")
    #include("LASS/lass.jl")    
    #include("auxiliar.jl")
    #include("tensao.jl")

    export main_otim


end 
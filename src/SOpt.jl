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
    using ForwardDiff


    # Carrega os modulos
    #include("main.jl")
    include("main_otim.jl")
    include("tensoes.jl")
    include("forcas.jl")
    include("restricao_global.jl")
    include("volume.jl")
    #include("LASS/lass.jl")    
    #include("auxiliar.jl")
    #include("tensao.jl")

    export main_otim


end 
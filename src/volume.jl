#############################################################################################################
#                             CÁLCULOS REFERENTES A RESTRIÇÃO DE VOLUME                                     #
#############################################################################################################

# Não seeeeeei onde deixar essas funções então vou por aqui mesmo :D
# Calculando o volume

function Volume(ne::Int64, dicionario_geometrias, dicionario_materiais, L::Vector, ρ::Vector, dados_elementos::Matrix{String})
    # Alocando o volume
    V = 0.0

    # loop pelos elementos
    for e=1:ne

        # Dados que vamos precisar:
        Ize, Iye, J0e, Ae, αe, Ee, Ge = LFrame.Dados_fundamentais(e, dados_elementos, dicionario_materiais, dicionario_geometrias)

        # Descobre ρ do elemento
        ρe = ρ[e]

        # Volume do elemento
        Ve = Ae*L[e]*ρe

        # Acumula o volume
        V += Ve

    end

    return V

end 

# Calculando a derivada do volume
function Derivada_volume(ne::Int64, dicionario_geometrias, dicionario_materiais, L::Vector{Float64}, dados_elementos::Matrix{String})
    
    # Aloca dV
    dv = zeros(ne)

    # Loop pelos elementos
    for e=1:ne

        # Dados que vamos precisar:
        Ize, Iye, J0e, Ae, αe, Ee, Ge = Dados_fundamentais(e, dados_elementos, dicionario_materiais, dicionario_geometrias)

        # dV fica:
        dv[e] = L[e]*Ae

    end

    return dv

end
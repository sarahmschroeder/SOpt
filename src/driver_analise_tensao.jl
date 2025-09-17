#
# Monta o sistema KU=F .... INFLUENCIADO POR ρ
#
function Monta_sistema(ρ::AbstractVector{T}, malha::LFrame.Malha) where T

    # Monta a matriz de rigidez global
    KG = LFrame.Monta_Kg(malha,ρ)

    # Número de nós 
    nnos = malha.nnos

    # Monta o vetor global de forças concentradas - não muda
    # e não precisamos dele exatamente agora...só para alocar 
    # no linsolve
    FG = zeros(6*nnos) #Monta_FG(forcas,nnos)

    # Modifica o sistema para considerar as restrições de apoios 
    KA, FA = LFrame.Aumenta_sistema(malha, KG, FG)

    # Cria um problema linear para ser solucionado pelo LinearSolve
    prob = LinearSolve.LinearProblem(KA,FA)
    linsolve = LinearSolve.init(prob,KLUFactorization())

    # Retorna o sistema linear 
    return linsolve

end


# 
#
#
#             Driver para receber um conjunto de variáveis de projeto e devolver o objetivo e a restrição global de tensões
#
#
# ρ -> vetor com ne variáveis de projeto
#
# bins <- Generate_bins(realizacoes, Nb) <-  realizacoes <- gera_distribuicoesforcas(malha, n_r=100, σ2=0.4)
#
# 
#


#
# 
#
function Driver_f(ρ::AbstractVector{T}, x::AbstractVector, malha::LFrame.Malha, forcas::AbstractMatrix) where T

    #
    # O cáculo da resposta aleatória não depende de alteração do ρ
    # Assim, podemos montar um problema linear e modificar somente
    # o r.h.s do KU = F
    linsolve = Monta_sistema(ρ,malha)

    # Cálculo da norma das tensões
    Realiza_norma_σe(x,malha,forcas,linsolve)

end



#
# Função que devolve a norma p das σe
#
# x são as variáveis aleatórias
#
function Realiza_norma_σe(x::AbstractVector,  malha::LFrame.Malha, forcas::AbstractMatrix,  linsolve, P=2.0)

    # Aplica as forças
    aplica_loads!(forcas, x)

    # Número de nós 
    nnos = malha.nnos

    # Monta o vetor de forças 
    F = Monta_FG(forcas,nnos)

    # Modifica o rhs
    linsolve.b[1:6*malha.nnos] .= F

    # Calcula o deslocamento 
    U = LinearSolve.solve!(linsolve)[1:6*malha.nnos]

    # Evita zeros em U
    #
    # XUNXÂÂÂÂÂÂO PARA EVITAR DIVISÕES POR ZERO NAS DF
    #
    U .+= 1E-12
      
    # Calcula tensoes equivalentes
    σe = tensao_equivalente(U, malha)

    # Calcula norma P
    norma = norm(σe,P)

    # Devolve
    return norma

end


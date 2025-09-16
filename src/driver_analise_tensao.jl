#
# Monta o sistema KU=F .... INFLUENCIADO POR ρ
#
function Monta_sistema(ρ, malha::LFrame.Malha)

    # Monta a matriz de rigidez global
    KG = LFrame.Monta_Kg(malha,ρ)

    # Monta o vetor global de forças concentradas - não muda
    FG = LFrame.Monta_FG(malha)

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
function Driver_f(ρ::AbstractVector{T}, x::AbstractVector, malha::LFrame.Malha) where T

    #
    # O cáculo da resposta aleatória não depende de alteração do ρ
    # Assim, podemos montar um problema linear e modificar somente
    # o r.h.s do KU=F
    linsolve = Monta_sistema(ρ,malha)

    # Cálculo da norma das tensões
    Realiza_norma_σe(x,malha,linsolve)

    # Dada a distribuição das variáveis de projeto, calcula a resposta da 
    # tensão equivalente com o LASS
    #Eσe, Varσe = Lass(bins,  x -> f(x))  

end





#
# Função que devolve a norma p das σe
#
# x são as variáveis aleatórias
#
function Realiza_norma_σe(x::AbstractVector,  malha, linsolve, P=2.0)

    # Aplica as forças
    aplica_loads!(malha, x)

    # Monta o vetor de forças 
    F = LFrame.Monta_FG(malha)

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


#=
# Função que devolve a derivada da norma p das σe
function Realiza_derivada_norma_σe(x::Vector,  malha, ρ::Vector, P=8.0)

    # Aplica as forças
    aplica_loads!(malha, x)

    # Calcula a resposta da estrutura
    U,_ = Analise3D(malha,false;ρ0=ρ)

    # Aloca os vetores Adjunto e da derivada parcial 
    Fλ = similar(U)
    ∂D = zeros(malha.ne)

    # 

    # Monta o vetor de carregamento adjunto e também o vetor das derivadas
    # parciais da norma em relação as variáveis de projeto
    #
    # Fλ, D <- .....
    # 
    #

    # Soluciona o problema adjunto 
    #
    # λ <- solução de K λ = Fλ
    #

    # Ultimos calculos e devolve a derivada

    # Calcula tensoes equivalentes
    #σe = tensao_equivalente(U, malha)

    # Calcula norma P
    #norma = norm(σe,P) # sum((σ^P))^(1/P)

    # Devolve
    return Dnorma

end


=#
############################################################################################################################
#                                          ROTINA DO DRIVER DE OTIMIZAÇÃO                                                  #
############################################################################################################################

function Driver(ρ::AbstractVector{T}, bins, r0::Float64, malha::LFrame.Malha, μ::Vector, sigma_y::Float64,
                m::Int64, ne,nnos,elems,dados_elementos,dicionario_materiais, 
                dicionario_geometrias,L,coord, loads,floads, apoios, mpc, deslocamentos, β = 3.0,
                opcao::String)


    #Faz a verificação da opção
    opcao in ["LA","dLA","g","U"] || error("Driver::opção $opcao inválida")


    ####################################### EQUILIBRIO ##############################################

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
        if opcao=="U"
            return linsolve
        end

    end


    ################################### FUNÇÃO OBJETIVO #################################################
    # Calcula o volume da estrutura
    V = Volume(ne,dicionario_geometrias,L,ρ, dados_elementos)


    ################################## RESTRIÇÃO DE TENSÃO #############################################

    #
    # Função que faz o cálculo da restrição aleatória
    #
    function Realiza_gσ(ρ::AbstractVector{T}, x::AbstractVector, malha::LFrame.Malha, forcas::AbstractMatrix, σ_Y) where T

        #
        # O cáculo da resposta aleatória não depende de alteração do ρ
        # Assim, podemos montar um problema linear e modificar somente
        # o r.h.s do KU = F
        linsolve = Monta_sistema(ρ,malha)

        # Cálculo da norma das tensões
        Calcula_gσ(x,malha,forcas,linsolve,σ_Y)

        # Devolve
        return gσ
    
    end

    #
    # Função que devolve a restrição global de tensão
    #
    # x são as variáveis aleatórias
    # 
    function Calcula_gσ(x::AbstractVector,  malha::LFrame.Malha, forcas::AbstractMatrix,  linsolve, σ_Y, P=2.0)

        # Aplica as forças
        aplica_loads!(forcas, x)

        # Número de nós 
        nnos = malha.nnos

        # Monta o vetor de forças 
        F = Monta_FG(forcas,nnos)

        # Modifica o rhs
        linsolve.b[1:6*malha.nnos] .= F

        # Calcula o deslocamento 
        sol = LinearSolve.solve!(linsolve)
        U = sol.u[1:6*malha.nnos]

        # Evita zeros em U
        #
        # XUNXÂÂÂÂÂÂO PARA EVITAR DIVISÕES POR ZERO NAS DF
        #
        U .+= 1E-12
        
        # Calcula tensoes equivalentes
        σe = tensao_equivalente(U, malha)

        # Calcula norma P
        norma = norm(σe,P)
        g = norma/σ_Y - 1


        return gσ
    

    end

    gσ = Calcula_gσ(ρ, x, malha, forcas, σ_Y)

    if op == "gσ"
        return gσ
    end

    # Objetivo:
    LA = (r0/2)*Heaviside(μ[i]/r0 + gσ)^2

    if opcao=="LA"
        return LA
    end


    ################################################# DERIVADAS ###################################################

    # Derivada da função objetivo
    dV = Derivada_volume(ne, dicionario_geometrias, L, dados_elementos)

    # Calcula a derivada da restrição: no nosso caso a dE e dVar
    funcaox(x) =  Realiza_gσ(ρ0, x, malha, forcas,σ_Y)

    # Dada a distribuição das variáveis de projeto, calcula a resposta da 
    # tensão equivalente com o LASS
    Egσ, Vargσ = Lass(bins,  x -> funcaox(x))  

    ##########################################################
    # Derivada em relação ao ρ (otimização) para um x fixo
    #########################################################
    function derivada(x,ρ0,forcas)

        # Substitui o valor de x 
        funcaoρ(ρ) =  Realiza_gσ(ρ, x, malha, forcas, σ_Y)

        # Calcula a derivada
        ForwardDiff.gradient(funcaoρ,ρ0)

    end

    #
    # Calcula a derivada em relação ao ρ usando o LASS
    #
    dEgσ, dVargσ = dLass(bins,  x-> funcaox(x), x -> derivada(x,ρ0,forcas), malha.ne)

    # Derivada da LA
    dLA = dV + (r0)*Heaviside(μ[1]/r0 + g[1])*(dE + β*dVar)

    if opcao == "dLA"
        return dLA
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


end


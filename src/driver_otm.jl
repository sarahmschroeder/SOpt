############################################################################################################################
#                                          ROTINA DO DRIVER DE OTIMIZAÇÃO                                                  #
############################################################################################################################




function Driver(ρ::AbstractVector{T}, bins, r0::Float64, malha::LFrame.Malha, μ::Vector,  σ_Y::Float64,
                m::Int64,dados_elementos,dicionario_materiais, 
                dicionario_geometrias,L, β, forcas::AbstractMatrix, opcao::String) where T


    #Faz a verificação da opção
    opcao in ["LA","dLA","gσ","U"] || error("Driver::opção $opcao inválida")

    # Número de elementos
    ne = malha.ne

    #
    # Prepara o problema de equilíbrio (sem definir as forças)
    # O linsolve não é modificado ao longo da simulação
    #
    linsolve =  Monta_linsolve(ρ, malha)

    ################################### FUNÇÃO OBJETIVO #################################################

    # Volume inicial
    #V0 = Volume(ne, dicionario_geometrias, dicionario_materiais, L, ρ, dados_elementos)
    V0 = 1.0

    # Calcula o volume da estrutura
    V = Volume(ne,dicionario_geometrias,dicionario_materiais,L,ρ, dados_elementos)

    ################################## RESTRIÇÃO DE TENSÃO #############################################

    # Calcula a derivada da restrição: no nosso caso a E e Var
    funcaox(x) =  Realiza_gσ(x, malha, forcas, linsolve, σ_Y)

    # Dada a distribuição das variáveis de projeto, calcula a resposta da 
    # tensão equivalente com o LASS
    Egσ, Vargσ = Lass(bins,  x -> funcaox(x))
    
    #@show V
    #@show  Egσ, Vargσ


    # A restrição robusta é 
    gr = Egσ + β*Vargσ

    #@show gr

    if opcao == "gσ"
        return gr
    end

    # Objetivo:
    LA = V/V0 + (r0/2)*Heaviside(μ[1]/r0 + gr)^2

    #@show V/V_ESCALA, (r0/2)*Heaviside(μ[1]/r0 + gr)^2

    if opcao=="LA"
        return LA
    end


    ################################################# DERIVADAS ###################################################

    # Derivada da função objetivo
    dV = Derivada_volume(ne, dicionario_geometrias, dicionario_materiais, L, ρ, dados_elementos )

    # Calcula a derivada em relação ao ρ usando o LASS
    dEgσ, dVargσ = dLass(bins,  x-> funcaox(x), x -> derivada(x,ρ,forcas, linsolve, σ_Y, malha), malha.ne)

    # Derivada da LA
    dLA = dV./V0 + (r0)*Heaviside(μ[1]/r0 + gr)*(dEgσ + β*dVargσ)

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


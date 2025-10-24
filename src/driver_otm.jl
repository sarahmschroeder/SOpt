############################################################################################################################
#                                          ROTINA DO DRIVER DE OTIMIZAÇÃO                                                  #
############################################################################################################################

function Driver(ρ::AbstractVector, bins, r0::Float64, malha::LFrame.Malha, μ::Vector,  σ_Y::Float64,
                m::Int64,dados_elementos,dicionario_materiais, 
                dicionario_geometrias,L, β, forcas::AbstractMatrix, opcao::String)


    #Faz a verificação da opção
    opcao in ["LA","dLA","gσ","U"] || error("Driver::opção $opcao inválida")

    # Número de elementos
    ne = malha.ne

    #
    # Prepara o problema de equilíbrio (sem definir as forças)
    # O linsolve não é modificado ao longo da simulação
    #
    # linsolve =  Monta_linsolve(ρ, malha)

    ################################### FUNÇÃO OBJETIVO #################################################

    # Volume inicial
    V0 = Volume(ne, dicionario_geometrias, dicionario_materiais, L, ones(ne), dados_elementos)

    # Calcula o volume da estrutura
    V = Volume(ne,dicionario_geometrias,dicionario_materiais,L,ρ, dados_elementos)

    ################################## RESTRIÇÃO DE TENSÃO #############################################

    # Alias para calcular a restrição de tensão em função de x 
    funcaox(x) =  Realiza_gσ(x, malha, forcas, σ_Y, ρ)

    # Dada a distribuição das variáveis de projeto, calcula a 
    # média e a variância da restrição de tensão equivalente 
    # utilizando o LASS
    Egσ, Vargσ = Lass(bins,  x -> funcaox(x))
 
    # Com isso, a restrição robusta é 
    gr = Egσ + β*Vargσ

    # Se solicitado, devolve a restrição 
    if opcao == "gσ"
        return gr
    end

    # Lagrangiano aumentado
    LA = V/V0 + (r0/2)*Heaviside(μ[1]/r0 + gr)^2
   
    #@show V/V0, (r0/2)*Heaviside(μ[1]/r0 + gr)^2

    # Se solicitado, retorna o LA 
    if opcao=="LA"
        return LA
    end


    ################################################# DERIVADAS ###################################################

    # Derivada da função objetivo - Volume 
    dV = Derivada_volume(ne, dicionario_geometrias, dicionario_materiais, L, ρ, dados_elementos )

    # Calcula a derivada do valor esperado e da variância em relação ao ρ usando o LASS
    dEgσ, dVargσ = dLass(bins,  x-> funcaox(x), x -> derivada(x,ρ,forcas, σ_Y, malha), malha.ne)

    # Derivada da função Lagrangiano aumentado
    dLA = dV./V0 + (r0)*Heaviside(μ[1]/r0 + gr)*(dEgσ + β*dVargσ)
   
    # Caso solicitado, retorna a derivada
    if opcao == "dLA"
        return dLA
    end
    
end


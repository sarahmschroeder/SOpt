#######################################################################################################
#                                         ROTINA PRINCIPAL                                            #
#######################################################################################################
function main_otim(arquivo,nr=10)

    # Tensão de escoamento 
    σ_Y = 4.5 # [MPa]

    # Define o arquivo
    malha = LFrame.Le_YAML(arquivo)

    # Número de elementos na malha
    ne = malha.ne

    # Vetor de variáveis de projeto  
    ρ0 = ones(ne)

    # Restrições laterais do problema
    ρmin = 1E-3*ones(ne)
    ρmax = ones(ne)

    # Copia das informações sobre forças concentradas 
    # para podermos utilizazar uma estrutura não mutável 
    # para a malha
    forcas = malha.loads

    # Comprimento
    L = malha.L

    # Recupera as intensidades originais das forças, conforme informado no yaml
    forcas0 = forcas[:,3]

    # Define o desvio padrão
    σ2=0.9

    # Vamos gerar as realizações para utilizar ao longo da otimização 
    # matriz com nforcas × nr
    realizacoes = gera_distribuicoesforcas(forcas,forcas0,nr, σ2)

    # Grava as realizações para estudo posterior
    writedlm("realizacoes.txt", realizacoes)

    # Número de amostras por variável (força)
    n_amostras = 1

    # Forma um vetor com os numeros de bins de cada variavel
    Nb = [n_amostras for i=1:size(realizacoes,1)] 
    
    # Gera bins
    bins = Generate_bins(realizacoes, Nb)

    # Número de iterações do procedimento de otimização
    niter = 5

    # Número de restrições
    m = 1

    # Penalização inicial
    r0 = 10.0

    # μ inicial
    μ = zeros(m)

    # Número de desvios para a restrição robusta 
    β = 0.0

    # Dicionarios LFrame
    dados_ele = malha.dados_elementos
    
    dicionario_mat = malha.dicionario_materiais 
    
    dicionario_geo = malha.dicionario_geometrias

    # DRIVER
    LA(ρ) = Driver(ρ, bins, r0, malha, μ, σ_Y,m,dados_ele,dicionario_mat, 
                dicionario_geo, L, β, forcas, "LA")

    dLA(ρ) = Driver(ρ, bins, r0, malha, μ, σ_Y,m,dados_ele,dicionario_mat, 
                    dicionario_geo, L, β, forcas, "dLA")

    restr(ρ) = Driver(ρ, bins, r0, malha, μ, σ_Y,m,dados_ele,dicionario_mat, 
                    dicionario_geo, L, β, forcas, "gσ")
    
    # Loop externo do LA
    # Loop de otimização que vai alterar os ρ para minimizar a função objetivo
    # e satistfazer a restrição de tensão incerta
    for iter=1:niter

        # Início do loop interno, de otimização, que vai devolver x*
        options = WallE.Init()
        options["NITER"] = 1_000
        output = WallE.Solve(LA,dLA,ρ0,ρmin,ρmax,options)

        # Recovering solution
        ρ = output["RESULT"]
        flag_converged = output["CONVERGED"]
        opt_norm = output["NORM"]
        @show flag_converged, opt_norm
        
        # Agora vamos precisar calcular a restrição atualizada
        g = restr(ρ)

        @show g

        # Atualiza a penalização
        r0 = r0*1.1

        # Atualiza os multiplicadores
        μ .= Heaviside.(μ .+ r0*g)

        #@show g, μ, g.*μ

        # Atualiza o ponto de ótimo
        ρ0 .= ρ

        # Critério de parada seria 
        if all(g.*μ.<=1E-6)
            println("Critério de parada atingido na iteração $iter ",g.*μ)
            break
        end
        
        
    end

    return ρ0

end

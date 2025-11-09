####################################################################################################################
#                                              MAIN PARA VALIDAÇÃO                                                 #
####################################################################################################################



function main(arquivo,n_r=200)

    # Vamos primeiro calcular as tensões, sem nada do LASS ainda

    # Define o arquivo
    malha = LFrame.Le_YAML(arquivo)

    # rho fixo 1
    ρ = ones(malha.ne)

    # Monta a matriz de rigidez global
    KG = LFrame.Monta_Kg(malha,ρ)

    # Número de nós 
    nnos = malha.nnos
    forcas = malha.loads

    # Monta o vetor de forças - Eq. 5
    F = Monta_FG(forcas,nnos)
    
    # Modifica o sistema para considerar as restrições de apoios 
    KA, FA = LFrame.Aumenta_sistema(malha, KG, F)

    #=# Cria um problema linear para ser solucionado pelo LinearSolve
    prob = LinearSolve.LinearProblem(KA,FA)
    linsolve = LinearSolve.init(prob,KLUFactorization())

    # Calcula o deslocamento e retorna
    sol = LinearSolve.solve!(linsolve)
    U = sol.u[1:6*malha.nnos] =#

    # Soluciona o problema utilizando o LFrame
    U, malha = Analise3D(arquivo)

    # Mensagem
    println("Iniciando o cálculo para validação das tensões presentes na estrutura descrita em $arquivo ...")

    # Pré-aloca os vetores
    vetor_tensoes = Float64[]
    vetor_tensoes_equivalentes = Float64[]


    # realiza um looping pelos elementos
    
    for ele=1: (malha.ne)
        for no=1:2
            for pto=0:1
                    println("ELEMENTO $ele, NÓ $no E PONTO $pto:")
                    (σ_N, τ, σ_M), σ = Tensao_val_elemento_no_ponto(ele,no,pto,malha,U)

                    #@show Esforcos_internos, [σ_N;τ;σ_M], σ
                    #println("\n")
                    # Armazena os valores no vetor
                    push!(vetor_tensoes, σ_N)
                    push!(vetor_tensoes, τ)
                    push!(vetor_tensoes, σ_M)

                    push!(vetor_tensoes_equivalentes, σ)
            end
        end  
    end

    
    # Agora, vamos fazer os histogramas com as distribuições

    # Copia das informações sobre forças concentradas 
    # para podermos utilizazar uma estrutura não mutável 
    # para a malha
    forcas = malha.loads

    # Comprimento
    L = malha.L

    # Recupera as intensidades originais das forças, conforme informado no yaml
    forcas0 = forcas[:,3]

    # Define o desvio padrão do angulo
    σ3 = 30
    
    # Gera as realizações de α
    realizacoes_alpha = gera_distribuicoesalpha(forcas, n_r, σ3)

    # Grava as realizações para estudo posterior
    writedlm("realizacoes_alpha.txt", realizacoes_alpha)

    # Define o desvio padrão da força
    σ2=0.9

    # Vamos gerar as realizações para utilizar ao longo da otimização 
    # matriz com nforcas × nr
    realizacoes = gera_distribuicoesforcas(forcas,forcas0,n_r, σ2)

    # Grava as realizações para estudo posterior
    writedlm("realizacoes.txt", realizacoes)

    ############ Plots das forças e angulos:
    # Número de forças 
    nload = size(malha.loads,1)

    # FORÇAS
    xf = readdlm("realizacoes.txt")

    p = plot(title = "Distribuição das Forças", xlabel = "Valor da força [N]", ylabel = "Frequência")

   for i = 1:nload
     histogram!(p, xf[i,:], label = "Carregamento $i")
    end

    display(p)

    savefig("figuras/histograma_forca.pdf")

    # angulos
    xα = readdlm("realizacoes_alpha.txt")

    p1 = plot(title = "Distribuição dos Ângulos", xlabel = "Valor do ângulo [°]", ylabel = "Frequência")

    for i = 1:nload
     histogram!(p1, xα[i,:], label = "Ângulo do carregamento $i")
    end

    savefig("figuras/histograma_alpha.pdf")

    ###############

    # TENSÕES 
    # Numero de elementos
    nele = malha.ne

    realizacoes_tot = vcat(realizacoes, realizacoes_alpha)


    # Ve o que acontece com as tensoes
    tensoes = distribui_tensoes(malha, realizacoes_tot)

    p2 = plot(title = "Distribuição das Tensões", xlabel = "Tensão [MPa]", ylabel = "Frequência")

   for i = 1:nele
        histogram!(p2, tensoes[i,:], label = "Elemento $i")
    end

    display(p2)

    savefig("figuras/histograma_tensao.pdf")
    
    return vetor_tensoes, vetor_tensoes_equivalentes
  
end

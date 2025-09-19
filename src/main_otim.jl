#
# Função principal por enquanto pra rodarmos os Testes
#
# arquivo <- yaml com os dados da malha
# n_r     <- número de realizações
#
function main_otim(arquivo,nr=100_000)

    # Tensão de escoamento 
    σ_Y = 200.0

    # Define o arquivo
    malha = LFrame.Le_YAML(arquivo)

    # Número de elementos na malha
    ne = malha.ne

    # Vetor de variáveis de projeto  
    ρ0 = 0.5*ones(ne)

    # Copia das informações sobre forças concentradas 
    # para podermos utilizazar uma estrutura não mutável 
    # para a malha
    forcas = malha.loads

    # Recupera as intensidades originais das forças, conforme informado no yaml
    forcas0 = forcas[:,3]

    # Vamos gerar as realizações para utilizar ao longo da otimização 
    # matriz com nforcas × nr
    realizacoes = gera_distribuicoesforcas(forcas,forcas0,nr)

    # Grava as realizações para estudo posterior
    writedlm("realizacoes.txt", realizacoes)

    # Número de amostras por variável (força)
    n_amostras = 5

    # Forma um vetor com os numeros de bins de cada variavel
    Nb = [n_amostras for i=1:size(realizacoes,1)] 
    
    # Gera bins
    bins = Generate_bins(realizacoes, Nb)

    
    # Inicializa r e μ
    
    # Loop externo do LA

    # Loop interno do LA  <- aqui vai usar as funções para calcular os E,Var, dE e dVar 

    # Atualiza μ e r 


    ########################################################################################
    # 
    #                        O QUE IMPORTA PARA O CÓDIGO DE OTIMIZAÇÃO
    #
    ########################################################################################
    #
    #
    # Driver para calcular a resposta para um dado x considerando o resto fixo
    #
    #
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

    #########################################################################################



end

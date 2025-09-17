#
# Função principal por enquanto pra rodarmos os Testes
#
# arquivo <- yaml com os dados da malha
# n_r     <- número de realizações
#
function main_otim(arquivo,nr=100_000)

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

    #
    #
    # Driver para calcular a resposta para um dado ρ
    #
    #
    
    #
    #
    # Driver para calcular a resposta para um dado x
    #
    #
    funcaox(x) =  Driver_f(ρ0, x, malha, forcas)

    # Derivada em relação ao ρ (otimização) para um x fixo
    function derivada(x,ρ0,forcas)

        # Substitui o valor de x 
        funcaoρ(ρ) =  Driver_f(ρ, x, malha, forcas)

        # Calcula a derivada
        ForwardDiff.gradient(funcaoρ,ρ0)

    end

    #for r in eachcol(realizacoes)
    #    @show funcaox(r)
    #end


    # Dada a distribuição das variáveis de projeto, calcula a resposta da 
    # tensão equivalente com o LASS
    Eσe, Varσe = Lass(bins,  x -> funcaox(x))  


    # Calcula a derivada 
    dEσe, dVarσe = dLass(bins,  x-> funcaox(x), x -> derivada(x,ρ0,forcas), malha.ne)  

    @show Eσe, Varσe, dEσe, dVarσe


end

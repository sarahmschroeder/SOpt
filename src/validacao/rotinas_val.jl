###############################################################################################
#                                  ROTINAS PARA O MAIN DE VALIDACAO                           #
###############################################################################################


# Função para tensao de VM PARA VALIDAÇÃO - aqui nao tem relaxação
function Tensao_val_elemento_no_ponto(ele,no,pto,malha,U::AbstractVector; verbose=true)

    # Testes de consistência
    no in [1;2]    || error("Tensao_elemento_no_ponto::nó deve ser 1 ou 2") 
    pto in [0;1]   || error("Tensao_elemento_no_ponto::pto deve ser 0 ou 1")

    # Obtem o vetor de forças nos nós do elemento 
    Fe, Ke, Re, gls = Forcas_elemento(ele,malha,U)

    if verbose
        println("Debug no ELEMENTO $ele, NÓ $no, PONTO $pto")
        println("Vetor Fe (Forcas_elemento): $Fe") 
    end

    # Recupera as Proprieades do elemento 
    Ize, Iye, J0e, Ae, αe, Ee, Ge, geo = LFrame.Dados_fundamentais(ele, malha.dados_elementos, 
                                                                   malha.dicionario_materiais, 
                                                                   malha.dicionario_geometrias)
    
    if verbose 
        println("Ize do LFrame: $Ize")
        println("Iye do LFrame: $Iye")
        println("J0e do LFrame: $J0e")
        println("Ae do LFrame: $Ae")
    end
        
    # O raio externo pode ser obtido com 
    re = sqrt(J0e/Ae + Ae/(2*pi))
    if verbose
       println("re calculado: $re")
    end

    # Dependendo do nó, temos os esforços internos
    if no==1
       # Sinal negativo para o primeiro nó porque estamos trabalhando com esforços
       N  = -Fe[1]
       T  = -Fe[4] 
       My = -Fe[5]
       Mz = -Fe[6]
    else 
       N  =  Fe[7]
       T  =  Fe[10] 
       My =  Fe[11]
       Mz =  Fe[12]   
    end
    if verbose
       println("N (esforço interno): $N, T (esforço interno): $T, My (esforço interno): $My, Mz (esforço interno): $Mz")
    end

    # O momento resultante é 
    Mr = sqrt(My^2 + Mz^2)

    if verbose
        println("Momento Resultante (Mr): $Mr")
    end

    # Podemos calcular as componentes de tensão diretamente:

    # Barra - Eq. 14
    σ_N = N/Ae

    # Eixo - Eq. 15
    τ = re*T/J0e

    # Flexão - Eq. 16
    σ_M = ((-1)^pto)*re*Mr/Ize
    
    if verbose 
        println("σ_N calculado (final): $σ_N")
        println("τ calculado (final): $τ")
        println("σ_M calculado (final): $σ_M")
        println("\n")
    end

    # Podemos calcular a tensão equivalente de von-Mises neste
    # ponto - Eq. 33
    σe = sqrt( (σ_N+σ_M)^2 + 3*τ^2 + 1E-6^2)


    # Retorna  a tensão equivalente 
    return (σ_N,τ,σ_M),σe 

end


function Realiza_norma(x::AbstractVector, malha::LFrame.Malha, forcas::AbstractMatrix, ρ::AbstractVector)

    # Cálculo da norma das tensões considerando a variável aleatória, o resultado disso vai ser o argumento (de tensão)
    # que vai pro LASS
    σ = Calcula_norma(x,malha,forcas, ρ)

    # Devolve
    return σ
    
end

function Calcula_norma(x::AbstractVector,  malha::LFrame.Malha, forcas::AbstractMatrix, ρ::AbstractVector, P=8.0)

    # Aplica as forças
    aplica_loads!(forcas, x)

    # Número de nós 
    nnos = malha.nnos

    # Monta o vetor de forças - Eq. 5
    F = Monta_FG(forcas,nnos)

    # Soluciona o problema de equilíbrio - Eq. 5
    U = Monta_linsolve(ρ,malha,F)
    
    # Evita zeros em U
    U .+= 1E-12
    
    # Calcula tensoes equivalentes - Eq. 33
    σe = tensao_equivalente(U, malha, ρ)

    # Calcula norma P - uso de acordo com a Eq. 78 na abordagem usada
    norma = norm(σe,P)

    # retorna o valor da restrição 
    return norma
   

end

#
# Gera todas as realizações de forças para usar no LASS
#
function gera_distribuicoesforcas1(malha, n_r=100; σ2=0.4)

    # Número de forças 
    nload = size(malha.loads,1)

    # Matriz n_load × n_r com as realizações 
    realiza = zeros(nload,n_r)
    
    # Loop pelas magnitudes da malha
    for i = 1:nload
        magn  = malha.loads[i,3]
        media = magn               # média é a magnitude que tá na malha
        desvio = sqrt(abs(σ2*media))         # como definir e onde?     
        realiza[i,:] .= rand(Normal(media, desvio),n_r) 
    end

    # Retorna a matriz com todas as realizações
    return realiza

end




# Função que ve o efeito das incertezas nas forças nas tensões
function distribui_tensoes(malha, realizacoes)

    # numero de elementos
    nele = malha.ne

    # Número de realizações
    n_r = size(realizacoes, 2)

    # Gera um vetor que guarda essas informações
    #
    # Eu entendo que aqui são as tensões equivalentes (escalar)
    #
    tensoes = zeros(4*nele, n_r) 

    # Pega os pares (nó, gdl) da malha — onde aplicar cada força
    # locais = [(linha[1], linha[2]) for linha in eachrow(malha.loads)]

    #@show locais 

    # Loop por cada realização de forças (cada vetor x)
    for j in 1:n_r

        # Aplica a realização j das forças
        for i in size(malha.loads,1) 
            valor_forca = realizacoes[i, j]
            malha.loads[i,3] = valor_forca
        end
        
        # Calcula a resposta da estrutura
        U,_ = Analise3D(malha,false)

        ρ = ones(malha.ne)

        # Calcula as tensões equivalentes pra essa realização
        σ_eq = tensao_equivalente(U, malha, ρ)

        # Salva o resultado da tensão dessa realização na matriz
        tensoes[:, j] .= σ_eq
    end

    # Retorna a matriz com todas as tensões por elemento e realização
    return tensoes
end


#### FUNÇÕES EM DESUSO

function roda_lass(malha,n_r,ρ) # DESUSO

    # Define distribuições das forças 
    dists = gera_distribuicoesforcas1(malha,n_r)

    # Grava as realizações 
    writedlm("realizaforcas.txt",dists)

    # Número de amostras por variável (força)
    n_amostras = 5

    # Forma um vetor com os numeros de bins de cada variavel
    Nb = [n_amostras for i=1:size(dists,1)]  # no exemplo validacao.yaml fica [50 50 50 50]
    
    # Gera bins
    bins = Generate_bins(dists, Nb)

    #@show bins

    # Roda LASS
    E, Var = Lass(bins, x -> Realiza_norma(x, malha, ρ))  

    println("Esperança da norma: $E")
    println("Variância da norma: $Var")

    return E, Var
end


# Função para extrair os nos e gdls do carregamento ignorando a magnitude 
# (vai ser definida pelo vetor x)
function extrai_loads(malha) # DESUSO

    # Cria um vetor de Tuples
    locais = Tuple{Int64,Int64}[]

    # Loop pelas linhas de malha.loads,
    # pegando só o nó e a direção
    for (no, gdl, _) in malha.loads
        push!(locais, (no, gdl))
    end

    # Retorna o vetor de tuples
    return locais 

end



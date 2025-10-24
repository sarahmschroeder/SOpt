####################################################################################################################
#                                      MAIN/DRIVER DO PROGRAMA PRINCIPAL                                           #
####################################################################################################################


#
#
# x = readdlm("realiza.txt")
#
# histogram(x[1,:])
#
#

#
#
#
# Função principal por enquanto pra rodarmos os Testes
#
# arquivo <- yaml com os dados da malha
# n_r     <- número de realizações
#

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


function main(arquivo,n_r=2)

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

    # Cria um problema linear para ser solucionado pelo LinearSolve
    prob = LinearSolve.LinearProblem(KA,FA)
    linsolve = LinearSolve.init(prob,KLUFactorization())

    # Calcula o deslocamento e retorna
    sol = LinearSolve.solve!(linsolve)
    U = sol.u[1:6*malha.nnos]

    # Define o desvio padrão
    σ2=0.4

    # Vamos gerar as realizações para utilizar ao longo da otimização 
    # matriz com nforcas × nr
    realizacoes = gera_distribuicoesforcas1(malha,n_r)

    # Grava as realizações para estudo posterior
    writedlm("realiza.txt", realizacoes)


    #roda_lass(malha,n_r,ρ)


    # Plotagem dos resultados
    # Número de forças 
    nload = size(malha.loads,1)

    # FORÇAS
    xf = readdlm("realiza.txt")

    # p = plot()
    # for i = 1:nload
    #       histogram!(p,xf[i,:])
    # end

    for i = 1:nload
        if i == 1
            histogram(xf[i,:])
        else
            histogram!(xf[i,:])
        end
    end

    # savefig("histograma_forca.pdf")

    # TENSÕES 
    # Numero de elementos
    nele = malha.ne

    # Ve o que acontece com as tensoes
    tensoes = distribui_tensoes(malha, realizacoes)
#=
    for i = 1:nele
        if i == 1
            display(histogram(tensoes[i,:]))
        else
            display(histogram!(tensoes[i,:]))
        end
    end
    =#



    # Soluciona o problema utilizando o LFrame
    U, malha = Analise3D(arquivo)# ; ρ0) # deixando o rho pra depois

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

    return vetor_tensoes, vetor_tensoes_equivalentes
  
end

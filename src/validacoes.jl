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
function main(arquivo,n_r=2)

    # Define o arquivo
    malha = LFrame.Le_YAML(arquivo)

    roda_lass(malha,n_r)


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
    realizacoes = gera_distribuicoesforcas(malha,n_r)

    # Ve o que acontece com as tensoes
    tensoes = distribui_tensoes(malha, realizacoes)

    for i = 1:nele
        if i == 1
            display(histogram(tensoes[i,:]))
        else
            display(histogram!(tensoes[i,:]))
        end
    end



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
                    Esforcos_internos, (σ_N, τ, σ_M), σ = Tensao_elemento_no_ponto(ele,no,pto,malha,U)

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

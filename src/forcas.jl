

# Função para aplicar esse carregamento desse vetor
#
# x é um vetor com as variáveis aleatórias ( magnitude
# das forças )
#
function aplica_loads!(malha::LFrame.Malha, x::AbstractVector)


    # Teste dimensao
    if size(malha.loads,1) != length(x)
        error("Dimensao de malha.loads tem que ser o mesmo de x")
    end

    # Aplica os carregamentos com os valores do vetor x
    malha.loads[:,3] .= x
    
end


#
# Gera todas as realizações de forças para usar no LASS
#
function gera_distribuicoesforcas(malha::LFrame.Malha, forcas0, nr, σ2=0.4)

    # Número de forças 
    nload = size(malha.loads,1)

    # Matriz n_load × n_r com as realizações 
    realiza = zeros(nload,nr)
    
    # Loop pelas magnitudes da malha
    for i in LinearIndices(forcas0)

        # Magnitude original da força
        media  = forcas0[i]

        # Variância da distribuição 
        variancia = sqrt(abs(σ2*media))   

        # Gera as realizações segundo uma distribuição normal 
        realiza[i,:] .= rand(Normal(media, variancia),nr) 

    end

    # Retorna a matriz com todas as realizações
    return realiza

end
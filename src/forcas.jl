################################################################################################################
#                                ROTINAS PARA A MONTAGEM DO VETOR DE FORÇAS, APLICAÇÃO                         #
#                              DOS LOADS INCERTOS E REALIZAÇÃO DOS CARREGAMENTOS INCERTOS                      #
################################################################################################################

# Função para montar o vetor de forças para o elemento de pórtico espacial. Condizente com Eq. 5
function Monta_FG(forcas::AbstractMatrix{T},nnos::Int64) where T
    

    # Aloca o vetor global
    FG = zeros(T,6*nnos)

    # Loop pelas informações dos carregamentos concentrados
    for i in axes(forcas,1)

        # Descobre o nó
        no = Int(forcas[i,1])

        # Descobre o gl(local)
        gl = Int(forcas[i,2])

        #Descobre o valor
        valor = forcas[i,3]

        # O grau de liberdade global
        glg = 6*(no-1)+gl

        # Sobrepoe no gl
        FG[glg] = FG[glg] + valor
    end

    # Retorna o vetor
    return FG
end

# Função para aplicar esse carregamento desse vetor
#
# x é um vetor com as variáveis aleatórias ( magnitude
# das forças )
#
function aplica_loads!(forcas::AbstractMatrix{T}, x::AbstractVector{T}) where {T}


    # Teste dimensao
    if size(forcas,1) != length(x)
        error("Dimensao de malha.loads tem que ser o mesmo de x")
    end

    # Aplica os carregamentos com os valores do vetor x
    forcas[:,3] .= x
    
end


#
# Gera todas as realizações de forças para usar no LASS
#
function gera_distribuicoesforcas(forcas::AbstractMatrix{T}, forcas0::AbstractVector{T}, nr, σ2=0.4) where T

    # Número de forças 
    nload = size(forcas,1)

    # Matriz n_load × n_r com as realizações 
    realiza = zeros(T,nload,nr)
    
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
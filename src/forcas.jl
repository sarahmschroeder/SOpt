################################################################################################################
#                                ROTINAS PARA A MONTAGEM DO VETOR DE FORÇAS, APLICAÇÃO                         #
#                              DOS LOADS INCERTOS E REALIZAÇÃO DOS CARREGAMENTOS INCERTOS                      #
################################################################################################################

########################################### APLICAÇÃO DAS INCERTEZAS ###########################################

# Função para aplicar esse carregamento desse vetor
#
# x é um vetor com as variáveis aleatórias ( magnitude
# das forças )
#

# Mudei o tipo de forcas pra aceitar um vetor dentro da matriz
function aplica_loads!(forcas::AbstractMatrix{Any}, x::AbstractVector{T}, y::AbstractVector{T}) where {T}


    # Teste dimensao
    if size(forcas,1) != length(x)
        error("Dimensao de malha.loads tem que ser o mesmo de x")
    end

    # Agora também temos os alphas. Precisamos calcular os cossenos e senos.

    Fx = x .* cos.(y)
    Fy = x .* sin.(y)
    Fz = zeros(length(y))

    # Cria uma matriz n×3 com as componentes:
    F = hcat(Fx, Fy, Fz) # concatenação horizontal, certo? cada linha de F tem as
    # três componentes da força para os x e y (realização de modulo e angulo) correspondentes

    # Atualiza a matriz 'forcas' pra usar essas componentes:
    # Forca = Fi[cos(α), sin(α), 0]
    # Cria os vetores [Fx, Fy, 0] para cada linha
    for i in 1:length(x)
        forcas[i, 3] = [Fx[i], Fy[i], Fz[i]]
    end

    
end






################################################### ÂNGULO ######################################################

function gera_distribuicoesalpha(forcas::AbstractMatrix{T}, α_vec::AbstractVector{T}, nr, σ2; deterministico=true) where T

    # Número de forças 
    nload = size(forcas,1)

    # Matriz n_load × n_r com as realizações 
    realiza_alpha = zeros(T,nload,nr)
    
    # Loop pelo numero de loads -> cada força vai ter uma distribuição de angulo
    for i=1:nload

        # Magnitude original do angulo. novamente, considerando que o padrão é zero
        media  = 0.0

        if deterministico
            # Pequena perturbação numérica para evitar Δ == 0 (LASS precisa de variação entre amostras)
            # Perturbação relativa muito pequena: não afeta o problema "praticamente determinístico"
            tol = 1e-12   # sujeirinha numérica 
            realiza_alpha[i,:] .= media                # determinístico
            realiza_alpha[i,1] = media - tol
            realiza_alpha[i,2] = media + tol
        else
            # Variância da distribuição 
            # (como a media é zero, e eu queria uma variação de +- 30°, coloquei isso direto na variancia)
            variancia = deg2rad(30.0)
            # Gera as realizações segundo uma distribuição normal 
            realiza_alpha[i,:] .= rand(Normal(media, variancia), nr)
        end


    end

    # Retorna a matriz com todas as realizações
    return realiza_alpha

end


################################################### MÓDULO ######################################################

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

#= COMENTANDO A FUNÇÃO ORIGINAL PARA NAO PERDER

function aplica_loads!(forcas::AbstractMatrix{T}, x::AbstractVector{T}) where {T}


    # Teste dimensao
    if size(forcas,1) != length(x)
        error("Dimensao de malha.loads tem que ser o mesmo de x")
    end

    # Aplica os carregamentos com os valores do vetor x
    forcas[:,3] .= x
    
end
=#


#
# Gera todas as realizações de forças para usar no LASS
#
function gera_distribuicoesforcas(forcas::AbstractMatrix{T}, forcas0::AbstractVector{T}, nr, σ2; deterministico=true) where T

    # Número de forças 
    nload = size(forcas,1)

    # Matriz n_load × n_r com as realizações 
    realiza = zeros(T,nload,nr)
    
    # Loop pelas magnitudes da malha
    for i in LinearIndices(forcas0)

        # Magnitude original da força 
        media  = forcas0[i]

        if deterministico
            # Pequena perturbação numérica para evitar Δ == 0 (LASS precisa de variação entre amostras)
            # Perturbação relativa muito pequena: não afeta o problema "praticamente determinístico"
            tol = 1e-12   # sujeirinha numérica 
            realiza[i,:] .= media                # determinístico
            realiza[i,1] = media - tol
            realiza[i,2] = media + tol
        else
            # Variância da distribuição 
            variancia = sqrt(abs(σ2*media))      # robusto padrão
            # Gera as realizações segundo uma distribuição normal 
            realiza[i,:] .= rand(Normal(media, variancia), nr)
        end


    end

    # Retorna a matriz com todas as realizações
    return realiza

end
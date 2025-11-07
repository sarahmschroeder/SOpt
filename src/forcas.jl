################################################################################################################
#                                ROTINAS PARA A MONTAGEM DO VETOR DE FORÇAS, APLICAÇÃO                         #
#                              DOS LOADS INCERTOS E REALIZAÇÃO DOS CARREGAMENTOS INCERTOS                      #
################################################################################################################

########################################### APLICAÇÃO DAS INCERTEZAS ###########################################

# Função para aplicar esse carregamento desse vetor
#
# x é um vetor com as variáveis aleatórias ( magnitude
# das forças e ângulo) da realização atual 
#
function aplica_loads(forcas, x::AbstractVector{T}) where {T}

    # Numero original de forças
    nload = size(forcas,1)
    
    # Teste dimensao
    if 2*nload != length(x)
        error("Verificar dimensões")
    end

    # Aloca uma nova matriz de forcas
    forcas2 = zeros(2*nload,3)

    # Loop por informação original 
    contador = 0
    for i=1:nload

        # Ponteiro para o ângulo 
        pα = nload + i

        # Nó da força 
        no = forcas[i,1]

        # gl original 
        gl = forcas[i,2]

        # Caso o gl seja Y
        if gl==2

            # Componente y 
            fy = x[i]*cosd(x[pα])

            # Componente x
            fx = x[i]*sind(x[pα])

        # Caso seja em X    
        elseif gl==1

            # Componente y 
            fy = x[i]*sind(x[pα])

            # Componente x
            fx = x[i]*cosd(x[pα])

        else

            error("aplica_loads!:: implementar para z")

        end

        # Incrementa o contador 
        contador += 1

        # Grava uma linha para a dir x
        forcas2[contador,:] = [no 1 fx]

        # Incrementa o contador 
        contador += 1

        # Grava uma linha para a dir y
        forcas2[contador,:] = [no 2 fy]

    end
    
    # Retorna uma nova matriz de forças, com o dobro de linhas 
    return forcas2
    
end


################################################### ÂNGULO ######################################################

function gera_distribuicoesalpha(forcas, nr, σ3; deterministico=false) 

    # Número de forças 
    nload = size(forcas,1)

    # Matriz n_load × n_r com as realizações 
    realiza_alpha = zeros(Float64, nload, nr)

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
            variancia = sqrt(abs(σ3))
            
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
    FG = zeros(Float64,6*nnos)

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

        # Se a força for um vetor [Fx, Fy, Fz], usa o componente correto
        if isa(valor, AbstractVector)
            FG[glg] += valor[gl]
        else
            FG[glg] += valor
        end

        # Sobrepoe no gl
        #FG[glg] = FG[glg] + valor
    end

    # Retorna o vetor
    return FG
end

# Função para aplicar esse carregamento desse vetor
#
# x é um vetor com as variáveis aleatórias ( magnitude
# das forças )



#
# Gera todas as realizações de forças para usar no LASS
#
function gera_distribuicoesforcas(forcas::AbstractMatrix{T}, forcas0::AbstractVector{T}, nr, σ2; deterministico=false) where T

    # Número de forças 
    nload = size(forcas,1)

    # Matriz n_load × n_r com as realizações 
    realiza = zeros(Float64,nload,nr)
    
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

################################################################################################
#                        ROTINAS PARA CÁLCULOS DE TENSÕES (e forças generalizadas)             #
################################################################################################


########################## ESTOCÁSTICAS ####################################3

#
# Função que faz o cálculo da restrição *aleatória*
#
function Realiza_gσ(x::AbstractVector, malha::LFrame.Malha, forcas::AbstractMatrix, linsolve, σ_Y, ρ) 

    # Cálculo da norma das tensões considerando a variável aleatória, o resultado disso vai ser o argumento (de tensão)
    # que vai pro LASS
    gσ = Calcula_gσ(x,malha,forcas,linsolve,σ_Y, ρ)

    # Devolve
    return gσ
    
end

#
# Função que devolve a restrição global de tensão
#
# x são as variáveis aleatórias
# 
function Calcula_gσ(x::AbstractVector,  malha::LFrame.Malha, forcas::AbstractMatrix,  linsolve, σ_Y, P=8.0, ρ)

    # Aplica as forças
    aplica_loads!(forcas, x)

    # Número de nós 
    nnos = malha.nnos

    # Monta o vetor de forças 
    F = Monta_FG(forcas,nnos)

    # Modifica o rhs
    linsolve.b[1:6*malha.nnos] .= F

    # Calcula o deslocamento 
    sol = LinearSolve.solve!(linsolve)
    U = sol.u[1:6*malha.nnos]
    
    # Evita zeros em U
    U .+= 1E-12
    
    # Calcula tensoes equivalentes
    σe = tensao_equivalente(U, malha, ρ)

    # Calcula norma P
    norma = norm(σe,P)

    # A restrição determinística é
    gd = norma/σ_Y - 1

    return gd
    

end

#################### DETERMINÍSTICAS ##########################



function tensao_equivalente(U::Vector{T}, malha, ρ) where T

    # Número de elementos na malha
    ne = malha.ne

    # Loop pelos elementos da malha, calculando as 
    # tensões ...
    σ_eq = zeros(T,4*ne) #Float64[]

    # Loops por elemento/nó/pto
    cont = 1
    for ele=1:ne
       for no=1:2
           for pto=0:1
               σe =  Tensao_eq_elemento_no_ponto(ele,no,pto,malha,U, ρ)    
               σ_eq[cont] = σe/1E6 
               cont += 1
            end
        end  
    end

    return σ_eq

end



function Tensao_eq_elemento_no_ponto(ele,no,pto,malha,U::Vector{TF}, ρ; verbose=false) where TF

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

    # Barra
    σ_N = N/Ae

    # Eixo
    τ = re*T/J0e

    # Flexão 
    σ_M = ((-1)^pto)*re*Mr/Ize
    
    if verbose 
        println("σ_N calculado (final): $σ_N")
        println("τ calculado (final): $τ")
        println("σ_M calculado (final): $σ_M")
        println("\n")
    end

    # Podemos calcular a tensão equivalente de von-Mises neste
    # ponto
    σe1 = sqrt( (σ_N+σ_M)^2 + 3*τ^2 + 1E-6^2)

    # Aplica relaxação
    p = 3.0
    q = q
    σe = σe1*ρ^(q-p)

    # Retorna  a tensão equivalente no ele, nó, pto
    return σe 

end



function Forcas_elemento(ele,malha::LFrame.Malha, U::AbstractVector{T}) where T
    
    # Recupera dados da estrutura malha
    conect = malha.conect
    coord = malha.coord
    dados_elementos = malha.dados_elementos
    dicionario_materiais = malha.dicionario_materiais
    dicionario_geometrias = malha.dicionario_geometrias

    # Recupera os dados do elemento
    Ize, Iye, J0e, Ae, αe, Ee, Ge, _ = LFrame.Dados_fundamentais(ele, dados_elementos, dicionario_materiais, 
                                                           dicionario_geometrias)
    
    
    # Descobre os dados do Elemento
    Le = malha.L[ele]
                                                                                                                    
    # Monta a matriz do elemento no sistema local
    Ke = LFrame.Ke_portico3d(Ee, Ize, Iye, Ge, J0e, Le, Ae)

    # Monta a matriz de rotação do elemento
    Re = LFrame.Rotacao3d(ele, conect, coord, αe)

    # Descobre os gls globais do elemento 
    gls = LFrame.Gls(ele,conect)   
 
    # Recupera os deslocamentos do elemento (ainda no sistema global)
    ug = U[gls]

    # Rotaciona para o sistema local
    ul = Re*ug

    # Multiplica pela rigidez (também no sistema local)
    fe = Ke*ul 

    # Devolve as forças generalizadas nos nós deste elemento
    # a rigidez, a matriz de rotação e também os gls
    return fe, Ke, Re, gls

end

#####################################################################################################
#                        MONTAGEM DE MATRIZES AUXILIARES E OUTRAS FUNÇÕES AUX                       #
#####################################################################################################

# Função Heaviside
Heaviside(a) = max(a,0.0)

#
# Monta o sistema KU=F .... INFLUENCIADO POR ρ
#
function Monta_linsolve(ρ::AbstractVector, malha::LFrame.Malha, F::AbstractVector)

    # Monta a matriz de rigidez global
    KG = LFrame.Monta_Kg(malha,ρ)

    # Número de nós 
    nnos = malha.nnos
    
    # Modifica o sistema para considerar as restrições de apoios 
    KA, FA = LFrame.Aumenta_sistema(malha, KG, F)

    # Cria um problema linear para ser solucionado pelo LinearSolve
    prob = LinearSolve.LinearProblem(KA,FA)
    linsolve = LinearSolve.init(prob,KLUFactorization())

    # Calcula o deslocamento e retorna
    sol = LinearSolve.solve!(linsolve)
    U = sol.u[1:6*malha.nnos]

end


# Derivada em relação ao ρ (otimização) para um x fixo
function derivada(z,ρ0,forcas, σ_Y, malha)

    # Substitui o valor de x 
    funcaoρ(ρ) =  Realiza_gσ(z, malha, forcas,  σ_Y, ρ)

    # Calcula a derivada
    ForwardDiff.gradient(funcaoρ,ρ0)

end

# Função que monta a matriz VM (Matriz utilizada no vetor de tensõesasumindo a ordem xx_N xy_T e xx_M)
function Matriz_VM()

    VM = [1.0 0.0 1.0 ;
          0.0 3.0 0.0 ;
          1.0 0.0 1.0 ]
    return VM

end

# Função que monta a matriz Pna (Matriz utilizada na derivada do vetor de tensões no elemento)
function Matriz_Pna()
    Pna =     [1/Ae    0            0;        # N
               0    re/J0e         0;        # T
               0     0    ((-1)^a)*(re/Ize)] # M
    return Pna
end

# Função que monta a matriz D (Matriz utilizada nas derivadas de tensão equivalente)
function Matriz_D()
   D =  [1     0       0       0;   # N
         0     1       0       0;   # T
         0     0     My/Mr   Mz/Mr] # M
    return D
end

function Matriz_Mn(no)
    if no == 1
        Mn = [-1   0   0   0   0   0   0   0   0   0   0   0; # N
               0   0   0  -1   0   0   0   0   0   0   0   0; # T
               0   0   0   0  -1   0   0   0   0   0   0   0; # My
               0   0   0   0   0  -1   0   0   0   0   0   0] # Mz
    else
        Mn = [0    0   0  0   0   0   1   0   0   0   0   0; # N
              0    0   0  0   0   0   0   0   0   1   0   0; # T
              0    0   0  0   0   0   0   0   0   0   1   0; # My
              0    0   0  0   0   0   0   0   0   0   0   1] # Mz
    end 
    return Mn
end
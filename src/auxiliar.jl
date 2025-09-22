#####################################################################################################
#                                      MONTAGEM DE MATRIZES AUXILIARES                              #
#####################################################################################################

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
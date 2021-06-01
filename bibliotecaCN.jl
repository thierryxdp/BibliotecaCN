using LinearAlgebra

function resolve_diagonal(D,b) #D é diagonal
    tamanho = length(b)
    x=zeros(tamanho,1)
    for i in 1:tamanho
        x[i] = b[i]/D[i,i]
    end
    return x
end

function resolve_triangular_superior(T,b) #T é triangular superior
    tamanho = length(b)
    x=zeros(tamanho,1)
    
    i = tamanho
    while (i > 0)
        x[i] = b[i]/T[i,i]
        k = tamanho
        for j in i:(tamanho - 1)
            x[i] -= T[i, k] * x[k]/T[i,i]
            k -= 1
        end
        i -= 1
    end
    return x
end

function resolve_triangular_inferior(T,b) #T é triangular superior
    tamanho = length(b)
    x=zeros(tamanho,1)
    
    i = 1
    while (i <= tamanho)
        x[i] = b[i]/T[i,i]
        for j in 1:(i-1)
            x[i] -= T[i, j] * x[j]/T[i,i]
        end
        i += 1
    end
    return x
end

function eliminacao_gaussiana(A,b) #A é uma matriz cheia ("matriz densa")
    n = length(b)
    for i = 1:n-1
        for j = i+1:n
            mult = A[j,i]/A[i,i]
            A[j,:] = A[j,:] - mult*A[i,:]
            b[j] = b[j] - mult*b[i]
        end
    end

    return A,b
end

function resolve_cheia(A,b)
    T,c= eliminação_gaussiana(A,b)  #O(n^3)
    x=resolve_triangular_superior(T,b)  #O(n^2)
    return x
end

function decomposição_LU(A) #A é uma matriz cheia ("matriz densa")
    tamanho = size(A)
    tamanho = tamanho[1]
    L=zeros(tamanho,tamanho)
    
    for i in 1:tamanho
        L[i,i] = 1
    end
    
    U = zeros(tamanho, tamanho)
    
    for i in 1:tamanho
        for j in 1:tamanho
            U[i, j] = A[i, j]
        end
    end
    
    for i in 1:tamanho-1
        for j in i+1:tamanho
            L[j,i] = U[j,i]/U[i,i]
            U[j, :] = U[j, :] -L[j,i]*U[i, :]
        end
    end

    return L,U 
end

function resolver_pós_LU(L,U,b)    #O(n^2) 
    y=resolve_triangular_inferior(L,b)
    x=resolve_triangular_superior(U,y)
    return x
end

function main()
    A=[2 4 3 ; 4 10 9 ; 6 18 22]
    b=[9 ; 23 ; 46]
    T,c = eliminacao_gaussiana(A,b,3)

    println(T)
    println(c)
end

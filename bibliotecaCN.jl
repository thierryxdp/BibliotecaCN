function eliminacao_gaussiana(A,b,n)
    for i = 1:n-1
        for j = i+1:n
            mult = A[j,i]/A[i,i]
            A[j,:] = A[j,:] - mult*A[i,:]
            b[j] = b[j] - mult*b[i]
        end
    end

    return A,b
end

function main()
    A=[2 4 3;4 10 9;6 18 22]
    b=[9; 23; 46]
    T,c = eliminacao_gaussiana(A,b,3)

    println(T)
    println(c)
end
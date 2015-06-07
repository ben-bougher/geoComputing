function Ricker(duration, dt, f)

    freq = f
    t = [-duration/2:dt:duration/2]

    pi2 = (pi ^ 2.0)
    fsqr = freq ^ 2.0
    tsqr = t .^ 2.0
    pft = pi2 * fsqr * tsqr
    A = (1 - (2 * pft)) .* exp(-pft)

    return A
end

function opNMO(n, fold, velocity, dt, dx)


    t = [0:dt:n*dt]
    x = [0:dx:fold*dx]
    x -= x[end] ./ 2.0

    NMO = [];

    for f = 1:fold

        mat = zeros(n,n)
        new_t = sqrt(t.^2 + x[f]^2 / velocity^2)

        idx = ceil(new_t ./ dt) + 1

        for i = [1:n]
            if idx[i] > n
                continue
            end 
            mat[idx[i], i] += 1
        end 

        if f== 1
            rep = zeros(fold)
            rep[f] = 1
            NMO = kron(sparse(mat), sparse(rep))
          
        else
            rep = zeros(fold)
            rep[f] = 1
            NMO = [NMO kron(sparse(mat), sparse(rep))]
        end 
    end 

    return NMO
end 

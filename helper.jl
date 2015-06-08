using PyCall
@pyimport matplotlib.pyplot as plt
using SpecialMatrices

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


function plot(data)
    
    plt.figure()
    plt.imshow(data, aspect="auto", cmap="Greys")
    
    plt.savefig("imshow.eps")
    
    plt.figure()
    for i in 1:size(data,2)
        plt.plot(data[:,i] + .3*(i-1), 1:size(data,1), "b") 
        plt.fill_betweenx(1:size(data,1), ones(size(data,1)) * .3*(i-1),
                          data[:,i] + .3*(i-1),
                          data[:,i] .> 0, edgecolor="b")
    end

    plt.title("Synthetic Angle Gather")
    plt.xlabel("Angle")
    plt.ylabel("Time")

    plt.xlim([-.3, size(data,2)*.3 + .3])
    plt.xticks([0,size(data,2)-1]*.3, [0,40])
    plt.yticks([])
    plt.tight_layout()
    plt.savefig("wiggles.eps")
end

function opToeplitz(vector)

    return sparse(full(Toeplitz(vector)))
end

function opGrad(n)
    # Returns a difference operator
    D = spdiagm(tuple(ones(n) * -1, ones(n) *1),[0,-1])
    return D'
end

function opAkiRichards(n, theta, vp_av, vs_av, rho_av)
    
    G = opGrad(n)
    
    theta = deg2rad(theta)

    A = (1/rho_av)*(1/2.)*(1- (4*((vs_av / vp_av)^2)) * sin(theta)^2)*G
    B = (1/vp_av)*(1/2)*(1+ (tan(theta)^2))*G
    C = ((1/vs_av) * 4 *((vs_av/vp_av)^2) * sin(theta)^2)*G

    return [A B C]

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

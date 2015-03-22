#=
The physics of seismic imaging is modelled by the elastic wave equation, which
in practice is reduced to an acoustic approximation. The cavelier
assumption that waves propagate through a heterogenous mechanical soup as
ripples in a pond has provided the physical footing for wavefield modelling, the
workhorse behind modern migration and inversion techniques.

This essay solves the acoustic wave-equation for a single frequency point
source, which is equivalent to solving the Helmholtz system:
 \nabla^2U + \omega^2mU = q(x)
where U is the wavefield, m is the slowness squared, \omega is the frequency,
and q is the source term

Although this example may sound trivial, superposition tells us that any waveform
can be built by a synthesis of single frequencies, and any wavefield can be
represented by a superposition of point sources. With a for loop
you can model any physical wavefield from this basic code snippet.

For the sake of simplicity, build a constant velocity model on a small 2-D
grid. To avoid floating point rounding errors, I will use km as units.
=#
    # make a 3 km by 3 km  constant velocity model
    n, dn = 40, .01 # 10 m spacing
    c = ones(n,n) * 1.5 # 1500 m/s
    w = 2.*pi*10 # 10 Hz source

#=
Now I need to make a discrete La Place operator. This will look a computational
voodoo, but keeping things in matrix form simplifies the problem. The discrete
LaPlacian operator can be represented as band diagonal matrix(ref), which we
extend to 2-D using a kronecker product(footnote).
=#    
    # Make the LaPlacian operator
    Dx = spdiagm(tuple(ones(n-1), -ones(n)*2, ones(n-1) ), [1,0,-1])
    Dx = kron(speye(n), Dx) / dn^2
    
    Dy = spdiagm(tuple(ones(n*n-n), -ones(n*n)*2, ones(n*n-n) ),
             [n,0,-n]) / dn^2
    LP = -(Dx + Dy)
    

#=
We can now solve the Helmholtz equation and see our modelled wavefield. Note
that solving this discrete PDE is solving a system a of linear equations. 
=#
    
    # Make the Helmholtz operator
    H = LP + spdiagm((w ./ m) .* 2)

    # reshape the factor back to our model dimensions
    U = reshape(u,n,n)

using PyCall
@pyimport matplotlib.pyplot as plt


    plt.imshow(U)
    plt.show()

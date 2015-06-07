#=
    Wavefield Modelling in 10 lines

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
Now I need to make a discrete La Place operator. This will look like computational
voodoo, but keeping things in matrix form simplifies the problem. Start with
the 1-D discrete La Place equation, and notice that this expression can be
written as a matrix-vector product. Using kronecker products, this is expanded
into a 2D Laplacian operator.
=#    
    # Make the LaPlacian operator
    Dx = spdiagm(tuple(ones(n-1), -ones(n)*2, ones(n-1) ), [1,0,-1])
    Dx = kron(speye(n), Dx) / dn^2
    
    Dy = spdiagm(tuple(ones(n*n-n), -ones(n*n)*2, ones(n*n-n) ),
             [n,0,-n]) / dn^2
    LP = -(Dx + Dy)
    

#=
Make a source wavefield consisting of a point source in the middle of the model.
=#
    # Put in a single source in the center of the model
    q = zeros(size(c))
    q[n/2,n/2] = 1.0
#=
Since matrices operate on vectors, we need to flatten our models into vectors.
=# 

    # Vectorize the model
    m = (1 ./ c)[:]
    q = q[:]

#=
We can now solve the Helmholtz equation and see our modelled wavefield.
=#
    
    # Make the Helmholtz operator
    H = LP + spdiagm((w ./ m) .* 2)

    # Hu = q solve for the wavefield
    u = H \q

using PyCall
@pyimport matplotlib.pyplot as plt

    # reshape the wavefield vector
    u = reshape(u, n,n)

    plt.imshow(u)
    plt.show()

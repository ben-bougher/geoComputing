###
Am=d

<<<<<<< HEAD
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
    
=======
In 52 things about geophysics, Bryan Russel wrote an excellent article
on Am=d. This simple equation is at the core of most seismic 
modelling and processing. This article demonstrates this relationship 
by forward modelling a seismic experiment using the form Am=d.

First and for most, I will put some clear definitions on what these
algebraic symbols represent. The "vector" m is the model of the earth,
which acts as inputs of the experiment. It is very counter-intuitive 
to think of earth model as 1D vector instead of a volume. We no longer
care about the physical space of the model, so we flatten the 
model into a vector of parameters. The matrix A describes the 
mechanics of our forward model, and transforms our earth model into a 
vector of observed data d.

I will do a simple forward model of three reflectors in constant
velocity medium. In this case m will the reflectivity of the earth.
###

m = zeros(1000); m[300,600,750] = 1.0
>>>>>>> 9e5abdcf6db90758243c48b9381921e4919cafb2

###
We will now build our operator A through a series of linear operators.
Convolution is a linear operation, meaning it can
be written as matrix vector product. I will do this explicitly using
a circulant matrix.
###
w = Ricker(10.0)
C = Circulant(w)

###
We want to make more than one measurement for our model, so we 
create a fold using another linear operator.
###
F = kron(speye(n),fold)

<<<<<<< HEAD
#=
We can now solve the Helmholtz equation and see our modelled wavefield. Note
that solving this discrete PDE is solving a system a of linear equations. 
=#
    
    # Make the Helmholtz operator
    H = LP + spdiagm((w ./ m) .* 2)
=======
###
>>>>>>> 9e5abdcf6db90758243c48b9381921e4919cafb2

In practice we create fold through offset, which because of the
geometry results in normal move out (NMO). Since NMO is just 
translating elements, we can also write it as a linear operator.
###

<<<<<<< HEAD
    # reshape the factor back to our model dimensions
    U = reshape(u,n,n)

using PyCall
@pyimport matplotlib.pyplot as plt


    plt.imshow(U)
    plt.show()
=======
N = opNMO(n, velocity)

###
We can cascade these operators to define our final opertor A:
###

A = N*F*C

### 
We can now generate our synthetic data using Am=d.
###
d = A*m;

###
We can reshape d into a gather and plot it to see the hyperbolic
arrivals of seismic wavelets.
###


>>>>>>> 9e5abdcf6db90758243c48b9381921e4919cafb2

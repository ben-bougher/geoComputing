#=
Am=d

In 52 things about geophysics, Bryan Russel wrote an excellent article
on Am=d. This simple equation is at the core of most seismic 
modelling and processing. This article demonstrates this relationship 
by forward modelling a seismic experiment using the linear Am=d.

First and for most, I will put some clear definitions on what these
algebraic symbols represent. The "vector" m is the model of the earth,
which acts as input parameters of the experiment. It is very counter-intuitive 
to think of earth model as 1D vector instead of a volume. It helps to think of
the earth model as a list of inputs into our equation, not as a geometrical
representation of the earth. The matrix A describes the 
mechanics of our forward model, and transforms our vector of earth parameters m
into a vector of observed data d.

This article demonstrates a simple forward model of three reflectors in
a constant velocity medium. In this case m will the reflectivity of the earth,
which for simplicity we will keep in the time domain.
=#
using SpecialMatrices
using PyCall
include("helper.jl")
    
@pyimport matplotlib.pyplot as plt
    
# Three reflectors
n = 1001
m = zeros(n)
m[300] = 1;
m[550] = 1;
m[750] = 1;
    
#=
We will now build our operator A through a series of linear operators. We start
by modelling the seismic experiment as a convloution between the reflectivity
and a wavelet. Convolution is a linear operation, meaning it can
be written as matrix vector product. I will do this explicitly using
a Toeplitz matrix, where each row is shifted version of the previous row.
=#

dt, dx, f = .001, 10.0, 10
duration = dt * (size(m,1)-1) * 2
    
w = Ricker(duration, dt, f)
C = full(Toeplitz(w))

#=
We want to make more than one measurement for our model, so we 
create a fold using another linear operator. 
=#
fold = 200 
F = kron(ones(fold),speye(n))

#=

In practice we create fold through offset, which because of the
geometry results in normal move out (NMO). Since NMO is just 
translating elements, we can also write it as a linear operator.
=#

N = opNMO(n, fold, velocity, dt, dx)

#=
We can cascade these operators to define our final opertor A:
=#

A = N*F*C

#= 
We can now generate our synthetic data using Am=d.
=#
    
d = A*m;

D = reshape(d, fold, size(m,1))'
#=
We can reshape d into a gather and plot it to see the hyperbolic
arrivals of seismic wavelets.
=#



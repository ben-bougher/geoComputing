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
into a vector of observed data d. Generally, m is the physical properties of the
earth and A is the physics that models the data.

This article demonstrates a simple forward model of an angle gather typically 
used in AVO analysis. In this case m will the elastic properties of the earth,
which for simplicity we will define in the time domain.

=#
using SpecialMatrices
using PyCall
include("helper.jl")
    
@pyimport matplotlib.pyplot as plt
    
#=
Build a 1D model with vp, vs, and rho.
=#  
n = 1001
vp, vs, rho = zeros(n), zeros(n), zeros(n)
vp[0:300] = 1500; vs[0:300] = 0; rho[0:300] = 1000.0;
vp[300:750] = 3500; vs[300:750] = 2000; rho[300:750] = 3000.0;
vp[750:end] = 5000; vs[750:end] = 2500; rho[750:end] = 3500.0;

# Add 10% some noise to the earth model
vp += randn(n)*.1*mean(vp)
vs += randn(n)*.1*mean(vs)
rho += randn(n)*.1*mean(rho)

#=
Using the aki-richards equation, the angle-dependent p-wave reflectivity
can be modelled by the linear equation:
    aki richards
Defining a gradient operator and lumping the constant terms the aki-richards
equation takes on a matrix form
=#

#=
Seismic reflection data can be modelled as a convolution of a wavelet with a
reflectivity series. Convolution is a linear operation, so by definition it
can be expressed in matrix form. A convolution operator has a special form
called a Toeplitz matrix, which I form with a 10 Hz Ricker wavelet.
=#
dt, f = .001, 10.0
duration = dt * n * 2
    
w = Ricker(duration, dt, f)
C = full(Toeplitz(w))




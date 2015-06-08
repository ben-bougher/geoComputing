#=
Am=d

In 52 things about geophysics, Bryan Russel wrote an excellent article
on Am=d, the simple matrix-vector product at the core of seismic 
modelling and processing. This article demonstrates this relationship 
by forward modelling a seismic experiment using the linear model Am=d.

First and for most, let's put clear definitions on what these
algebraic symbols represent. The "vector" m is the model of the earth,
which acts as the input parameters of the experiment. If it seems counter-intuitive 
to represent the earth as a 1D vector instead of a slice or volume, think of
the earth model as a list of inputs into our experiment, not as a geospatial
representation. The matrix A describes the mechanics of our 
forward model, and transforms our vector of earth parameters m
into a vector of observed data d. Generally, m is the physical 
properties of the earth and A is the physics that models data.

This article demonstrates a simple forward model of an angle gather typically 
used in AVO analysis. In this case m is the elastic properties of the earth,
which for simplicity we will define in the time domain.

Using the aki-richards equation, the angle-dependent p-wave reflectivity
can be modelled by the linear equation:
    aki richards
where blah blah blah

Defining a gradient operator matrix
 G matrix
and lumping the constant terms the aki-richards
equation can be written as a matrix vector product.

(aki richards equation)

We model the physics of seismic reflection data as a convolution of a wavelet with a
reflectivity series. Convolution is a linear operation, so by definition it
can also be expressed in matrix form. This operator is a special type of matrix,
called a Toeplitz matrix, where each row is a shifted copy of the previous row.

(toeplitz example)

Combining the Aki-Richards operator with the convolution operator we
get the forward modelling operator A. 
A = C*R
Our angle gather can now be synthesized by the
matrix vector product Am=d. 
=#

# helper functions available on github
include("helper.jl")
    
# Build a 1D 2-layer earth model
n = 250
vp, vs, rho = zeros(n+1), zeros(n+1), zeros(n+1)
vp[1:n/2] = 2750; vs[1:n/2] = 1600; rho[1:n/2] = 2150.0;
vp[n/2:end] = 3000; vs[n/2:end] = 2200; rho[n/2:end] = 2500.0;

# stack the parameters into a flat vector
m = [rho, vp, vs]

# Make the Aki Richards reflectivity operator for angles [0-40 deg]
theta = [0:2:40]
R = [opAkiRichards(n, angle, mean(vp), mean(vs), mean(rho))
     for angle in theta]
R = cat(1, R...)

# Make a 40 Hz Ricker wavelet
dt, f = .001, 40.0
duration = dt * (n-1) * 2  
w = Ricker(duration, dt, f)

# Make a single Toeplitz convolution matrix
C0 = opToeplitz(w)

# Expand the convolution matrix to match the offset dimensions
# using a kronecker matrix product. 
C =  kron(speye(size(theta,1)), C0)

# Make the forward modelling operator A
A = C*R

# forward model the data
d = A*m

# reshape the data vector into the proper dimensions
d = reshape(d, n, size(theta,1))

# Plot the gather
plot(d)





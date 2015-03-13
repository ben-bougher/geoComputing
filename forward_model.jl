###
Am=d

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

###

In practice we create fold through offset, which because of the
geometry results in normal move out (NMO). Since NMO is just 
translating elements, we can also write it as a linear operator.
###

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



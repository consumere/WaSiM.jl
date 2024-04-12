using Distributions

# Define a Weibull distribution with shape parameter a and 
#scale parameter b
#First, the shape parameter provides information about the shape 
#of the distribution. Second, the scale parameter describes 
#how stretched or shrunk the distribution should be.

a = 1.5
b = 2.0
wbl = Weibull(a, b)

# Get the parameters of the Weibull distribution
params(wbl)  # This will return a tuple of the shape and scale parameters (a, b)
# Get the shape parameter of the Weibull distribution
shape(wbl)   # This will return the shape parameter (a)
# Get the scale parameter of the Weibull distribution
scale(wbl)   # This will return the scale parameter (b)


#fn = "D:/Wasim/Goldbach/meteo/q15min_goldbach.txt"
#da = xddf(fn)
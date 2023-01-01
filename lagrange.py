# find the Sun-Earth Lagrange point L1, using Newton's method

from numpy import arange, zeros, ones, linspace
from pylab import plot, show, xlim, ylim, xlabel, ylabel

# physical parameters
G = 6.674e-11 # Newton's gravitational constant in SI units (m^3 kg^-1 s^-2_
Me = 5.974e24 # Mass of the Earth in SI units (kg) 
Ms = 1.989e30 # Mass of the Sun in SI units (kg)

R = 1.496e11  # 1AU = mean Earth-Sun distance (m)
w = 1.99e-7 # Earth's orbital frequency (rad s^-1)

#the function whose root we want to find (force along the Sun-Earth line)
def f(x):
    return G*Ms/x/x - G*Me/(R-x)/(R-x) - w*w*x

#finding the derivative of this function we get:

#-2*G*Ms/x^3 - 2MeG/(R-X)^3 - w^2

def f1(x):
    return -2*G*Ms/x**3 - 2*Me*G/(R-x)**3 - w**2

#x = linspace(1, 1.4e11, 100)
rs = arange(1e6, 1e12, 1.0e5)
ef = zeros(len(rs))
y =  []


def Newton(x):
    #define our estimate as per newtons formula.
    x2 = x - f(x)/f1(x)

    while(abs((x2 - x)) > 0.00000001):
#Changing our guess to our estimate thereby changing the next estimate until the error in the location of the root is small.
        x = x2
        x2 = x - f(x)/f1(x)
    return (x2)

print(Newton(1.4e11))

#plotting to find a good estimate for newtons method
for i in range(len(rs)):
    r = rs[i]
    ef[i] = f(r)
    
plot(rs, ef)
ylim(-1,1)
show()

#Give it some time it works :)

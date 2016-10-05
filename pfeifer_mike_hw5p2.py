# import commands to be used later in code
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# Define the different functions that we will integrate.
def fun1(x):
    val = x**5
    return val
def fun2(x):
     val = 1/(x**2 + 1)
     return val
def fun3(x):
    val = (x*np.sin(x))/(1 + np.cos(x)*np.cos(x))
    return val

# Write code that will integrate the functions above by midpoing approximation
def midpoint(func, a, b, n):
    ans = 0
    step = (b-a)/n
# Start at first mid point
    a = a + step/2
# End integration at midpoing before bound b
# Iteratively sum the value at the midpoints multiplies by the step size
    while a <= b:
        value = func(a)
        z = step*value
        a = a + step
        ans = z + ans     
    return ans
    
# Write code that will integrate the functions above by simpsons approximation
def simpsons(func, a, b, n):
    step = (b-a)/n
# The multiplication of the calculated function values will alternate between
# 2 and 4. This count will allow us to alternate with an if statement
    count = 1
    f1 = 0
    f2 = 0
# set constants and solve at the end points
    f0 = func(a)
    f3 = func(b)
# take first step to start the simpsons rule
    x = a + step
    while x < b:
# if count is even then multiply by two
        if count%2 == 0 :
            f1 = 2 * func(x) + f1
# if count is odd then multiply by four
        else :
            f2 = 4 * func(x) + f2
        count = count + 1
        x += step
# Solve using the simpsons rule
    ans = (step/3)*(f0 + f1 + f2 + f3)
    return ans

# Write code that will integrate the functions above by gauss legendre approximation
def gauss(func, a, b, n):
    ans = 0
# Find the roots of the polynomial using legendre function
    x,w = np.polynomial.legendre.leggauss(n)
# Plug in values of x so we may find f(t)
    t = 0.5*(x+1)*(b-a)+a
# sum the results and multiply by step size
    ans = sum(w * func(t))* (0.5*(b-a))
    return ans

# Graph the error on a log log plot for the first function for all integration methods
errmid = [abs(midpoint(fun1,-1,1,i)) for i in range (1,101)]
errsim = [abs(simpsons(fun1,-1,1,i)) for i in range (1,101)]
errgau = [abs(gauss(fun1,-1,1,i)) for i in range (1,101)]

plt.Figure
plt.figure(figsize = (14,9))
plt.yscale('log')
plt.xscale('log')
plt.xlabel("n")
plt.ylabel("Error")
plt.title ("Error in Function 1")
plt.xlim(1, 100)
plt.ylim(-.1, 1 )

plt.plot(range (1,101), errmid, 'k-', range (1,101), errsim, 'r-', range (1,101), errgau, 'b-')
# It is found that the simpsons and gauss legendre methods have no error when calculating the first function

# Graph the error on a log log plot for the second function for all integration methods
errmid = [abs(2*np.arctan(5) - midpoint(fun2,-5,5,i)) for i in range (1,101)]
errsim = [abs(2*np.arctan(5) - simpsons(fun2,-5,5,i)) for i in range (1,101)]
errgau = [abs(2*np.arctan(5) - gauss(fun2,-5,5,i)) for i in range (1,101)]

plt.Figure
plt.figure(figsize = (14,9))
plt.yscale('log')
plt.xscale('log')
plt.xlabel("n")
plt.ylabel("Error")
plt.title ("Error in Function 2")
plt.xlim(1, 100)
plt.ylim(-.1, 10 )

plt.plot(range (1,101), errmid, 'k-', range (1,101), errsim, 'r-', range (1,101), errgau, 'b-')

# Graph the error on a log log plot for the third function for all integration methods
errmid = [abs((np.pi**2/4)-midpoint(fun3,0,np.pi,i)) for i in range (1,101)]
errsim = [abs((np.pi**2/4)-simpsons(fun3,0,np.pi,i)) for i in range (1,101)]
errgau = [abs((np.pi**2/4)-gauss(fun3,0,np.pi,i)) for i in range (1,101)]

plt.Figure
plt.figure(figsize = (14,9))
plt.yscale('log')
plt.xscale('log')
plt.xlabel("n")
plt.ylabel("Error")
plt.title ("Error in Function 3")
plt.xlim(1, 100)
plt.ylim(-.1, 10 )

plt.plot(range (1,101), errmid, 'k-', range (1,101), errsim, 'r-', range (1,101), errgau, 'b-')



    
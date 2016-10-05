# import commands to be used later in code
from __future__ import division
from fractions import Fraction
from scipy.optimize import minimize
import re
import matplotlib.pyplot as plt
import numpy as np

# open the file for atomic data
atom = open('atomic_masses.txt', 'r')
# rest mass of neutron and hydrogen atom
mn = 1.0086649
mh = 1.00782503207

# search the text file for the data to be used
num = r'Atomic Number = .*'
relmass = r'Relative Atomic Mass .*'
massnum = r'Mass Number = .*'

# form a matrix to extract the data for each element
m = [num, relmass,massnum]
# make some empty lists/matricies that will be used 
give = []
binding = []
ex=[]

# use a for loop to extract the data for each isotope
for i in range(0,3):
    catalog = m[i]
    atm = open('atomic_masses.txt', 'r')
    elem = []
    for line in atm:
        s = str(line)
        n = re.search(catalog, s)
        if n:
           res = n.group(0).split(' = ')
           rest = float(res[1].translate(None,"()[]#")) # Filter out any unwanted data or symbols
           elem.append(rest)
    give.append(elem)

# assign variables to each of the lists of atomic data
for z, ma, a in map(None,give[0],give[1],give[2]):
    bind = (-ma + z*mh + (a-z)*mn)* 931.494
    bpn = bind/ma  #normalize the binding energy per nucleon
# binding energy for each isotope    
    binding.append(bpn)
    
# Form a matrix A of the atomic data and the solution matrix b to find coefficients
b = np.array([binding]).T
# Values of x that will be used to form matrix A
x1 = np.asarray(give[2])
x2 = np.asarray([a**Fraction('2/3') for a in give [2]])
x3 = np.asarray([(z**2)/(a**Fraction('1/3')) for a,z in map(None,give [2],give[0])])
x4 = np.asarray([((a-2*z)**2)/a for a,z in map(None,give [2],give[0])])
x5 = np.asarray([((-1)**z%2 + (-1)**((a-z)%2))/(a**Fraction('1/2')) for a,z in map(None,give [2],give[0])])

A = np.column_stack((x1,x2,x3,x4,x5))
# using the transpose and origional matricies we may find a solution c for the coefficients
AtA = (A.T).dot(A)

Atb = (A.T).dot(b)

c = np.linalg.solve(AtA, Atb)
# Plug the values of c back in to the equation to solve by least-squares approximation
lstsqr = c[0]*x1 + c[1]*x2 + c[2]*x3 + c[3]*x4 + c[4]*x5

# Format the graph that will be used later (this is here so I could check my work on the least-squares approx.)
plt.figure(figsize = (14,9))
plt.xlabel("Atomic Number (Z)")
plt.ylabel("Binding Energy (MeV)")
plt.xlim(0, 300)
plt.ylim(0, 10)
# print the least-squares coefficients and errors
print "The coefficients for the least-squares approximation are given as" 
print c.T
print "The maximum and average error of the least-squares approximation are"
eq = (c[0]*x1 + c[1]*x2 + c[2]*x3 + c[3]*x4 + c[4]*x5)
err1 = max([abs(float(a - h)) for a,h in map(None, b, eq)])
err3 = sum([abs(float(a - h)) for a,h in map(None, b, eq)])/3178
print "%f and %f" % (err1,err3)
print "respectively"

# Solve by max-min approximation
# Finding the error for the max-min approach
def maxmin(n):
# Finding average error
    eq = (x1*n[0] + x2* n[1] + x3*n[2] + x4* n[3] + x5*n[4])
    err2 = [float(a - c) for a,c in map(None, b, eq)]
    return np.linalg.norm(err2, ord = 1)
def merr(n):
# Finding maximum error
    eq = (x1*n[0] + x2* n[1] + x3*n[2] + x4* n[3] + x5*n[4])
    err2 = [abs(float(a - c)) for a,c in map(None, b, eq)]
    return max(err2)
# minimize the total error to obtain a solution 
x0 = np.array([float(a) for a in c])
sol = minimize(maxmin,c)
sol2 = sol.x
mxmn = sol2[0]*x1 + sol2[1]*x2 + sol2[2]*x3 + sol2[3]*x4 + sol2[4]*x5
err2 = (merr(sol2))
err4 = (maxmin(sol2))/3178
# Print the solution and the errors
print "\n"
print "The coefficients for the max-min approximation are given as" 
print sol2
print "The maximum and average error of the least squares approximation are"
print "%f and %f" % (err2,err4)
print "respectively"

# Plot all of the data
plt.plot(give[1], binding, 'k.', give[1], lstsqr, 'r.', give[1], mxmn, 'b.')
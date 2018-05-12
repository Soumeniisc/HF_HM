#https://matplotlib.org/examples/color/colormaps_reference.html
#here we are calculating Ginzburg landau coeffiecents for IHM next nearest neighbour hopping on square lattice at half filling , using formula mentioned in soumen et al, PRB2015. 
#created on 07/05/18
#TODO: in the k integration on should fix the chemical potential. Here I have used wrong expression. please look at the .cpp

from __future__ import division
import matplotlib.pyplot as plt
from scipy import *
import numpy as np
import scipy.linalg as sla
from matplotlib import rc
rc('text', usetex=True)

rc('xtick', labelsize=20)
rc('ytick', labelsize=20)

N=500

t = 0.5
delta = 0.5
t2 = 0.2
dim = 2.0
cos_ = np.zeros(2*N)

for i in range(2*N):
	kx = np.pi*(i-N)/N + np.pi/(2*N)
	cos_[i] = np.cos(kx) 


def function(X,U):
	ms = X
	g = ms*U/2.
	mu = U/2.0
	energy = 0.0
	counter = 0.0
	for i in range(N):
		for j in range(i+1):
			if (j!=i):
				weight = 8.0
			else:
				weight = 4.0

			E_k =  np.sqrt((2*t*(cos_[i] + cos_[j]))**2 + g**2) 
			lambdan = -E_k + U/2.0 -mu
			lambdap = +E_k + U/2.0 -mu
			if (lambdan < 0.0):
				energy = energy + weight/E_k
				counter = counter + weight*1.
			if (lambdap < 0.0):
				energy = energy - weight/E_k
				counter = counter + weight*1.
	#print counter/(2*N)**2
	return U*energy/(2*N)**2 - 2.0

def broyden_good(x,U, f_equations, J_equations, tol=10e-10, maxIters=50):
    steps_taken = 0
    #print "steps_taken,ms",steps_taken,x
    f = f_equations(x,U)
    J = J_equations(x,U)
    while abs(f) > tol and steps_taken < maxIters: # norm calculate length of a vector
 	print "steps_taken:", steps_taken, " ms: ", x, "J:", J , "f:",f, "s=-f/J:", -f/J
        #s = sla.solve(J,-1*f)  # numpy.linalg.solve(a, b). It does a*x = b and return x.
	s = -f/J
 	alpha = 0.4
	s = alpha*s
	
	#temporary_increament
	xp =  x + s
	newf = f_equations(xp,U)
	look_alpha = 1
	while(abs(newf)>abs(f) and look_alpha< 20):
		alpha = alpha/2.0
		print "alpha,s,----------------------------",alpha
		s = alpha*s
		xp = x + s # notice that his alpha paramter is very important otherwise it could go to infintity.
        	newf = f_equations(xp,U)
		look_alpha = look_alpha + 1
	""" 
	alpha=-alpha_
	look_alpha = 1
	while(abs(newf)>abs(f) and look_alpha< 20):
		alpha = alpha/2.0
		print "alpha,s,----------------------------",alpha
		xp = x + alpha*s # notice that his alpha paramter is very important otherwise it could go to infintity.
        	newf = f_equations(xp,U)
		look_alpha = look_alpha + 1
	"""
        x = x + s # notice that his alpha paramter is very important otherwise it could go to infintity.
        newf = f_equations(x,U)
        z = newf - f
 	
	# Updating J by Broyden's approximation method.
        J = J + ((z - J*s)*s) / s**2
 
        f = newf
        steps_taken += 1
     
    return steps_taken, x  


#def Js(x,y):
#    return np.array([[1,2],
#             [2, 16]])

def Js(x,U):
    dx = 0.01
    fx = function(x,U)
    fxdx = function(x+dx,U)
    if (fxdx > fx):
	Jack = (fx-fxdx)/dx
    else:
	Jack = -(fx-fxdx)/dx
    
    return Jack  # CHOOSING THE INTIAL JACOBIAN OR SLOP IS VERY IMPORTANT OTHERWISE SOLUTION WILL BLOWUP AT FIRST ITERATION.


U_list = [0.08,0.09,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0]
#U_list = [0.2,0.3,1.0]
#U_list = [0.06,0.07]
data = open("AFM_t20_HM_python_broyden.dat","w")
print >> data, "#U, mu, x, no of iteration, energy condition"
x0 = 0.00298635064876
for l,U in enumerate(U_list):
	tol = 10.0** -11
	print "U:----------------------------",U
	maxIters = 10000
	#x0 = np.array([0.1])
	
	mu = U/2.
	n,x = broyden_good(x0, U, function, Js, tol, maxIters )
	f = function(x,U)
	x0 = x
	print "iterations: ", n, "ms:",x
	print >> data, U, mu, x, n, f


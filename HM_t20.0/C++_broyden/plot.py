import numpy as np
import matplotlib.pyplot as plt
from numpy import *

t2 = 0
t = 0.5
SCF = loadtxt("delta0t2%s.txt"%t2)
broyden = loadtxt("broyden_AFM2.dat")


plt.plot(np.array(SCF[:,0])/t,np.array(SCF[:,6]),'-o',label="SCF")
plt.plot(np.array(broyden[:,0])/t,abs(np.array(broyden[:,2])),'-o',label="broyden")
plt.legend()
plt.savefig("ms_U_t2%s.eps"%t2, format="eps")
plt.show()

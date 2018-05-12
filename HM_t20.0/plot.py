import numpy as np
import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rc
#rc('text', usetex=True)

#rc('xtick', labelsize=20)
#rc('ytick', labelsize=20)

t2 = 0
t = 0.5
SCF = loadtxt("AFM_t20_HM_CPP_self_consistency.txt")
cppbroyden = loadtxt("AFM_t20_HM_CPP_broyden.dat")
python_broyden = loadtxt("AFM_t20_HM_python_broyden.dat")

plt.title("HF  ms of Hubbard Model as a function U from different solver")
plt.plot(np.array(SCF[:,0])/t,np.array(SCF[:,6]),'-o',label="cppSCF")
plt.plot(np.array(cppbroyden[:,0])/t,abs(np.array(cppbroyden[:,2])),'-o',label="cppbroyden")
plt.plot(np.array(python_broyden[:,0])/t,abs(np.array(python_broyden[:,2])),'-*',label="python_broyden")
plt.legend()
plt.xlabel("U/t")
plt.ylabel("ms")
plt.savefig("ms_U_t2%s.eps"%t2, format="eps")
plt.show()

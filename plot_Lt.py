import gasgiant as gg
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def derivs(f,x,M):
	S = f[0]
	t = 10.0**x
	return -2.303*t*S/tS(M,S)

def cooling_curve(M,S_0):
	t = np.linspace(2.0,10.0,100)
	results = odeint(derivs,S_0,t,args=(M,))
	lum = [LL(M,S)/gg.Lsun for S in results[:,0]]
	return 10.0**t,lum

# read in the models
models = gg.read_planet_models()
LL = models[0]
tS = models[1]

# now plot some curves
# same as Fig 4 lower panel in Marleau et al. 2014
t,L = cooling_curve(1.0,11.45)
plt.plot(t,L)
t,L = cooling_curve(3.0,11.45)
plt.plot(t,L)
t,L = cooling_curve(10.0,11.45)
plt.plot(t,L)

t,L = cooling_curve(1.0,9.45)
plt.plot(t,L)
t,L = cooling_curve(3.0,9.45)
plt.plot(t,L)
t,L = cooling_curve(10.0,9.45)
plt.plot(t,L)

# plotting stuff
plt.xlabel(r'$\mathrm{Time\ (yr)}$')
plt.ylabel(r'$\mathrm{Luminosity\ (L_\odot)}$')
plt.xlim((6e5,3e9))
plt.ylim((1e-8,1e-3))
plt.xscale('log')
plt.yscale('log')
plt.show()

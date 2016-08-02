from __future__ import print_function
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import sys
import gasgiant as gg
import seaborn as sns
#sns.set_style("ticks",{'xtick.direction': u'in','ytick.direction': u'in'})
#sns.set_context("paper",font_scale=1.1)
fig = plt.figure(figsize=(8.5,5))
plt.rc('text', usetex=False)
plt.rc('font', family='serif')

def get_S(L,Starget,Ttop,Ptop,Mdot):
	x,x,x,x,x,x,x,Sconv,x,x,x = gg.make_envelope(L*4e33,Ptop,Ttop,Mdot)
	print('-----> ',L,Sconv,Starget)
	return Sconv-Starget

lum,tS, Teff, rads = gg.read_planet_models()


#plt.style.use('seaborn-white')
ax = fig.add_subplot(3,3,1)

if len(sys.argv) < 5:
	print("Arguments: L(LSun) Mdot(ME/yr) Ptop(bars) Ttop(K) plot_flag mass radius Starget")
	exit()

Lin = float(sys.argv[1])*4e33
Mdot = float(sys.argv[2])*2e20      # g/s corresponds to 1e-2 MEarth/yr
Ptop = float(sys.argv[3])*1e6
Ttop = float(sys.argv[4])
plot_flag = int(sys.argv[5])
mass = float(sys.argv[6])
radius = float(sys.argv[7])
Starget = float(sys.argv[8])

gg.radius = 7e9 * radius
gg.mass = 2e30 * mass
gg.grav = 6.67e-8 * gg.mass/gg.radius**2

gg.kappa_min = 0.0

#Mdot = 0.01 *2e20
#Ptop = 0.1 * 1e6
#plot_flag = 0
#L = 6.2e-5*4e33
#Mdot = 0.01 *2e20
#Ptop = 0.1 * 1e6
#Ttop =1200.0

#Temps=[2000.0,1800.0,1600.0,1400.0,1200.0,1000.0]

#gg.tau_cutoff=2.0/3.0
#gg.tau_cutoff = -1.0

#for Ttop in Temps:
#Starget = 10.5
if Ltop<0.0:
	Ltop = 4e33*brentq(get_S,1e-5,1e-3,rtol=1e-5,args=(Starget,Ttop,Ptop,Mdot))
else:
	Ltop=Lin
P, T, Larr, dell, Sarr, Prcb, Trcb, Sconv,Lrcb,tau,deladarr = gg.make_envelope(Ltop,Ptop,Ttop,Mdot)
Ttop = T[0]
Ptop = P[0]
Lrcb=Larr[-1]
print("Ttop: L, Srcb, Prcb, Trcb, Lrcb, kappa_rcb = ",Ttop,Ltop/4e33,Sconv,Prcb,Trcb,Lrcb/4e33,gg.opacity(Prcb,Trcb))


i=0
while tau[i]<2.0/3.0 and i<len(tau)-1:
	i+=1
print("At the photosphere: T=", T[i], " P=",P[i]);
print("Teff = ", (Larr[0]/(4.0*3.1415*gg.radius**2*gg.sigmaSB))**0.25)


#if not plot_flag:
#	f = open('models.dat', 'a')
#	f.write("%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n" % (Mdot/2e20,Ptop/1e6,Ttop,Ltop/4e33,Sconv,Prcb,Trcb,Lrcb/4e33,gg.opacity(Prcb,Trcb),gg.entropy(Ptop,Ttop),gg.mass,gg.radius))
#	f.close()


print("---------------------------------------------------")
print("top entropy = ",gg.entropy(Ptop,Ttop))
print("for S=",Starget," L=",lum(mass,Starget)/3.83e33, " tS=", tS(mass,Starget))
print("with accretion, cooling time = ", lum(mass,Starget)*tS(mass,Starget)/Lrcb)
print("cooling time/accretion time = ", (lum(mass,Starget)*tS(mass,Starget)/Lrcb) / (10.0*318.0*2e20/Mdot))
print("---------------------------------------------------")


if plot_flag:

	f = open('profile.dat','w')
	for i in range(len(P)):
		f.write("%lg %lg %lg %lg %lg %lg %lg %lg\n" % (P[i],T[i],Larr[i],Sarr[i],gg.opacity(P[i],T[i]),dell[i],gg.delad(P[i],T[i]),tau[i]))
	f.close()

	plt.plot(P,T,'k')
	plt.plot(P[dell>=deladarr],T[dell>=deladarr],'ko')
	
	ax.set_yscale('log')
	ax.set_xscale('log')
	plt.xlabel('$P\ (\mathrm{erg\ cm^{-3}})$')
	plt.ylabel('$T\ (\mathrm{K})$')

	ax = fig.add_subplot(3,3,2)
	ax.set_xscale('log')
	#plt.ylim([8e28,1e30])
	plt.plot(P,Larr/4e33,'k')
	#plt.plot(P2,Larr2,'r--')
	ax.set_yscale('linear')
	plt.xlabel('$P\ (\mathrm{erg\ cm^{-3}})$')
	plt.ylabel('$L\ (cgs)$')
	#plt.ylabel('$L\ (L_{\odot})$')
	#plt.savefig('LS.pdf')

	ax = fig.add_subplot(3,3,3)
	ax.set_xscale('log')
	plt.plot(P,Sarr,'k')
	#plt.plot(P2,Sarr2,'r--')
	ax.set_yscale('linear')
	plt.xlabel('$P\ (\mathrm{erg\ cm^{-3}})$')
	plt.ylabel('$S\ (k_b/m_p)$')
	#plt.ylabel('$L\ (L_{\odot})$')


	ax = fig.add_subplot(3,3,4)
	ax.set_xscale('log')
	kappa = [gg.opacity(P1,T1) for P1,T1 in zip(P,T)]
	plt.plot(P,kappa,'k')
	#kappa2 = [opacity(P0,T0) for P0,T0 in zip(P2,T2)]
	#plt.plot(P2,kappa2,'r')
	ax.set_yscale('log')
	plt.xlabel('$P\ (\mathrm{erg\ cm^{-3}})$')
	plt.ylabel('$\kappa (\mathrm{cm^2\ g^{-1}})$')

	ax = fig.add_subplot(3,3,5)
	ax.set_xscale('log')
	delad = [gg.delad(P1,T1) for P1,T1 in zip(P,T)]
	delrad = [gg.calculate_dell_rad(P1,T1,L1) for P1,T1,L1 in zip(P,T,Larr)]
	plt.plot(P,dell-deladarr,'k')
	#plt.plot(P,delrad,'r')
	#kappa2 = [opacity(P0,T0) for P0,T0 in zip(P2,T2)]
	#plt.plot(P2,kappa2,'r')
	#ax.set_yscale('log')
	plt.xlabel('$P\ (\mathrm{erg\ cm^{-3}})$')
	plt.ylabel(r'$\nabla-\nabla_{\rm ad}$')

	ax = fig.add_subplot(3,3,6)
	ax.set_xscale('log')
	#cp = [gg.CP(P1,T1)/(gg.kB/gg.mp) for P1,T1 in zip(P,T)]
	plt.plot(P,tau,'k')
	#ax.set_yscale('log')
	plt.xlabel('$P\ (\mathrm{erg\ cm^{-3}})$')
	plt.ylabel(r'$\tau$')

	ax = fig.add_subplot(3,3,7)
	ax.set_xscale('log')
	X = [gg.get_XH(P1,T1) for P1,T1 in zip(P,T)]
	#rho = [gg.density(P1,T1) for P1,T1 in zip(P,T)]
	#plt.plot(P,rho,'k')
	plt.plot(P,X,'k')
	ax.set_yscale('linear')
	plt.xlabel('$P\ (\mathrm{erg\ cm^{-3}})$')
	#plt.ylabel(r'$\rho$')
	plt.ylabel(r'$X_H$')


	ax = fig.add_subplot(3,3,8)
	ax.set_xscale('log')
	cp = [gg.CP(P1,T1)/(gg.kB/gg.mp) for P1,T1 in zip(P,T)]
	plt.plot(P,cp,'k')
	ax.set_yscale('log')
	plt.xlabel('$P\ (\mathrm{erg\ cm^{-3}})$')
	plt.ylabel(r'$C_P$')


	ax = fig.add_subplot(3,3,9)
	ax.set_xscale('log')
	plt.plot(P,deladarr,'k')
	ax.set_yscale('linear')
	plt.xlabel('$P\ (\mathrm{erg\ cm^{-3}})$')
	plt.ylabel(r'$\nabla_{ad}$')

	
	#plt.tight_layout()
	plt.show()
	#plt.savefig('profile.pdf')


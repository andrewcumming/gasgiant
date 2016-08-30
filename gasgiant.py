from __future__ import print_function
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import odeint
from scipy import interpolate

def entropy(P,T):
	return entropy_mesa(P,T)
	#return full_entropy(P,T)
	#rho = density(P,T)
	# Gabriel's formula for S
	#return 14.4 + 1.0*np.log(T/1540.0)-0.44*np.log(rho/2.5e-10)

#Calculated from Eq A.23a of Gabriel's thesis
def get_XH(P,T):
	sigma = 2  # 2 for a symmetric molecule
	Trot = 85.4 #rotational temperature in kelvin
	Tvib = 6210.0
	zr = T/(sigma*Trot)*(1.0+Trot/(3.0*T))
	zv = np.exp(-Tvib/(2.0*T))/(1.0-np.exp(-Tvib/T))
	nprime = P/(kB*T)
	Tc = 52439.0
	xi = (1.489e-20)/(T**1.5)*zr*zv*np.exp(Tc/T)
	H_frac = -1.0/(2*xi*nprime) + np.sqrt(1.0/(2*xi*nprime)**2+1.0/(xi*nprime))
	return H_frac

#calculated using formulas from Gabriel's thesis
def full_entropy(P,T):
	
	#Individual Entropies per baryon in units of Kb
	#S = [A, B , C] 
	#S = A + B*log10(T/1000k) + c*log10(P/1Bar) 
	con = np.log10(np.e)
	SH = np.array([16.15,2.5/con,-1.0/con])
	SH2 = np.array([9.97,(7.0/4.0)/con,-0.5/con])
	SHe = np.array([4.56,(5.0/8.0)/con,-0.25/con])

	Tvib = 6210.0
	xvib = np.exp(-Tvib/T)
	Svib = (Tvib/T)*xvib/(1.0-xvib) - np.log(1.0-xvib)
	#Svib = 0.0
	
	SH2 = SH2 + 0.5*Svib*np.array([1,0,0])	

	#H2 number fraction (relative to all hydrogen molecules)
	XH2 = 1-get_XH(P,T) 
	
	#H2 and H mass fractions
	mf_H = (1.0-YHe)*(1.0-XH2)/(1.0+XH2)
	mf_H2 = 1.0-mf_H-YHe
	
	#average number of baryons per particle mu, eq A.28
	u = 1.0/((1.0-YHe)/(1.0+XH2)+YHe/4.0)
	
	#partial number fractions eq A.30
	x_He = YHe*u/4.0
	x_H = mf_H*u
	x_H2 = mf_H2*u/2.0
	#x_Htot = (1.0-YHe)*u/(1+XH2)

	#mixing entropy eq A.29
	eps = 1e-6
	#S_mix = (1.0/u)*(-x_Htot*np.log(max(x_Htot,eps)) - x_He*np.log(max(x_He,eps)))
	S_mix = (1.0/u)*(-x_H*np.log(max(x_H,eps)) - x_H2*np.log(max(x_H2,eps)) - x_He*np.log(max(x_He,eps)))
	
	S = mf_H*SH + mf_H2*SH2 + YHe*SHe + S_mix*np.array([1,0,0])
	
	Stot = S[0] + S[1]*np.log10(T/1000.0) + S[2]*np.log10(P*1e-6)
	
	return Stot   #+0.3

def entropy_rhoT(rho,T):
	# Gabriel's formula for S
	return 14.4 + 1.0*np.log(T/1540.0)-0.44*np.log(rho/2.5e-10)

def delad(P,T):
	f = 0.001 #fractional difference for deltas below
	#taking derivatives of entropy equation
	return -(entropy(P*(1+f),T)-entropy(P*(1-f),T))/(entropy(P,T*(1+f))-entropy(P,T*(1-f))) 
	#return 0.3056    # the value corresponding to Gabriel's formula for S

def pressure(rho,T):
	P = rho * kB * T /mp
	XH2 = 1-get_XH(P,T)
#	XH2=1.0
	P*= ((1.0-YHe)/(1.0+XH2) + YHe*0.25)
	return P

def find_pressure(rho,T,Ptarget):
	return pressure(rho,T)-Ptarget

def density(P,T):
	if T<100.0:
		T=100.0
	f1 = find_pressure(1e-10,T,P)
	f2 = find_pressure(100.0,T,P)
	if (f1<0.0 and f2<0.0) or (f1>0.0 and f2>0.0):
		print("same sign: ", T,P,f1,f2)
	rho = brentq(find_pressure,1e-10,100.0,rtol=1e-6,args=(T,P))
	return rho

def opacity(P,T):
	if T<100.0:
		T=100.0
		print("negative temperature!!")
	opac = 10.0**kappa(np.log10(T),logR(P,T))[0][0]
	if T<1600.0:
		opac = opac + kappa_min
	else:
		if T<1700.0:
			opac = opac + kappa_min*(1700.0-T)/100.0
	return opac

def logR(P,T):
	return np.log10(density(P,T))-3.0*np.log10(T)+18.0

def entropy_mesa(P,T):
	W = logW(P,T)
	if W>-2.9 or W<-20.0:
		print("W out of bounds for P,T=%lg,%lg" %(P,T))
		return full_entropy(P,T)
	S=S_mesa(W,np.log10(T))[0]/(kB/mp)
	return S[0]

def logW(P,T):
	return np.log10(P)-4.0*np.log10(T)


def read_planet_models():
	M_grid = np.array([], dtype='float64')
	R_grid = np.array([], dtype='float64')
	S_grid = np.array([], dtype='float64')
	L_grid = np.array([], dtype='float64')
	tS_grid = np.array([], dtype='float64')
	Teff_grid = np.array([], dtype='float64')
	for line in open('data/grid3S'):
		data = line.split()
		if len(data)>0:
			if data[0]!='#':   # data starts at line 10
				M_grid = np.append(M_grid,float(data[2]))
				R_grid = np.append(R_grid,float(data[3]))
				S_grid = np.append(S_grid,float(data[4]))
				L_grid = np.append(L_grid,float(data[12]))
				tS_grid = np.append(tS_grid,float(data[8]))
				Teff_grid = np.append(Teff_grid,float(data[6]))
	return interpolate.interp2d(M_grid,S_grid,L_grid,kind='linear'),\
		interpolate.interp2d(M_grid,S_grid,tS_grid,kind='linear'),\
		interpolate.interp2d(M_grid,S_grid,Teff_grid,kind='linear'),\
		interpolate.interp2d(M_grid,S_grid,R_grid,kind='linear')
		

def read_eos_file(filename):
	logT_grid = np.array([], dtype='float64')
	logW_grid = np.array([], dtype='float64')
	entropy_grid = np.array([], dtype='float64')
	count = 0
	T_flag = 1
	Svec = np.array([],dtype='float64')
	logT_grid = np.array([],dtype='float64')
	for line in open('data/'+filename):
		data = line.split()	
		if count == 4:
			logW_grid = np.append(logW_grid,float(data[0]))
		if count>6 and count<313:
			if T_flag:
				logT_grid = np.append(logT_grid,float(data[0]))
			Svec = np.append(Svec,float(data[3]))
		if count == 313:
			count = 1
			#print(Tvec[0],Svec[0])
			entropy_grid = np.append(entropy_grid,Svec)
			T_flag = 0
			Svec = np.array([],dtype='float64')
		count+=1
	nT = len(logT_grid)   # number of temperatures
	nW = len(logW_grid)   # number of log R values
	entropy_grid = entropy_grid.reshape(nW,nT)
	return logT_grid, logW_grid, entropy_grid

def read_eos():
	logT_grid, logW_grid, entropy_grid = read_eos_file('mesa-eosPT_02z80x.data')
	logT_grid, logW_grid, entropy_grid2 = read_eos_file('mesa-eosPT_02z60x.data')
	XX = 1.0-YHe
	entropy_grid = 10.0**entropy_grid
	entropy_grid2 = 10.0**entropy_grid2	
	entropy_grid = entropy_grid2 + (entropy_grid-entropy_grid2)*(XX-0.6)/0.2
	return interpolate.RectBivariateSpline(logW_grid, logT_grid, entropy_grid)

def read_opacity_file(opacity_option = 0):
	if opacity_option == 1:
		filename = 'lowT_fa05_gs98_z2m2_x70.data'
	else:
		filename = 'lowT_Freedman11_z0.02.data'
	count = 1
	logT_grid = np.array([], dtype='float64')
	logR_grid = np.array([], dtype='float64')
	kappa_grid = np.array([], dtype='float64')
	for line in open('data/'+filename):
		data = line.split()
		if count==6:
			# logR values given on line 6
			logR_grid = np.append(logR_grid,np.array([float(x) for x in data]))
		if count>7 and data!=[]:  # data starts at line 8, also reject the last line which is empty
			logT_grid = np.append(logT_grid,float(data[0]))
			kappa_grid = np.append(kappa_grid,[float(x) for x in data[1:]])		
		count+=1
	nT = len(logT_grid)   # number of temperatures
	nR = len(logR_grid)   # number of log R values
	kappa_grid = kappa_grid.reshape(nT,nR)
	return interpolate.RectBivariateSpline(logT_grid, logR_grid, kappa_grid)

def CP(P,T):
	S1 = entropy(P,T)
	T2 = T*1.001
	S2 = entropy(P,T2)
	dSdT = (S2-S1)/(T2-T)
	CP=T*dSdT*kB/mp	
	#return 1.472*kB/mp
	return CP
	
def calculate_dell_rad(P,T,flux):
	dell_rad = (P/T)* flux * 3.0 * opacity(P,T) / (4.0*arad*clight*T**3*grav)
	return dell_rad
	
def find_del_minus_delad(del_minus_delad,Ftarget,AA,BB,CC,dlad):
	Frad = AA*(dlad+del_minus_delad)
	del_minus_delad_eff = del_minus_delad  + 0.5*BB**2 - BB*(0.25*BB*BB+del_minus_delad)**0.5
	Fconv = CC*del_minus_delad_eff**1.5
	return Fconv + Frad - Ftarget
		
def calculate_dell(P,T,L,tau):

	flux = L/(4.0*np.pi*radius**2)
	dell_rad = calculate_dell_rad(P,T,flux)
	dlad = delad(P,T)
	
	if dell_rad < dlad:   # convectively stable
		dell = dell_rad
	else:    # convectively unstable
		AA = 4.0*arad*clight*T**4*grav/ (3.0*opacity(P,T)*P)
		CC = CP(P,T)*T*(P*density(P,T)/32.0)**0.5
		BB = 0.75*(AA/CC)*tau**2/(1.0+0.5*tau**2)
		del_minus_delad = brentq(find_del_minus_delad,0.0,20.0,rtol=1e-3,args=(flux,AA,BB,CC,dlad,))
		dell = dlad + del_minus_delad
		
	return dell

def derivs(f, P, Mdot):
	T = f[0]
	L = f[1]
	tau = f[2]
	if L<1e-10: 
		L=1e-10

	# entropy derivatives
	S1 = entropy(P,T)
	T2 = T*1.001
	S2 = entropy(P,T2)
	dSdT = (S2-S1)/(T2-T)

	P2 = P*1.001
	S2 = entropy(P2,T)
	dSdP = (S2-S1)/(P2-P)

	dell = calculate_dell(P,T,L,tau)
	if tau<0.0:
		dell = 0.0
	dTdP = (T/P)*dell
	dLdP = Mdot*T*(kB/mp)*(dSdP + dTdP*dSdT)
	dtaudP = opacity(P,T)/grav
							
	return dTdP, dLdP, dtaudP


def photosphere_solve(P,T):
	return P - 2.0*grav/(3.0*opacity(P,T))


def make_envelope(L, Psurf, Tsurf, Mdot):
	tau0=0.0
	if Psurf < 0.0:
		if Tsurf<0.0:
			Tsurf = (L/(4.0*np.pi*radius**2*sigmaSB))**0.25
		Psurf = brentq(photosphere_solve,1e3,1e8,rtol=1e-6,args=(Tsurf,))
		tau0=2.0/3.0
		print("Surface T,P=",Tsurf,Psurf)
	P = 10.0**np.linspace(np.log10(Psurf),8.0,400)
	#print(Tsurf,L,tau0,Mdot)
	results = odeint(derivs, [Tsurf,L,tau0], P, args=(Mdot,))#, atol=1e-10, rtol=1e-10)
	Tarr = results[:,0]
	Larr = results[:,1]
	tauarr = results[:,2]

	# locate the convective boundary
	#dell = calculate_dell

	dell = np.array([])
	deladarr = np.array([])
	Sarr = np.array([])
	for i in range(len(Tarr)):
		Sarr = np.append(Sarr, entropy(P[i],Tarr[i]))
		deladarr = np.append(deladarr, delad(P[i],Tarr[i]))
		if tauarr[i]>0.0:
			dell = np.append(dell, calculate_dell(P[i],Tarr[i],Larr[i],tauarr[i]))
		else:
			dell = np.append(dell, 0.0)
	Prad = P[dell<deladarr]
	Trad = Tarr[dell<deladarr]
	Lrad = Larr[dell<deladarr]
	if len(Prad)>0:
		Prcb = Prad[-1]
		Trcb = Trad[-1]
		Lrcb = Lrad[-1]
		Sconv = entropy(Prad[-1],Trad[-1])
	else:
		Prcb = P[-1]
		Trcb = Tarr[-1]
		Sconv = entropy(Prcb,Trcb)
		Lrcb = Larr[-1]

	Sconv = Sarr[-1]
	Lrcb = Larr[-1]
	return P, Tarr, Larr, dell, Sarr, Prcb, Trcb, Sconv, Lrcb, tauarr, deladarr


def ram_pressure(Mdot):
	return 1.9e4*(Mdot/2e18)*(mass/2e30)**0.5/(radius/7e9)**2.5

# Constants
kB = 1.38e-16
mp = 1.67e-24
arad = 7.5657e-15
clight = 3e10
sigmaSB = 5.67e-5
Lsun=3.828e33

RJ = 7.15e9
MJ = 1.898e30
MEperyear = 1.893e20

# planet properties
radius = 1.0*RJ
mass = 1.0*MJ
grav = 6.67e-8 * mass/radius**2
YHe = 0.243
#YHe = 0.40

# Read in opacity data
kappa = read_opacity_file(opacity_option=0)
kappa_min = 0.0
S_mesa = read_eos()

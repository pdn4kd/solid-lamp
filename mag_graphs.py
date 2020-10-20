import numpy as np
from matplotlib import pyplot as plt

Teff = 3000 # K
k = 1.380649e-23 # J/K, per CODATA 2018
h = 6.62607015e-34 #J*s, per CODATA 2018
c = 299792458 # m/s, has been the defined value since like 1981
σ = 5.67037442e-08 # W/m², per CODATA 2018
f0 = 3630.78 # Jy for V == 0. Derived from S(Jy) = 1e23 * 10**-(AB+48.6)/2.5 = 10**(3.56-AB/2.5) = 10**((8.9-AB)/2.5). Not used directly, but is listed for comparison purposes. Sources: Oke&Gunn 1983, Fukugita et al 1996. This may implicitly be making Vega magnitude 0 in V, rather than eg: 0.03.

area = 100/(σ*Teff**4) # area in m² of our 100 W emitter
distance = np.linspace(400e3,10000e3,961) # 400 to 10,000 km in 10 km steps

def planck(λ,T):
	B = (2*h*c**2/λ**5) / (np.exp(h*c/(λ*k*T)-1))
	return B

def bandwidth(λeff,Δλ):
	Δν = c*Δλ/(λeff**2-Δλ**2/4)
	return Δν

# λeff and Δλ in Å for both Johnson-Cousins UBVRI and SDSS ugriz, as given by Bessel 2005's article in the Annual Review of Astronomy and Astrophysics
# name	λeff	Δλ
#	U	3663	650
#	B	4361	890
#	V	5448	840
#	R	6407	1580
#	I	7980	1540
#	u'	3596	570
#	g'	4639	1280
#	r'	6122	1150
#	i'	7439	1230
#	z'	8896	1070

Johnson = [["U",1e-10*3663,1e-10*650],["B",1e-10*4361,1e-10*890],["V",1e-10*5448,1e-10*840],["R",1e-10*6407,1e-10*1580],["I",1e-10*7980,1e-10*1540]]
SDSS = [["u'",1e-10*3596,1e-10*570],["g'",1e-10*4639,1e-10*1280],["r'",1e-10*6122,1e-10*1150],["i'",1e-10*7439,1e-10*1230],["z'",1e-10*8896,1e-10*1070]]

#λmin and λmax for 50% transmission in 2MASS. Taken from https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec3_1b.html
#  band	optics and camera	system + atmosphere
#  J	1.11 - 1.40 um 	1.12 - 1.36* um
#  H	1.51 - 1.79 um 	1.51 - 1.78 um
#  Ks	2.00 - 2.31 um 	2.02 - 2.30 um 

twoMASS = [["J",1e-6*1.255,1e-6*0.29],["H",1e-6*1.65,1e-6*0.28],["Ks",1e-6*2.155,1e-6*0.155]]

fig, ax = plt.subplots(2, 2, sharex=True, sharey=False)
fig.suptitle("Johnson UBVRI Fluxes vs Distance (km)")
ax[0,0].set_prop_cycle(color=['purple', 'blue', 'green', 'red', 'black'])
ax[0,1].set_prop_cycle(color=['purple', 'blue', 'green', 'red', 'black'])
ax[1,0].set_prop_cycle(color=['purple', 'blue', 'green', 'red', 'black'])
ax[1,1].set_prop_cycle(color=['purple', 'blue', 'green', 'red', 'black'])
for band in Johnson:
	name,λeff,Δλ = band
	Δν = bandwidth(λeff,Δλ)
	xs = np.linspace(λeff-Δλ,λeff+Δλ,1000)
	ys = [planck(x,Teff) for x in xs]
	power = np.trapz(ys,xs)
	luminosity = power*area
	flux = [luminosity/(4*np.pi*d**2) for d in distance]
	fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
	Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
	ABmag = [8.9-2.5*np.log10(J) for J in Jy]
	ax[0,0].plot(distance*1e-3, flux, label=name)
	ax[0,0].set_title("Flux (W/m²)")
	ax[0,1].plot(distance*1e-3, fλ, label=name)
	ax[0,1].set_title("Fλ (W/m²/Å)")
	ax[1,0].plot(distance*1e-3, Jy, label=name)
	ax[1,0].set_title("Fν (Jy)")
	ax[1,1].plot(distance*1e-3, ABmag, label=name)
	ax[1,1].set_title("ABmag")
ax[0,0].legend()
ax[0,0].set_xlim(min(distance)*1e-3,max(distance)*1e-3)
fig.set_size_inches(8, 6.8)
fig.tight_layout()
#fig.savefig("Filters_Johnson.png", dpi=100)

fig, bx = plt.subplots(2, 2, sharex=True, sharey=False)
fig.suptitle("SDSS ugriz Fluxes vs Distance (km)")
bx[0,0].set_prop_cycle(color=['purple', 'green', 'red', 'gray', 'black'])
bx[0,1].set_prop_cycle(color=['purple', 'green', 'red', 'gray', 'black'])
bx[1,0].set_prop_cycle(color=['purple', 'green', 'red', 'gray', 'black'])
bx[1,1].set_prop_cycle(color=['purple', 'green', 'red', 'gray', 'black'])
for band in SDSS:
	name,λeff,Δλ = band
	Δν = bandwidth(λeff,Δλ)
	xs = np.linspace(λeff-Δλ,λeff+Δλ,1000)
	ys = [planck(x,Teff) for x in xs]
	power = np.trapz(ys,xs)
	luminosity = power*area
	flux = [luminosity/(4*np.pi*d**2) for d in distance]
	fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
	Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
	ABmag = [8.9-2.5*np.log10(J) for J in Jy]
	bx[0,0].plot(distance*1e-3, flux, label=name)
	bx[0,0].set_title("Flux (W/m²)")
	bx[0,1].plot(distance*1e-3, fλ, label=name)
	bx[0,1].set_title("Fλ (W/m²/Å)")
	bx[1,0].plot(distance*1e-3, Jy, label=name)
	bx[1,0].set_title("Fν (Jy)")
	bx[1,1].plot(distance*1e-3, ABmag, label=name)
	bx[1,1].set_title("ABmag")
bx[0,0].legend()
bx[0,0].set_xlim(min(distance)*1e-3,max(distance)*1e-3)
fig.set_size_inches(8, 6.8)
fig.tight_layout()
#fig.savefig("Filters_Sloan.png", dpi=100)

fig, cx = plt.subplots(2, 2, sharex=True, sharey=False)
fig.suptitle("2MASS JHK Fluxes vs Distance (km)")
cx[0,0].set_prop_cycle(color=['blue', 'black', 'red'])
cx[0,1].set_prop_cycle(color=['blue', 'black', 'red'])
cx[1,0].set_prop_cycle(color=['blue', 'black', 'red'])
cx[1,1].set_prop_cycle(color=['blue', 'black', 'red'])
for band in twoMASS:
	name,λeff,Δλ = band
	Δν = bandwidth(λeff,Δλ)
	xs = np.linspace(λeff-Δλ,λeff+Δλ,1000)
	ys = [planck(x,Teff) for x in xs]
	power = np.trapz(ys,xs)
	luminosity = power*area
	flux = [luminosity/(4*np.pi*d**2) for d in distance]
	fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
	Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
	ABmag = [8.9-2.5*np.log10(J) for J in Jy]
	cx[0,0].plot(distance*1e-3, flux, label=name)
	cx[0,0].set_title("Flux (W/m²)")
	cx[0,1].plot(distance*1e-3, fλ, label=name)
	cx[0,1].set_title("Fλ (W/m²/Å)")
	cx[1,0].plot(distance*1e-3, Jy, label=name)
	cx[1,0].set_title("Fν (Jy)")
	cx[1,1].plot(distance*1e-3, ABmag, label=name)
	cx[1,1].set_title("ABmag")
cx[0,0].legend()
cx[0,0].set_xlim(min(distance)*1e-3,max(distance)*1e-3)
fig.set_size_inches(8, 6.8)
fig.tight_layout()
#fig.savefig("Filters_2MASS.png", dpi=100)

plt.show()

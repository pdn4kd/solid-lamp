'''Flux calc specifically modified to compare with Vega'''
import numpy as np
import csv

Teff = 9602 # K
k = 1.380649e-23 # J/K, per CODATA 2018
h = 6.62607015e-34 #J*s, per CODATA 2018
c = 299792458 # m/s, has been the defined value since like 1981
σ = 5.67037442e-08 # W/m², per CODATA 2018
f0 = 3630.78 # Jy for V == 0. Derived from S(Jy) = 1e23 * 10**-(AB+48.6)/2.5 = 10**(3.56-AB/2.5) = 10**((8.9-AB)/2.5). Not used directly, but is listed for comparison purposes. Sources: Oke&Gunn 1983, Fukugita et al 1996. This may implicitly be making Vega magnitude 0 in V, rather than eg: 0.03.

area = 40.12*3.828e26/(σ*Teff**4) # area in m² of Vega, roughly (~40 L_sun)
distance = 7.68*149597870700*648000/np.pi # distance of ~7.68 pc in m (au value via IAU 2012 and 2015 resolutions)

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

with open('vega.csv', 'w', newline='') as csvfile:
	writer = csv.writer(csvfile, dialect='unix', delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	writer.writerow(['name'] + ['λeff'] + ['Δλ'] + ['power'] + ['luminosity'] + ['flux'] + ['fλ'] + ['Δν'] + ['Jy'] + ['ABmag'])
	
	for band in Johnson:
		name,λeff,Δλ = band
		xs = np.linspace(λeff-Δλ,λeff+Δλ,1000)
		ys = [planck(x,Teff) for x in xs]
		power = np.trapz(ys,xs)
		luminosity = power*area
		flux = luminosity/(4*np.pi*distance**2)
		fλ = flux/Δλ
		Δν = bandwidth (λeff,Δλ)
		Jy = flux/Δν * 1e26 # 1 Jy = 1e-26 W/m²/Hz
		ABmag = 8.9-2.5*np.log10(Jy)
		writer.writerow([name] + [λeff] + [Δλ] + [power] + [luminosity] + [flux] + [fλ] + [Δν] + [Jy] + [ABmag])
	
	for band in SDSS:
		name,λeff,Δλ = band
		xs = np.linspace(λeff-Δλ,λeff+Δλ,1000)
		ys = [planck(x,Teff) for x in xs]
		power = np.trapz(ys,xs)
		luminosity = power*area
		flux = luminosity/(4*np.pi*distance**2)
		fλ = flux/Δλ
		Δν = bandwidth (λeff,Δλ)
		Jy = flux/Δν * 1e26 # 1 Jy = 1e-26 W/m²/Hz
		ABmag = 8.9-2.5*np.log10(Jy)
		writer.writerow([name] + [λeff] + [Δλ] + [power] + [luminosity] + [flux] + [fλ] + [Δν] + [Jy] + [ABmag])
	
	for band in twoMASS:
		name,λeff,Δλ = band
		xs = np.linspace(λeff-Δλ,λeff+Δλ,1000)
		ys = [planck(x,Teff) for x in xs]
		power = np.trapz(ys,xs)
		luminosity = power*area
		flux = luminosity/(4*np.pi*distance**2)
		fλ = flux/Δλ
		Δν = bandwidth (λeff,Δλ)
		Jy = flux/Δν * 1e26 # 1 Jy = 1e-26 W/m²/Hz
		ABmag = 8.9-2.5*np.log10(Jy)
		writer.writerow([name] + [λeff] + [Δλ] + [power] + [luminosity] + [flux] + [fλ] + [Δν] + [Jy] + [ABmag])

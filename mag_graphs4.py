import numpy as np
from matplotlib import pyplot as plt

# covers UBRVI, ugriz, JHK, NGRST WFI bands in relevant cases for a 3000 K source at isotropic and limited diffraction angle for first 4 sets of graphs. Then for 532 nm / 1064 nm lasers for the second 2 sets. (for those, 100% of emitted light is assumed to be completely captured by the filters in question)

Teff = 3000 # K
k = 1.380649e-23 # J/K, per CODATA 2018
h = 6.62607015e-34 #J*s, per CODATA 2018
c = 299792458 # m/s, has been the defined value since like 1981
σ = 5.67037442e-08 # W/m², per CODATA 2018
f0 = 3630.78 # Jy for V == 0. Derived from S(Jy) = 1e23 * 10**-(AB+48.6)/2.5 = 10**(3.56-AB/2.5) = 10**((8.9-AB)/2.5). Not used directly, but is listed for comparison purposes. Sources: Oke&Gunn 1983, Fukugita et al 1996. This may implicitly be making Vega magnitude 0 in V, rather than eg: 0.03.

wattage = 100 #emitter net output in watts
area = wattage/(σ*Teff**4) # area in m² of our 100 W emitter
distance_close = ["close", np.linspace(1000e3,180000e3,717)] # 1000 to 180,000 km in 250 km steps, corresponding to the distances as een from the ground for an astrostationary orbit
distance_far = ["far", np.linspace(1.3e9,1.7e9,41)] # 1.3 million to 1.7 million km in 10,000 km steps, corresponding roughly to ground distances for an ESL-2 halo orbit, depending on details

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

#center wavelength and wavelngth range for WFIRST/NGRST instruments. Taken from https://wfirst.ipac.caltech.edu/sims/Param_db.html Phase C (July 29 2020)
#Instrument	band	min-max	mid (μm)
#WFI R062	0.480-0.760	0.620
#WFI Z087	0.760-0.977	0.869
#WFI Y106	0.927-1.192	1.060
#WFI J129	1.131-1.454	1.293
#WFI H158	1.380-1.774	1.577
#WFI F184	1.683-2.000	1.842
#WFI W146	0.927-2.000	1.464

WFI = [["R062",0.620e-6,0.140e-6],["Z087",0.8685e-6,0.1085e-6],["Y106",1.0595e-6,0.1325e-6],["J129",1.2925e-6,0.1615e-6],["H158",1.577e-6,0.197e-6],["F184",1.8415e-6,0.1585e-6],["W146",1.4635e-6,0.5365e-6]]

# LSST data derived from https://github.com/lsst/throughputs (as of 2020-11-04)
#LSST = 

fivethirtytwo = [["V",5448e-10,840e-10],["R",6407e-10,1580e-10],["g'",4639e-10,1280e-10],["r'",6122e-10,1150e-10],["R062",0.620e-6,0.140e-6],["W146",1.4635e-6,0.5365e-6]] #filters that overlap with a 532 nm laser
tensixtyfour = [["W146",1.4635e-6,0.5365e-6]] # None of the ground based bands I have currently overlap with 1064 nm laser aside perhaps from the wings of J. Maybe some sort of Y?
fourfortyfive = [["B",4361e-10,890e-10],["V",5448e-10,840e-10],["g'",4639e-10,1280e-10]]
sixthirtyeight = [["V",5448e-10,840e-10],["R",6407e-10,1580e-10],["r'",6122e-10,1150e-10],["i'",7439e-10,1230e-10],["R062",0.620e-6,0.140e-6]]
sixeighty = [["R",6407e-10,1580e-10],["I",7980e-10,1540e-10],["r'",6122e-10,1150e-10],["i'",7439e-10,1230e-10],["R062",0.620e-6,0.140e-6]]
eightoheight = [["I",7980e-10,1540e-10],["i'",7439e-10,1230e-10],["z'",8896e-10,1070e-10],["F184",1.8415e-6,0.1585e-6]]
dispersion = [4*np.pi,np.pi,1e-6] #angle of dispersion in steradians

for distance in [distance_close, distance_far]:
	for angle in dispersion:
		fig, ax = plt.subplots(2, 2, sharex=True, sharey=False)
		fig.suptitle("Johnson UBVRI Fluxes for a "+str(wattage)+" W, "+str(Teff)+" K blackbody ("+str(angle)[0:5]+" steradians)")
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
			flux = [luminosity/(angle*d**2) for d in distance[1]]
			fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
			Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
			ABmag = [8.9-2.5*np.log10(J) for J in Jy]
			ax[0,0].plot(distance[1]*1e-6, flux, label=name)
			ax[0,0].set_title("Flux (W/m²)")
			ax[0,1].plot(distance[1]*1e-6, fλ, label=name)
			ax[0,1].set_title("Fλ (W/m²/Å)")
			ax[1,0].plot(distance[1]*1e-6, Jy, label=name)
			ax[1,0].set_title("Fν (Jy)")
			ax[1,1].plot(distance[1]*1e-6, ABmag, label=name)
			ax[1,1].set_title("ABmag")
		ax[0,0].legend()
		ax[0,0].set_xlim(min(distance[1])*1e-6,max(distance[1])*1e-6)
		ax[1,0].set_xlabel("Distance (Mm)")
		ax[1,1].set_xlabel("Distance (Mm)")
		ax[0,0].set_yscale("log")
		ax[0,1].set_yscale("log")
		ax[1,0].set_yscale("log")
		fig.set_size_inches(8, 6.8)
		fig.tight_layout()
		#fig.savefig("Filters_Johnson_"+distance[0]+"_"+str(angle)[0:5]+".png", dpi=100)
		
		fig, bx = plt.subplots(2, 2, sharex=True, sharey=False)
		fig.suptitle("SDSS ugriz Fluxes for a "+str(wattage)+" W, "+str(Teff)+" K blackbody ("+str(angle)[0:5]+" steradians)")
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
			flux = [luminosity/(angle*d**2) for d in distance[1]]
			fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
			Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
			ABmag = [8.9-2.5*np.log10(J) for J in Jy]
			bx[0,0].plot(distance[1]*1e-6, flux, label=name)
			bx[0,0].set_title("Flux (W/m²)")
			bx[0,1].plot(distance[1]*1e-6, fλ, label=name)
			bx[0,1].set_title("Fλ (W/m²/Å)")
			bx[1,0].plot(distance[1]*1e-6, Jy, label=name)
			bx[1,0].set_title("Fν (Jy)")
			bx[1,1].plot(distance[1]*1e-6, ABmag, label=name)
			bx[1,1].set_title("ABmag")
		bx[0,0].legend()
		bx[0,0].set_xlim(min(distance[1])*1e-6,max(distance[1])*1e-6)
		bx[1,0].set_xlabel("Distance (Mm)")
		bx[1,1].set_xlabel("Distance (Mm)")
		bx[0,0].set_yscale("log")
		bx[0,1].set_yscale("log")
		bx[1,0].set_yscale("log")
		fig.set_size_inches(8, 6.8)
		fig.tight_layout()
		#fig.savefig("Filters_Sloan_"+distance[0]+"_"+str(angle)[0:5]+".png", dpi=100)
		
		fig, cx = plt.subplots(2, 2, sharex=True, sharey=False)
		fig.suptitle("2MASS JHK Fluxes for a "+str(wattage)+" W, "+str(Teff)+" K blackbody ("+str(angle)[0:5]+" steradians)")
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
			flux = [luminosity/(angle*d**2) for d in distance[1]]
			fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
			Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
			ABmag = [8.9-2.5*np.log10(J) for J in Jy]
			cx[0,0].plot(distance[1]*1e-6, flux, label=name)
			cx[0,0].set_title("Flux (W/m²)")
			cx[0,1].plot(distance[1]*1e-6, fλ, label=name)
			cx[0,1].set_title("Fλ (W/m²/Å)")
			cx[1,0].plot(distance[1]*1e-6, Jy, label=name)
			cx[1,0].set_title("Fν (Jy)")
			cx[1,1].plot(distance[1]*1e-6, ABmag, label=name)
			cx[1,1].set_title("ABmag")
		cx[0,0].legend()
		cx[0,0].set_xlim(min(distance[1])*1e-6,max(distance[1])*1e-6)
		cx[1,0].set_xlabel("Distance (Mm)")
		cx[1,1].set_xlabel("Distance (Mm)")
		cx[0,0].set_yscale("log")
		cx[0,1].set_yscale("log")
		cx[1,0].set_yscale("log")
		fig.set_size_inches(8, 6.8)
		fig.tight_layout()
		#fig.savefig("Filters_2MASS_"+distance[0]+"_"+str(angle)[0:5]+".png", dpi=100)
		
		
		fig, dx = plt.subplots(2, 2, sharex=True, sharey=False)
		fig.suptitle("NGRST WFI Fluxes for a "+str(wattage)+" W, "+str(Teff)+" K blackbody ("+str(angle)[0:5]+" steradians)")
		dx[0,0].set_prop_cycle(color=['red', 'purple', 'black', 'blue', 'green', 'orange', 'gray'])
		dx[0,1].set_prop_cycle(color=['red', 'purple', 'black', 'blue', 'green', 'orange', 'gray'])
		dx[1,0].set_prop_cycle(color=['red', 'purple', 'black', 'blue', 'green', 'orange', 'gray'])
		dx[1,1].set_prop_cycle(color=['red', 'purple', 'black', 'blue', 'green', 'orange', 'gray'])
		for band in WFI:
			name,λeff,Δλ = band
			Δν = bandwidth(λeff,Δλ)
			xs = np.linspace(λeff-Δλ,λeff+Δλ,1000)
			ys = [planck(x,Teff) for x in xs]
			power = np.trapz(ys,xs)
			luminosity = power*area
			flux = [luminosity/(angle*d**2) for d in distance[1]]
			fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
			Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
			ABmag = [8.9-2.5*np.log10(J) for J in Jy]
			dx[0,0].plot(distance[1]*1e-6, flux, label=name)
			dx[0,0].set_title("Flux (W/m²)")
			dx[0,1].plot(distance[1]*1e-6, fλ, label=name)
			dx[0,1].set_title("Fλ (W/m²/Å)")
			dx[1,0].plot(distance[1]*1e-6, Jy, label=name)
			dx[1,0].set_title("Fν (Jy)")
			dx[1,1].plot(distance[1]*1e-6, ABmag, label=name)
			dx[1,1].set_title("ABmag")
		dx[0,0].legend()
		dx[0,0].set_xlim(min(distance[1])*1e-6,max(distance[1])*1e-6)
		dx[1,0].set_xlabel("Distance (Mm)")
		dx[1,1].set_xlabel("Distance (Mm)")
		dx[0,0].set_yscale("log")
		dx[0,1].set_yscale("log")
		dx[1,0].set_yscale("log")
		fig.set_size_inches(8, 6.8)
		fig.tight_layout()
		#fig.savefig("Filters_WFI_"+distance[0]+"_"+str(angle)[0:5]+".png", dpi=100)


wattage = 1
efficiency = 100 #%
for distance in [distance_close, distance_far]:
	for angle in dispersion:
		fig, ex = plt.subplots(2, 2, sharex=True, sharey=False)
		fig.suptitle("Fluxes for a "+str(wattage)+" W, 532 nm LED ("+str(angle)[0:5]+" steradians)")
		for band in fivethirtytwo:
			name,λeff,Δλ = band
			Δν = bandwidth(λeff,Δλ)
			luminosity = wattage*efficiency/100
			flux = [luminosity/(angle*d**2) for d in distance[1]]
			fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
			Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
			ABmag = [8.9-2.5*np.log10(J) for J in Jy]
			ex[0,0].plot(distance[1]*1e-6, flux, label=name)
			ex[0,0].set_title("Flux (W/m²)")
			ex[0,1].plot(distance[1]*1e-6, fλ, label=name)
			ex[0,1].set_title("Fλ (W/m²/Å)")
			ex[1,0].plot(distance[1]*1e-6, Jy, label=name)
			ex[1,0].set_title("Fν (Jy)")
			ex[1,1].plot(distance[1]*1e-6, ABmag, label=name)
			ex[1,1].set_title("ABmag")
		ex[0,0].legend()
		ex[0,0].set_xlim(min(distance[1])*1e-6,max(distance[1])*1e-6)
		ex[1,0].set_xlabel("Distance (Mm)")
		ex[1,1].set_xlabel("Distance (Mm)")
		ex[0,0].set_yscale("log")
		ex[0,1].set_yscale("log")
		ex[1,0].set_yscale("log")
		fig.set_size_inches(8, 6.8)
		fig.tight_layout()
		#fig.savefig("Filters_532nm_"+distance[0]+"_"+str(angle)[0:5]+".png", dpi=100)
		
		fig, fx = plt.subplots(2, 2, sharex=True, sharey=False)
		fig.suptitle("Fluxes for a "+str(wattage)+" W, 1064 nm LED ("+str(angle)[0:5]+" steradians)")
		for band in tensixtyfour:
			name,λeff,Δλ = band
			Δν = bandwidth(λeff,Δλ)
			luminosity = wattage*efficiency/100
			flux = [luminosity/(angle*d**2) for d in distance[1]]
			fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
			Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
			ABmag = [8.9-2.5*np.log10(J) for J in Jy]
			fx[0,0].plot(distance[1]*1e-6, flux, label=name)
			fx[0,0].set_title("Flux (W/m²)")
			fx[0,1].plot(distance[1]*1e-6, fλ, label=name)
			fx[0,1].set_title("Fλ (W/m²/Å)")
			fx[1,0].plot(distance[1]*1e-6, Jy, label=name)
			fx[1,0].set_title("Fν (Jy)")
			fx[1,1].plot(distance[1]*1e-6, ABmag, label=name)
			fx[1,1].set_title("ABmag")
		fx[0,0].legend()
		fx[0,0].set_xlim(min(distance[1])*1e-6,max(distance[1])*1e-6)
		fx[1,0].set_xlabel("Distance (Mm)")
		fx[1,1].set_xlabel("Distance (Mm)")
		fx[0,0].set_yscale("log")
		fx[0,1].set_yscale("log")
		fx[1,0].set_yscale("log")
		fig.set_size_inches(8, 6.8)
		fig.tight_layout()
		#fig.savefig("Filters_1064nm_"+distance[0]+"_"+str(angle)[0:5]+".png", dpi=100)
		
		fig, gx = plt.subplots(2, 2, sharex=True, sharey=False)
		fig.suptitle("Fluxes for a "+str(wattage)+" W, 445 nm LED ("+str(angle)[0:5]+" steradians)")
		for band in fourfortyfive:
			name,λeff,Δλ = band
			Δν = bandwidth(λeff,Δλ)
			luminosity = wattage*efficiency/100
			flux = [luminosity/(angle*d**2) for d in distance[1]]
			fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
			Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
			ABmag = [8.9-2.5*np.log10(J) for J in Jy]
			gx[0,0].plot(distance[1]*1e-6, flux, label=name)
			gx[0,0].set_title("Flux (W/m²)")
			gx[0,1].plot(distance[1]*1e-6, fλ, label=name)
			gx[0,1].set_title("Fλ (W/m²/Å)")
			gx[1,0].plot(distance[1]*1e-6, Jy, label=name)
			gx[1,0].set_title("Fν (Jy)")
			gx[1,1].plot(distance[1]*1e-6, ABmag, label=name)
			gx[1,1].set_title("ABmag")
		gx[0,0].legend()
		gx[0,0].set_xlim(min(distance[1])*1e-6,max(distance[1])*1e-6)
		gx[1,0].set_xlabel("Distance (Mm)")
		gx[1,1].set_xlabel("Distance (Mm)")
		gx[0,0].set_yscale("log")
		gx[0,1].set_yscale("log")
		gx[1,0].set_yscale("log")
		fig.set_size_inches(8, 6.8)
		fig.tight_layout()
		#fig.savefig("Filters_445nm_"+distance[0]+"_"+str(angle)[0:5]+".png", dpi=100)
		
		fig, hx = plt.subplots(2, 2, sharex=True, sharey=False)
		fig.suptitle("Fluxes for a "+str(wattage)+" W, 638 nm LED ("+str(angle)[0:5]+" steradians)")
		for band in sixthirtyeight:
			name,λeff,Δλ = band
			Δν = bandwidth(λeff,Δλ)
			luminosity = wattage*efficiency/100
			flux = [luminosity/(angle*d**2) for d in distance[1]]
			fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
			Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
			ABmag = [8.9-2.5*np.log10(J) for J in Jy]
			hx[0,0].plot(distance[1]*1e-6, flux, label=name)
			hx[0,0].set_title("Flux (W/m²)")
			hx[0,1].plot(distance[1]*1e-6, fλ, label=name)
			hx[0,1].set_title("Fλ (W/m²/Å)")
			hx[1,0].plot(distance[1]*1e-6, Jy, label=name)
			hx[1,0].set_title("Fν (Jy)")
			hx[1,1].plot(distance[1]*1e-6, ABmag, label=name)
			hx[1,1].set_title("ABmag")
		hx[0,0].legend()
		hx[0,0].set_xlim(min(distance[1])*1e-6,max(distance[1])*1e-6)
		hx[1,0].set_xlabel("Distance (Mm)")
		hx[1,1].set_xlabel("Distance (Mm)")
		hx[0,0].set_yscale("log")
		hx[0,1].set_yscale("log")
		hx[1,0].set_yscale("log")
		fig.set_size_inches(8, 6.8)
		fig.tight_layout()
		#fig.savefig("Filters_638nm_"+distance[0]+"_"+str(angle)[0:5]+".png", dpi=100)
		
		fig, ix = plt.subplots(2, 2, sharex=True, sharey=False)
		fig.suptitle("Fluxes for a "+str(wattage)+" W, 680 nm LED ("+str(angle)[0:5]+" steradians)")
		for band in sixeighty:
			name,λeff,Δλ = band
			Δν = bandwidth(λeff,Δλ)
			luminosity = wattage*efficiency/100
			flux = [luminosity/(angle*d**2) for d in distance[1]]
			fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
			Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
			ABmag = [8.9-2.5*np.log10(J) for J in Jy]
			ix[0,0].plot(distance[1]*1e-6, flux, label=name)
			ix[0,0].set_title("Flux (W/m²)")
			ix[0,1].plot(distance[1]*1e-6, fλ, label=name)
			ix[0,1].set_title("Fλ (W/m²/Å)")
			ix[1,0].plot(distance[1]*1e-6, Jy, label=name)
			ix[1,0].set_title("Fν (Jy)")
			ix[1,1].plot(distance[1]*1e-6, ABmag, label=name)
			ix[1,1].set_title("ABmag")
		ix[0,0].legend()
		ix[0,0].set_xlim(min(distance[1])*1e-6,max(distance[1])*1e-6)
		ix[1,0].set_xlabel("Distance (Mm)")
		ix[1,1].set_xlabel("Distance (Mm)")
		ix[0,0].set_yscale("log")
		ix[0,1].set_yscale("log")
		ix[1,0].set_yscale("log")
		fig.set_size_inches(8, 6.8)
		fig.tight_layout()
		#fig.savefig("Filters_680nm_"+distance[0]+"_"+str(angle)[0:5]+".png", dpi=100)
		
		fig, fx = plt.subplots(2, 2, sharex=True, sharey=False)
		fig.suptitle("Fluxes for a "+str(wattage)+" W, 808 nm LED ("+str(angle)[0:5]+" steradians)")
		for band in eightoheight:
			name,λeff,Δλ = band
			Δν = bandwidth(λeff,Δλ)
			luminosity = wattage*efficiency/100
			flux = [luminosity/(angle*d**2) for d in distance[1]]
			fλ = [f/Δλ *1e10 for f in flux] # so units are W/m²/Å
			Jy = [f/Δν * 1e26 for f in flux] # 1 Jy = 1e-26 W/m²/Hz
			ABmag = [8.9-2.5*np.log10(J) for J in Jy]
			fx[0,0].plot(distance[1]*1e-6, flux, label=name)
			fx[0,0].set_title("Flux (W/m²)")
			fx[0,1].plot(distance[1]*1e-6, fλ, label=name)
			fx[0,1].set_title("Fλ (W/m²/Å)")
			fx[1,0].plot(distance[1]*1e-6, Jy, label=name)
			fx[1,0].set_title("Fν (Jy)")
			fx[1,1].plot(distance[1]*1e-6, ABmag, label=name)
			fx[1,1].set_title("ABmag")
		fx[0,0].legend()
		fx[0,0].set_xlim(min(distance[1])*1e-6,max(distance[1])*1e-6)
		fx[1,0].set_xlabel("Distance (Mm)")
		fx[1,1].set_xlabel("Distance (Mm)")
		fx[0,0].set_yscale("log")
		fx[0,1].set_yscale("log")
		fx[1,0].set_yscale("log")
		fig.set_size_inches(8, 6.8)
		fig.tight_layout()
		#fig.savefig("Filters_808nm_"+distance[0]+"_"+str(angle)[0:5]+".png", dpi=100)
		
plt.show()

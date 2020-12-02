import numpy as np
from matplotlib import pyplot as plt
import csv

# covers UBRVI, ugriz, JHK, NGRST WFI bands in relevant cases for a 2700/3000/6500 K source at isotropic and limited diffraction angle for first 4 sets of graphs. Then for 532 nm / 1064 nm lasers for the second 2 sets. (for those, 100% of emitted light is assumed to be completely captured by the filters in question)
# photons/m² instead of f_λ


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

# LSST/NGRST/JWST data taken from wavelength table 2, but using mean and width. Min and max values *will not* match up!
# Also, some values in the table are of questionable accuracy. Best estimate of what's going on used.
LSST = [["u",0.3694e-6,0.0473e-6],["g",0.4841e-6,0.1253e-6/2],["r",0.6358e-6,0.1207e-6/2],["i",0.7560e-6,0.1175e-6/2],["z",0.8701e-6,0.0998e-6],["y",0.9749e-6,0.0872e-6]]
#Mostly WFI, but last two are GRS and PRS
NGRST = [["F062",0.620e-6,0.280e-6/2],["F087",0.869e-6,0.217e-6/2],["F106",1.060e-6,0.265e-6/2],["F129",1.293e-6,0.323e-6/2],["F158",1.577e-6,0.394e-6/2],["F184",1.842e-6,0.317e-6/2],["W146",1.464e-6,1.030e-6/2],["G150",1.465e-6,0.930e-6/2],["P120",0.975e-6,1.05e-6/2]]
# So much NIRCam
JWST = [["F070W",0.7088e-6,0.1213e-6],["F090W",0.9083e-6,0.1773e-6],["F115W",1.1624e-6,0.2055e-6],["F140M",1.4074e-6,0.1367e-6],["F150W",1.5104e-6,0.2890e-6],["F162M",1.6297e-6,0.1626e-6],["F182M",1.8494e-6,0.2251e-6],["F200W",2.0028e-6,0.4190e-6],["F210M",2.0982e-6,0.2055e-6],["F250M",2.5049e-6,0.1783e-6],["F277W",2.7845e-6,0.6615e-6],["F300M",2.9940e-6,0.3256e-6],["F335M",3.3675e-6,0.3389e-6],["F356W",3.5935e-6,0.7239e-6],["F360M",3.6298e-6,0.3585e-6],["F410M",4.0887e-6,0.4263e-6],["F430M",4.2829e-6,0.2295e-6],["F444W",4.4394e-6,1.0676e-6],["F460M",4.6316e-6,0.2309e-6],["F480M",4.8213e-6,0.3141e-6]]

fivethirtytwo = [["V",5448e-10,840e-10],["R",6407e-10,1580e-10],["g'",4639e-10,1280e-10],["r'",6122e-10,1150e-10],["R062",0.620e-6,0.140e-6],["W146",1.4635e-6,0.5365e-6]] #filters that overlap with a 532 nm laser
tensixtyfour = [["Y106",1.05950e-6,0.1325e-6],["W146",1.4635e-6,0.5365e-6]] # None of the ground based bands I have currently overlap with 1064 nm laser aside perhaps from the wings of J. Maybe some sort of Y?
fourfortyfive = [["B",4361e-10,890e-10],["V",5448e-10,840e-10],["g'",4639e-10,1280e-10]]
sixthirtyeight = [["V",5448e-10,840e-10],["R",6407e-10,1580e-10],["r'",6122e-10,1150e-10],["i'",7439e-10,1230e-10],["R062",0.620e-6,0.140e-6]]
sixeighty = [["R",6407e-10,1580e-10],["I",7980e-10,1540e-10],["r'",6122e-10,1150e-10],["i'",7439e-10,1230e-10],["R062",0.620e-6,0.140e-6]]
eightoheight = [["I",7980e-10,1540e-10],["i'",7439e-10,1230e-10],["z'",8896e-10,1070e-10]]
dispersion = [4*np.pi,np.pi,1e-6] #angle of dispersion in steradians

list0 = Johnson + SDSS + twoMASS + WFI #hand assembled
list1 = LSST+NGRST+JWST #from table 2
wavelengths = np.reshape(np.linspace(100e-9,5500e-9,5401),(5401,1)) # 100 - 5500 nm in 10 nm chunks
overlap = np.reshape(np.zeros(5401),(5401,1))

#Johnson, SDSS, 2MASS, WFI, LSST, NGRST, JWST
list = np.concatenate((wavelengths,np.zeros((5401,7))),axis=1)
#list2 = np.concatenate((Johnson,SDSS,twoMASS,WFI,LSST,NGRST,JWST),axis=1)
for k in np.arange(0,list.shape[0]): #row
	λ = list[k][0]
	for name,λeff,Δλ in Johnson:
		if (λ >= (λeff-Δλ)) and (λ <= (λeff+Δλ)):
			list[k][1] += 1
	for name,λeff,Δλ in SDSS:
		if (λ >= (λeff-Δλ)) and (λ <= (λeff+Δλ)):
			list[k][2] += 1
	for name,λeff,Δλ in twoMASS:
		if (λ >= (λeff-Δλ)) and (λ <= (λeff+Δλ)):
			list[k][3] += 1
	for name,λeff,Δλ in WFI:
		if (λ >= (λeff-Δλ)) and (λ <= (λeff+Δλ)):
			list[k][4] += 1
	for name,λeff,Δλ in LSST:
		if (λ >= (λeff-Δλ)) and (λ <= (λeff+Δλ)):
			list[k][5] += 1
	for name,λeff,Δλ in NGRST:
		if (λ >= (λeff-Δλ)) and (λ <= (λeff+Δλ)):
			list[k][6] += 1
	for name,λeff,Δλ in JWST:
		if (λ >= (λeff-Δλ)) and (λ <= (λeff+Δλ)):
			list[k][7] += 1
with open('overlaps.csv', 'w') as csvfile: #all filter sets
	spamwriter = csv.writer(csvfile)
	spamwriter.writerow(["wavelength,Johnson,SDSS,2MASS,WFI,LSST,NGRST,JWST"])
	for line in list:
		spamwriter.writerow(line)

list = np.concatenate((wavelengths,overlap),axis=1)
for k in np.arange(0,list.shape[0]):
	λ = list[k][0]
	for name,λeff,Δλ in list1:
		if (λ >= (λeff-Δλ)) and (λ <= (λeff+Δλ)):
			list[k][1] += 1
with open('list1.csv', 'w') as csvfile: # subset
	spamwriter = csv.writer(csvfile)
	spamwriter.writerow(["wavelength,number"])
	for line in list:
		spamwriter.writerow(line)

'''
for Teff in Teffs:
	area = wattage/(σ*Teff**4) # area in m² of our 100 W emitter
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
				photons = [planck(x,Teff)*x/h/c for x in xs]
				photon_band = np.trapz(photons,xs)
				photon_flux = [photon_band*area/(angle*d**2) for d in distance[1]]
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
			#fig.savefig("Filters_Johnson_"+str(Teff)+"K_"+distance[0]+"_"+str(angle)[0:5]+".png", dpi=100)
			
'''
#plt.show()

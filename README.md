# solid-lamp
Informal repository of ORCAS related calculations

Currently contains fluxes.py, vega.py, and their default outputs (fluxes.csv and vega.csv).

Both of these are approximate flux calculations in a number of bands (Johnson UBVRI, SDSSS ugriz, 2MASS JHK), giving results in terms of fλ (W/m²/m), fν (Janskys), and ABmag. fluxes.py is currently set to a 100 W, 3000 K blackbody 250 km away, while vega.py uses Vega (as values taken from Wikipedia). Informally, the vega check suggests that accuracy is no better than 0.05 mag, and may be 0.1 mag or worse.

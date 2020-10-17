# solid-lamp
Informal repository of ORCAS related calculations

Currently contains fluxes.py, vega.py, their default outputs (fluxes.csv and vega.csv), and mag_graphs.py.

Both of fluxes.py and vega.py are approximate flux calculations in a number of bands (Johnson UBVRI, SDSSS ugriz, 2MASS JHK), giving results in terms of fλ (W/m²/m), fν (Janskys), and ABmag. fluxes.py is currently set to a 100 W, 3000 K blackbody 250 km away, while vega.py uses Vega (as values taken from Wikipedia). Atmospheric effects are not considered. Informally, the vega check suggests that accuracy is no better than 0.05 mag, and may be 0.1 mag or worse.

mag_graphs.py generates graphs of fluxes in the bands with the same assumptions, except over distances from 100 to 1200 km.

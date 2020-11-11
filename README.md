# solid-lamp
Informal repository of ORCAS related calculations

Currently contains a photon count estimate (photon_counts.ods), fluxes.py, vega.py, their default outputs (fluxes.csv and vega.csv), and mag_graphs.py.

Both of fluxes.py and vega.py are approximate flux calculations in a number of bands (Johnson UBVRI, SDSSS ugriz, 2MASS JHK), giving results in terms of fλ (W/m²/m), fν (Janskys), and ABmag. fluxes.py is currently set to a 100 W, 3000 K blackbody 250 km away, while vega.py uses Vega (as values taken from Wikipedia). Atmospheric effects are not considered. Informally, the vega check suggests that accuracy is no better than 0.05 mag, and may be 0.1 mag or worse.

mag_graphs.py generates graphs of fluxes in the bands with the same assumptions, except over distances from 100 to 1200 km.


These files are being altered with new versions increasing the appended numbers. Currently (2020-11-11), the most recent is mag_graphs4.py, which considers 2 kinds of sources:
* 100 W, 3000 K blackbody
* Various monochromatic LEDs. For these, the bandwidth is defined by the detector only. (So a broader band will result in lower flux, even though it recieves the same number of watts)

It also looks at two distance ranges: 
* "close" at 1 - 180 Mm (comparable to a highly elliptical orbit as seen from the ground or halo orbit as seen by another spacecraft at the same Lagrange point)
* "far" at 1.3 - 1.7 Gm (L2 halo orbit as seen from the ground)

Current detector bands are: Johnson UBVRI, SDSS ugriz, 2MASS JHK, and WFIRST WFI. These are representative of most systems, but a number of actual instruments that would find this useful for calibration (eg: VRO ugrizy)

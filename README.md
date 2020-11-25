# solid-lamp
Informal repository of ORCAS related calculations

Currently contains a photon count estimate (photon_counts.ods), fluxes.py, vega.py, their default outputs (fluxes.csv and vega.csv), and mag_graphs.py.

Both of fluxes.py and vega.py are approximate flux calculations in a number of bands (Johnson UBVRI, SDSSS ugriz, 2MASS JHK), giving results in terms of fλ (W/m²/m), fν (Janskys), and ABmag. fluxes.py is currently set to a 100 W, 3000 K blackbody 250 km away, while vega.py uses Vega (as values taken from Wikipedia). Atmospheric effects are not considered. Informally, the vega check suggests that accuracy is no better than 0.05 mag, and may be 0.1 mag or worse.

mag_graphs.py generates graphs of fluxes in the bands with the same assumptions, except over distances from 100 to 1200 km.

divergence.py considers beam divergence for visibility "all over" the earth. That is, for a given zenith angle (default is 60 degrees / airmass of 2), what sort of beam divergence (as a half-angle) is needed for an evenly distributed beam to be visible. Atmospheric refraction is not considered.


These files are being altered with new versions increasing the appended numbers. Currently (2020-11-11), the most recent is mag_graphs4.py, which considers 2 kinds of sources:
* 100 W, 3000 K blackbody
* Various monochromatic LEDs. For these, the bandwidth is defined by the detector only. (So a broader band will result in lower flux, even though it recieves the same number of watts)

It also looks at two distance ranges: 
* "close" at 1 - 180 Mm (comparable to a highly elliptical orbit as seen from the ground or halo orbit as seen by another spacecraft at the same Lagrange point)
* "far" at 1.3 - 1.7 Gm (L2 halo orbit as seen from the ground)

Current detector bands are: Johnson UBVRI, SDSS ugriz, 2MASS JHK, and WFIRST WFI. These are representative of most systems, but a number of actual instruments that would find this useful for calibration (eg: VRO ugrizy)


As of 2020-11-24, mag_graphs5.py expands this to 3 blackbodies (2700 K, 3000 K, 6500 K) to compare with typical LED sources that simulate them, as well as replacing f_λ with photons/m²/s. This risks giving innacurate impressions, since you're not getting a true 6500 K blackbody without looking at a star!

# Comparison with physical sources
Incandecent bulbs are more or less comparable to blackbodies over a large range. LED sources (eg: Gamma Scientific SpectralLED RS-7-1 can fake it over a somewhat narrower range. In particular, it focuses on 380 to 1000 nm, which may present issues with eg: Johnson U band, and 2MASS JHK. This source is also insufficiently rugged for spacecraft operations, and has a bunch of precision angular limitations that would have to be considered.

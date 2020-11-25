'''
Estimating maximum beam divergence as a function of distance from Earth for maximum viewing. Narrower is of course possible, if somewhat useage limiting.
'''
import numpy as np
from matplotlib import pyplot as plt

r_earth = 6378100.0 # Earth radius in meters, per IAU 2015 Resoultion B 3
horizon_angle = np.pi/6 # 30 degrees in radians
distance_close = np.linspace(1000e3,180000e3,717) # 1000 to 180,000 km in 250 km steps, corresponding to the distances as een from the ground for an astrostationary orbit
angle_close = [np.arctan(np.cos(horizon_angle)*r_earth/x) for x in distance_close]
distance_far = np.linspace(1.3e9,1.7e9,41) # 1.3 million to 1.7 million km in 10,000 km steps, corresponding roughly to ground distances for an ESL-2 halo orbit, depending on details
angle_far = [np.arctan(np.cos(horizon_angle)*r_earth/x) for x in distance_far]

fig, ax = plt.subplots(1,2)
ax[0].plot(distance_close, angle_close)
ax[1].plot(distance_far, angle_far)

ax[0].set_xlabel("distance (m)")
ax[1].set_xlabel("distance (m)")
ax[0].set_ylabel("Beam divergence half-angle (radians)")
fig.suptitle("Requirement for visibility at 2 airmass / 60 degree zenith angle")
ax[0].set_xlim(min(distance_close),max(distance_close))
ax[0].set_ylim(0,np.pi/2)
ax[1].set_xlim(min(distance_far),max(distance_far))
ax[1].set_ylim(0.0032,0.0043)
ax[0].set_xscale("log")
fig.set_size_inches(8, 5)
fig.tight_layout()
#fig.savefig("Beam_divergence_distance.png", dpi=100)
plt.show()

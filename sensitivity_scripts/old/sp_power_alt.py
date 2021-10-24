import os
import hickle
import numpy as np
import matplotlib.pyplot as plt


def sp(p, h):
    os.system("/home/daniel/anaconda3/envs/thesis_code/bin/python3.8 ../matlab_engine_sizing_new.py " + str(p) + " "
             + str(h) + " 0.3 0.2 0")
    res = hickle.load("sensitivity_output.hkl")
    return p / sum(res[:4]) / 1e3


# Power - altitude
ps = np.linspace(5e5, 1e7, 10)
hs = np.linspace(2e3, 1e4, 10)
ps, hs = np.meshgrid(ps, hs)
spv = np.vectorize(sp)
sps = spv(ps, hs)

plt.imshow(sps, extent=[ps.min()/1e6, ps.max()/1e6, hs.min()/1e3, hs.max()/1e3], vmin=sps.min(), vmax=sps.max(),
           aspect="auto", interpolation="bicubic", origin="lower")
cbar = plt.colorbar()
cbar.set_label("Specific power in kW/kg")
plt.xlabel("Power in MW")
plt.ylabel("Altitude in km")
plt.savefig("sp_power_alt.pdf")

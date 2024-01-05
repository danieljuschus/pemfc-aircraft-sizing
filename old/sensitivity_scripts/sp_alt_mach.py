# Effect of altitude and Mach number on specific power - Fig. 3.14
import os
import hickle
import numpy as np
import matplotlib.pyplot as plt


def sp(m, h):
    p = 1e6
    os.system("/home/daniel/anaconda3/envs/thesis_code/bin/python3.8 ../matlab_engine_sizing_new.py " + str(p) + " "
              + str(h) + " " + str(m) + " 0.2 0")
    res = hickle.load("sensitivity_output.hkl")
    return p / sum(res[:4]) / 1e3


# Power - altitude
hs = np.linspace(2e3, 1e4, 10)
ms = np.linspace(0.05, 0.5, 10)
ms, hs = np.meshgrid(ms, hs)
spv = np.vectorize(sp)
sps = spv(ms, hs)

plt.imshow(sps, extent=[ms.min(), ms.max(), hs.min()/1e3, hs.max()/1e3], vmin=sps.min(), vmax=sps.max(),
           aspect="auto", interpolation="bicubic", origin="lower")
cbar = plt.colorbar()
cbar.set_label("Specific power in kW/kg")
plt.xlabel("Mach number")
plt.ylabel("Altitude in km")
plt.savefig("sp_alt_mach.pdf")

import os
import hickle
import numpy as np
import matplotlib.pyplot as plt


def sp(p, m):
    os.system("/home/daniel/anaconda3/envs/thesis_code/bin/python3.8 ../matlab_engine_sizing_new.py "
              + str(p) + " 5000 " + str(m) + " 0.2 0")
    res = hickle.load("sensitivity_output.hkl")
    return p / sum(res[:4]) / 1e3


# Power - altitude
ps = np.linspace(1e5, 1e6, 10)
ms = np.linspace(0.05, 0.6, 10)
ps, ms = np.meshgrid(ps, ms)
spv = np.vectorize(sp)
sps = spv(ps, ms)

plt.imshow(sps, extent=np.array([ps.min()/1e3, ps.max()/1e3, ms.min(), ms.max()]), vmin=sps.min(), vmax=sps.max(),
           aspect="auto", interpolation="bicubic", origin="lower")
cbar = plt.colorbar()
cbar.set_label("Specific power in kW/kg")
plt.xlabel("Power in kW")
plt.ylabel("Mach number")
plt.savefig("sp_power_mach.pdf")

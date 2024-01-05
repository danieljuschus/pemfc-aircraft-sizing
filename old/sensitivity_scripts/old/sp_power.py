import os
import hickle
import matplotlib.pyplot as plt
import numpy as np

ps = np.linspace(1e4, 1e6, 10)
sps, sp_ss, sp_cs, sp_hs, sp_hxs = [], [], [], [], []
for p in ps:
    os.system("/home/daniel/anaconda3/envs/thesis_code/bin/python3.8 ../matlab_engine_sizing_new.py "
              + str(p) + " 5000 0.3 0.2 0")
    res = hickle.load("sensitivity_output.hkl")
    sp_s, sp_c, sp_h, sp_hx, _ = [p/mi/1e3 for mi in res]
    sp_ss.append(sp_s)
    sp_cs.append(sp_c)
    sp_hs.append(sp_h)
    sp_hxs.append(sp_hx)
    sps.append(p/sum(res[:4])/1e3)

plt.plot(ps, sps)
plt.show()


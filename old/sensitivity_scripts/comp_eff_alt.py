# Effect of altitude and stack and compressor mass at different compressor efficiencies - Fig 4.8
import numpy as np
import matplotlib.pyplot as plt
import os
import hickle

alts = np.linspace(1, 1.1e4, 10)

comp_effs = np.arange(0.6, 1, 0.1)

# ms, m_ss, m_cs, m_hs, m_hxs, etas = [], [], [], [], [], []
# for alt in alts:
#     os.system("/home/daniel/anaconda3/envs/thesis_code/bin/python3.8 ../matlab_engine_sizing_new.py 1e6 "
#               + str(alt) + " 0.3 0.2 0 1")
#     res = hickle.load("sensitivity_output.hkl")
#     m_s, m_c, m_h, m_hx, eta = res
#     m_ss.append(m_s)
#     m_cs.append(m_c)
#     m_hs.append(m_h)
#     m_hxs.append(m_hx)
#     etas.append(eta)
#     ms.append(m_s + m_c + m_h + m_hx)
#
# hickle.dump([ms, m_ss, m_cs, m_hs, m_hxs], "comp_eff_alt_100.hkl")
cmapin = np.linspace(0, 100, len(comp_effs))
for i, comp_eff in enumerate(comp_effs):
    comp_eff = str(int(round(comp_eff, 1)*100))
    res = hickle.load("comp_eff_alt_" + comp_eff + ".hkl")
    plt.plot(alts[1:]/1e3, res[1][1:], label=r"Stack ($\eta_{is,comp}$=" + comp_eff + "%)")
for i, comp_eff in enumerate(comp_effs):
    comp_eff = str(int(round(comp_eff, 1) * 100))
    res = hickle.load("comp_eff_alt_" + comp_eff + ".hkl")
    plt.plot(alts[1:]/1e3, res[2][1:], label=r"Compressor ($\eta_{is,comp}$=" + comp_eff + "%)")
plt.grid()
plt.legend()
plt.xlabel("Altitude in km")
plt.ylabel("Component mass in kg")
plt.savefig("figures/comp_eff_alt.pdf")
plt.show()

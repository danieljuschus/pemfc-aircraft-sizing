# Effect of altitude on FC system without compressor - Fig 4.7
import os
import hickle
import numpy as np
import matplotlib.pyplot as plt

alts = np.linspace(0, 10000, 10)


def masses(h):
    p = 1e6
    os.system("/home/daniel/anaconda3/envs/thesis_code/bin/python3.8 ../matlab_engine_sizing_new.py " + str(p) + " "
              + str(h) + " 0.3 0.2 0 0")
    res = hickle.load("sensitivity_output.hkl")
    m_s, m_c, m_h, m_hx, eta = res
    return m_s, m_c, m_h, m_hx, eta


ms, m_ss, m_cs, m_hs, m_hxs, etas = [], [], [], [], [], []
for alt in alts:
    m_s, m_c, m_h, m_hx, eta = masses(alt)
    m_ss.append(m_s)
    m_cs.append(m_c)
    m_hs.append(m_h)
    m_hxs.append(m_hx)
    etas.append(eta)
    ms.append(m_s + m_c + m_h + m_hx)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(alts, m_ss, label="Stack mass")
ax1.plot(alts, m_cs, label="Compressor mass")
ax1.plot(alts, m_hs, label="Humidifier mass")
ax1.plot(alts, m_hxs, label="Heat exchanger mass")
ax1.plot(alts, ms, label="Total FC system mass")
ax2.plot(alts, etas, "k--", label="FC system efficiency")
ax1.grid()
# leg = ax1.legend(title="Mass of")
# leg.remove()
ax1.set_xlabel(r"Altitude in m")
ax1.set_ylabel("Mass in $kg$")
ax2.set_ylabel("Cell efficiency")
# ax1.set_zorder(1)
# ax2.yaxis.label.set_color("g")
# ax2.spines["right"].set_edgecolor("g")
# ax2.tick_params(axis='y', colors="g")
# ax2.legend()
# ax2.add_artist(leg)
fig.legend()
plt.tight_layout()
plt.show()
fig.savefig("figures/alt_no_comp.pdf")

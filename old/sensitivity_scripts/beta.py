# Effect of compressore pressure ration on FC system - Fig 4.6
import os
import hickle
import numpy as np
import matplotlib.pyplot as plt

betas = np.linspace(1.5, 2.5, 10)
h = 5000


def sp(beta, h):
    p = 1e6
    os.system("/home/daniel/anaconda3/envs/thesis_code/bin/python3.8 ../matlab_engine_sizing_new.py " + str(p) + " "
              + str(h) + " 0.3 0.2 " + str(beta) + " 1")
    res = hickle.load("sensitivity_output.hkl")
    m_s, m_c, m_h, m_hx, eta = res
    return m_s, m_c, m_h, m_hx, eta


ms, m_ss, m_cs, m_hs, m_hxs, etas = [], [], [], [], [], []
for betai in betas:
    m_s, m_c, m_h, m_hx, eta = sp(betai, h)
    m_ss.append(m_s)
    m_cs.append(m_c)
    m_hs.append(m_h)
    m_hxs.append(m_hx)
    etas.append(eta)
    ms.append(m_s + m_c + m_h + m_hx)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(betas, m_ss, label="Stack mass")
ax1.plot(betas, m_cs, label="Compressor mass")
ax1.plot(betas, m_hs, label="Humidifier mass")
ax1.plot(betas, m_hxs, label="Heat exchanger mass")
ax1.plot(betas, ms, label="Total FC system mass")
ax2.plot(betas, etas, "k--", label="FC system efficiency")
ax1.grid()
# leg = ax1.legend(title="Mass of")
# leg.remove()
ax1.set_xlabel(r"Compressor pressure ratio $\beta$")
ax1.set_ylabel("Mass in $kg$")
ax2.set_ylabel("FC system efficiency")
# ax1.set_zorder(1)
# ax2.yaxis.label.set_color("g")
# ax2.spines["right"].set_edgecolor("g")
# ax2.tick_params(axis='y', colors="g")
# ax2.legend()
# ax2.add_artist(leg)
fig.legend()  # legend moved manually in Inkscape for report
plt.tight_layout()
plt.show()
fig.savefig("figures/beta.pdf")

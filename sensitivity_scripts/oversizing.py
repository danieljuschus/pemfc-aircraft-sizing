# Effect of oversizing factor on FC system - Fig 4.9
import os
import hickle
import numpy as np
import matplotlib.pyplot as plt

ovs = np.arange(0, 0.6, 0.1)

ms, m_ss, m_cs, m_hs, m_hxs, etas = [], [], [], [], [], []
for ov in ovs:
    os.system("/home/daniel/anaconda3/envs/thesis_code/bin/python3.8 ../matlab_engine_sizing_new.py 3e6 7000 0.3 "
              + str(ov) + " 0 1")
    m_s, m_c, m_h, m_hx, eta = hickle.load("sensitivity_output.hkl")
    m_ss.append(m_s)
    m_cs.append(m_c)
    m_hs.append(m_h)
    m_hxs.append(m_hx)
    etas.append(eta)
    ms.append(sum(hickle.load("sensitivity_output.hkl")))

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(ovs, m_ss, label="Stack mass")
ax1.plot(ovs, m_cs, label="Compressor mass")
ax1.plot(ovs, m_hs, label="Humidifier mass")
ax1.plot(ovs, m_hxs, label="Heat exchanger mass")
ax1.plot(ovs, ms, label="Total FC system mass")
ax2.plot(ovs, etas, "k--", label="FC system efficiency")
ax1.grid()
ax1.set_xlabel("Oversizing factor")
ax1.set_ylabel("Mass in $kg$")
ax2.set_ylabel("FC system efficiency")
# ax2.yaxis.label.set_color("g")
# ax2.spines["right"].set_edgecolor("g")
# ax2.tick_params(axis='y', colors="g")
# ax1.legend(title="Mass of")
fig.legend()
plt.tight_layout()
plt.show()
fig.savefig("figures/oversizing.pdf")

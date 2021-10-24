# Effect of payload on entire aircraft - not used in report
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

data = loadmat("../aircraft_sensitivities_data/payload_res.mat")
pls = data["pls_percent"][0]
m_bat = data["res_bat"][0]
m_fuel = data["res_fuel"][0]
m_mtom = data["res_mtom"][0]
m_oem = data["res_oem"][0]
m_fc = data["res_fc"][0]

plt.plot(pls, m_bat/m_bat[0], marker=4, label="Battery", markersize=10)
plt.plot(pls, m_fuel/m_fuel[0], marker=5, label="Fuel", markersize=10)
plt.plot(pls, m_mtom/m_mtom[0], marker=6, label="MTOM", markersize=10)
plt.plot(pls, (m_oem+m_bat)/(m_oem[0]+m_bat[0]), marker=7, label="OEM", markersize=10)
plt.plot(pls, m_fc/m_fc[0], marker="o", label="FC units", markersize=7)
plt.xlabel("Fraction of original payload")
plt.ylabel("Fraction of original mass")
plt.grid()
plt.legend()
plt.savefig("figures/payload_ac.pdf")
plt.show()

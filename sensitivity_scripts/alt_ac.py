# Effect of altitude on entire aircraft - Fig 4.15
import matplotlib.pyplot as plt
from scipy.io import loadmat

data = loadmat("../aircraft_sensitivities_data/alt_res.mat")
alts = data["alts"][0][:-1]/1e3
m_bat = data["res_bat"][0][:-1]
m_fuel = data["res_fuel"][0][:-1]
m_mtom = data["res_mtom"][0][:-1]
m_oem = data["res_oem"][0][:-1]
m_fc = data["res_fc"][0][:-1]

plt.plot(alts, m_bat/2882, marker=4, label="Battery", markersize=10)
plt.plot(alts, m_fuel/956, marker=5, label="Fuel", markersize=10)
plt.plot(alts, m_mtom/38870, marker=6, label="MTOM", markersize=10)
plt.plot(alts, (m_oem+m_bat)/31265, marker=7, label="OEM", markersize=10)
plt.plot(alts, m_fc*2/3949, marker="o", label="FC units", markersize=7)
plt.xlabel("Cruise altitude in km")
plt.ylabel("Fraction of original mass")
plt.grid()
plt.legend()
plt.savefig("figures/alt_ac.pdf")
plt.show()

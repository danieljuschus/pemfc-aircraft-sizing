# Effect of range on entire aircraft - Fig 4.12
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

data = loadmat("../aircraft_sensitivities_data/range_res.mat")
ranges = data["ranges_percent"][0]
m_bat = data["res_bat"][0]
m_fuel = data["res_fuel"][0]
m_mtom = data["res_mtom"][0]

mid = np.argwhere(ranges == 1)

plt.plot(ranges, m_bat/m_bat[mid][0], marker=4, label="Battery", markersize=10)
plt.plot(ranges, m_fuel/m_fuel[mid][0], marker=5, label="Fuel", markersize=10)
plt.plot(ranges, m_mtom/m_mtom[mid][0], marker=6, label="MTOM", markersize=10)
plt.xlabel("Fraction of original range")
plt.ylabel("Fraction of original mass")
plt.grid()
plt.legend()
plt.savefig("figures/range_ac.pdf")
plt.show()

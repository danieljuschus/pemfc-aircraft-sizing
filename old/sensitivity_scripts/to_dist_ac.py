# Effect of take-off distance on entire aircraft - Fig 4.14
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

# 1/2: final climb phi 0.3/0.35
data1 = loadmat("../aircraft_sensitivities_data/to_dist_res1.mat")
data2 = loadmat("../aircraft_sensitivities_data/to_dist_res2.mat")
dist = data1["dist"][0]
m_bat1 = data1["res_bat"][0]
m_fuel1 = data1["res_fuel"][0]
m_mtom1 = data1["res_mtom"][0]
m_bat2 = data2["res_bat"][0]
m_fuel2 = data2["res_fuel"][0]
m_mtom2 = data2["res_mtom"][0]


plt.plot(dist/1000, m_bat1/2882, marker="o", label=r"Battery, $\phi_{climb,end} = 0.3$", markersize=7)
plt.plot(dist/1000, m_fuel1/956, marker=5, label=r"Fuel, $\phi_{climb,end} = 0.3$", markersize=10)
plt.plot(dist/1000, m_mtom1/38870, marker=6, label=r"MTOM, $\phi_{climb,end} = 0.3$", markersize=10)
plt.plot(dist/1000, m_bat2/2882, marker="s", label=r"Battery, $\phi_{climb,end} = 0.35$", markersize=6)
plt.plot(dist/1000, m_fuel2/956, marker=5, label=r"Fuel, $\phi_{climb,end} = 0.3$5", markersize=10)
plt.plot(dist/1000, m_mtom2/38870, marker=6, label=r"MTOM, $\phi_{climb,end} = 0.35$", markersize=10)
plt.xlabel("Take-off distance in km")
plt.ylabel("Fraction of original mass")
plt.grid()
plt.legend()
plt.savefig("figures/to_dist_ac.pdf")
plt.show()

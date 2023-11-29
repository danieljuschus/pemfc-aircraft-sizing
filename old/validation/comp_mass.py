# Compressor mass validation - Fig 3.17
from compressor_performance import compressor_performance_model
from compressor_mass import compressor_mass_model
import numpy as np
import matplotlib.pyplot as plt

powers = np.linspace(1e4, 4.5e5, 10)
ms = []
p_comps = []
for po in powers:
    g, p, _, mdot = compressor_performance_model(po, 0.8, 1.75, 100000, 288, 0.0000181206, 1e5)
    #print(mdot)
    ms.append(compressor_mass_model(g, p))
    p_comps.append(p/1000)

ref_data = np.genfromtxt("compressor_data.csv", skip_header=1, usecols=(4, 5), delimiter=",").transpose()
ref_data = np.delete(ref_data, np.argwhere(np.isnan(ref_data))[:, 1], 1)
ref_data = ref_data[:, np.argsort(ref_data[0])]
pol = np.poly1d(np.polyfit(*ref_data, 1))

plt.plot(p_comps, ms, "-k", label="Own model")
plt.plot(*ref_data, "xk", label="Reference data")
plt.plot(ref_data[0], pol(ref_data[0]), "--k", label="Reference data linear fit")
plt.grid()
plt.legend()
plt.xlabel("Compressor power in kW")
plt.ylabel("Compressor mass in kg")
# plt.savefig("figures/comp_valid.pdf")
plt.show()

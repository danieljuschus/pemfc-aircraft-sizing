# Experimental voltage loss model for operation at high altitude - Fig 3.2
import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("../../../../volt_pres.csv", delimiter=",", skip_header=1)
pfs = np.linspace(data[:, 0][0], data[:, 0][-1], 50)
vfs = [-0.022830*pf**4 + 0.230982*pf**3 - 0.829603*pf**2 + 1.291515*pf + 0.329935 for pf in pfs]

plt.plot(data[:, 0], data[:, 2], "kX", label="Experimental data", markersize=10)
plt.plot(pfs, vfs, "k-", label="Polynomial fit")
plt.xlabel("Fraction of ISA sea level pressure")
plt.ylabel("Fraction of reference voltage")
plt.legend()
plt.grid()
plt.savefig("figures/voltage_pressure_exp_model.pdf")
plt.show()

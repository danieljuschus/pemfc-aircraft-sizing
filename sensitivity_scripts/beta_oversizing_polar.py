# Effect of changes in beta and oversizing factor on polarization curve - Fig 4.10
import matplotlib.pyplot as plt
import numpy as np
from stack_functions import cell_voltage, cell_model

pres1 = 54000  # 5 km
pres2 = 1.5*pres1
temp = 273.15 + 80

js = np.arange(0, 20000)
vs1 = [cell_voltage(ji, pres1, pres1, temp) for ji in js]
vs2 = [cell_voltage(ji, pres2, pres2, temp) for ji in js]
volt_cell_1, power_dens_cell_1, _ = cell_model(pres1, pres1, temp, oversizing=0.2)
volt_cell_2, power_dens_cell_2, _ = cell_model(pres1, pres1, temp, oversizing=0.4)
ps1 = js * vs1
ps2 = js * vs2

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(js, vs1, label=r"Voltage, $\beta = 0$")
ax1.plot(js, vs2, label=r"Voltage, $\beta = 1.5$")
ax2.plot(js, ps1, color="tab:blue", linestyle="--", label=r"Power density, $\beta = 0$")
ax2.plot(js, ps2, color="tab:orange", linestyle="--", label=r"Power density, $\beta = 1.5$")
ax2.plot(power_dens_cell_1 / volt_cell_1, power_dens_cell_1, "rX",
         label="Design point, oversizing factor of 0.2", ms=15, mew=0.001)
ax2.plot(power_dens_cell_2 / volt_cell_2, power_dens_cell_2, "r+",
         label="Design point, oversizing factor of 0.4", ms=15, mew=5)
ax1.grid()
ax1.set_xlabel("Current density in $A/m^2$")
ax1.set_ylabel("Cell voltage in $V$")
ax2.set_ylabel("Power density in $W/m^2$")
fig.legend(loc=(0.25, 0.15))
plt.tight_layout()
plt.savefig("figures/beta_oversizing_pol_curve.pdf")
plt.show()

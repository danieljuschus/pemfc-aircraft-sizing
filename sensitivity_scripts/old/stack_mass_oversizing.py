from stack_functions import stack_model, cell_model
from heat_exchanger import heat_exchanger_model
from ambiance import Atmosphere
import numpy as np
import matplotlib.pyplot as plt

pres = 1e5
temp = 273.15 + 80
power = 1e6
n_ss = 2
v_req = 800
mach_cr = 0.3
h_cr = 7000
atm_cr = Atmosphere(h_cr)
c_cr = atm_cr.speed_of_sound[0]  # cruise speed of sound in m
v_cr = mach_cr*c_cr
p_cr = atm_cr.pressure[0]  # static pressure at cruise altitude
p_cr_tot = p_cr*(1+0.4/2*mach_cr**2)**(1.4/0.4)  # total pressure at cruise in Pa
t_cr = atm_cr.temperature[0]  # static temperature at cruise altitude in K
t_cr_tot = t_cr * (1 + 0.4 / 2 * mach_cr ** 2)  # total temperature at cruise in K
mu_cr = atm_cr.dynamic_viscosity[0]  # dynamic viscosity at cruise altitude in Pa s
rho_cr = atm_cr.density[0]  # air density at cruise altitude in kg/m3

mu_f = 0.95

os = np.arange(0, 0.6, 0.1)

m_ss, m_hxs, etas = [], [], []
for o in os:
    volt_cell, power_dens_cell, eta = cell_model(pres, pres, temp, o)
    m_ss.append(stack_model(n_ss, v_req, volt_cell, power, power_dens_cell))
    m_hxs.append(heat_exchanger_model(power, volt_cell, temp, mu_f, v_cr, mach_cr, p_cr_tot, t_cr_tot, rho_cr, mu_cr))
    etas.append(eta)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(os, m_ss, label="Stacks")
ax1.plot(os, m_hxs, "r", label="Heat exchanger")
ax2.plot(os, etas, "g")
ax1.grid()
ax1.set_xlabel("Oversizing factor")
ax1.set_ylabel("Component mass in $kg$")
ax2.set_ylabel("Cell efficiency")
ax2.yaxis.label.set_color("g")
ax2.spines["right"].set_edgecolor("g")
ax2.tick_params(axis='y', colors="g")
ax1.legend(title="Component mass")
plt.tight_layout()
plt.show()
fig.savefig("figures/oversizing_comp_mass.pdf")

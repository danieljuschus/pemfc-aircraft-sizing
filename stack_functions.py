import numpy as np
from math import log, sqrt, exp, pi, isclose


def stack_model(n_stacks_series, volt_req, volt_cell, power_req, power_dens_cell):
    """
    Calculate mass of stack(s) in the FC unit.

    :param n_stacks_series: Number of stacks in series in FC unit
    :param volt_req: Voltage to be delivered by FC unit in V
    :param volt_cell: Nominal cell voltage in V
    :param power_req: Electrical power to be delivered by stacks (bigger than propulsive output power of FC unit)
    :param power_dens_cell: Nominal cell power density in W/m^2
    :return: mass of stack(s) in kg
    """
    # # constants
    # bipolar plate
    t_bp = 2e-4  # m
    rho_bp = 8e3  # kg/m3 - SS304L

    # endplate
    t_ep = 2.5e-2  # m
    rho_ep = 8e3  # kg/m3 - same as bp

    # bolts
    n_bolt = 10  # see Dey 2019
    rho_bolt = 8e3

    # MEA
    rho_mea = 0.2  # kg/m2 - Kadyk 2018

    # # calculations - per stack
    n_cells = volt_req / volt_cell / n_stacks_series  # number of cells per stack
    area_cell = power_req / n_stacks_series / n_cells / power_dens_cell  # area of single cell

    m_bp = n_cells*t_bp*rho_bp*area_cell  # mass of bipolar plates of one stack
    m_ep = 2*t_ep*rho_ep*area_cell  # mass of endplates of one stack
    m_mea = rho_mea*area_cell*n_cells  # mass of MEA of one stack

    d_bolt = sqrt(area_cell/600)  # diameter of bolts - see Dey 2019
    l_bolt = 2*t_ep + n_cells*t_bp  # length of bolts
    m_bolts = n_bolt*pi*(d_bolt/2)**2*l_bolt*rho_bolt  # mass of bolts of one stack

    m_tot = m_bp+m_ep+m_bolts+m_mea  # mass of a single stack
    # print("BP: {}, EP: {}, Bolts: {}, MEA: {}".format(m_bp, m_ep, m_bolts, m_mea    ))

    return m_tot*n_stacks_series


def mass_flow_stack(power_stack, volt_cell):
    """
    Compute mass flow of air required by fuel cell stack with given power output and cell voltage.

    :param power_stack: Electrical output power of stack in W
    :param volt_cell: Cell voltage in V
    :return: mass flow in kg/s
    """
    stoich = 2  # assumed stoichiometry
    return 3.58e-7 * stoich * power_stack / volt_cell


def cell_model(pres_air, pres_h, cell_temp, oversizing=0.1):
    """
    Find nominal operating point of cell, i.e. the pair of current density and voltage at which the cell is supposed to 
    operate at cruise.
    
    :param pres_air: Inlet air pressure in Pa
    :param pres_h: Inlet hydrogen pressure in Pa
    :param cell_temp: Operational temperature of cell in K
    :param oversizing: oversizing factor (e.g. 0.1 means that the maximum power of the cell lies 10 percent above
    the chosen operating point) - default 0.1
    :returns:
        - volt_cell: cell voltage at chosen point in V
        - power_dens_cell: power density at chosen point in W/m^2
        - eta_cell: cell efficiency
    """
    js = np.arange(0, 20000)  # array of possible current densities in A/m2
    ps = [cell_voltage(j, pres_air, pres_h, cell_temp) * j for j in js]  # power density for each j in W/m2
    j_op = (1 - oversizing) * js[np.nanargmax(ps)]  # current density of nominal operating point
    volt_cell = cell_voltage(j_op, pres_air, pres_h, cell_temp)  # cell voltage resulting from chosen point
    power_dens_cell = volt_cell*j_op  # power density at chosen point
    eta_cell = volt_cell / 1.482  # cell efficiency - HHV
    # eta_cell = volt_cell / 1.229  # cell efficiency - LHV
    return volt_cell, power_dens_cell, eta_cell


def cell_voltage(j, pres_air, pres_h, cell_temp):
    """
    Compute cell voltage from the current density (i.e. a certain operating point) and the pressure of the reactants.

    :param j: Current density in A/m^2
    :param pres_air: Inlet air pressure in Pa
    :param pres_h: Inlet hydrogen pressure in Pa
    :param cell_temp: Operational temperature of cell in K
    :returns:
        - Cell voltage in V
    """
    # constants - mostly based on textbook by O'Hayre - see report
    pot_rev0 = 1.229  # V - reversible potential at standard state for H2-O2 reaction - O'Hayre
    pres0 = 1  # atm - reference pressure
    temp0 = 289.15  # K - reference temperature
    farad = 96485.3329  # C/mol - Faraday constant
    r_gas = 8.314  # J/K/mol - gas constant
    r_air = 287.058  # J/kg/K - specific gas constant for air
    trans = 0.3  # transfer coefficient at cathode
    j_leak = 100  # A/m^2 - leakage current density
    res = 1e-6  # ohm/m^2 - area specific resistance
    mass_trans_c = 0.5  # V - mass transport loss constant
    j_lim = 2e4  # A/m^2 - limiting current density

    pot_rev = pot_rev0 - 44.34 / (2 * farad) * (cell_temp - temp0) + r_gas * cell_temp / (2 * farad) * log(
        pres_h * sqrt(pres_air*0.21))  # reversible potential, as function of operating temp and air pressure - O'Hayre

    # cathode exchange current density as function of temperature and pressure (Barbir 2012) - not used
    # compare with Dicks page 118
    # j_0 = 1 * (pres_air*0.21/101325)**1 * exp(-66e3/(r_gas*cell_temp)*(1-cell_temp/298.15))

    # constant exchange current density
    j_0 = 1  # O'Hayre

    # Determine voltage
    if j_lim - j - j_leak > 0:  # otherwise undefined
        volt = pot_rev - r_gas * cell_temp / (2 * trans * farad) * log((j + j_leak) / j_0) - res * j - \
               mass_trans_c * log(j_lim / (j_lim - j - j_leak))  # operational cell voltage, including all losses
    else:
        volt = np.nan
    if volt < 0:  # don't include negative voltages
        volt = np.nan

    if not isclose(pres_air, 101325):
        # correction in case cathode inlet pressure is not sea level - derived from Sondergaard and Pratt et al.
        pf = pres_air/101325
        volt = volt*(-0.022830*pf**4 + 0.230982*pf**3 - 0.829603*pf**2 + 1.291515*pf + 0.329935)

    return volt


def plot_polarization_curve():
    """
    Plot polarization curve with design point.
    """
    import matplotlib.pyplot as plt
    pres = 1e5
    temp = 273.15+80

    js = np.arange(0, 20000)
    vs = [cell_voltage(ji, pres, pres, temp) for ji in js]
    volt_cell, power_dens_cell, _ = cell_model(pres, pres, temp)
    ps = js*vs

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(js, vs, "k-", label="Voltage")
    ax2.plot(js, ps, "k--", label="Power density")
    ax2.plot(power_dens_cell/volt_cell, power_dens_cell, "rX", label="Design point", markersize=15)
    ax1.grid()
    ax1.set_xlabel("Current density in $A/m^2$")
    ax1.set_ylabel("Cell voltage in $V$")
    ax2.set_ylabel("Power density in $W/m^2$")
    fig.legend(loc=(0.25, 0.15))
    plt.tight_layout()
    plt.savefig("validation/figures/pol_curve.pdf")
    plt.show()


if __name__ == "__main__":
    plot_polarization_curve()

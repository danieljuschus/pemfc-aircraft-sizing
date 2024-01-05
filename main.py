# Sizing script for a fuel cell (FC) system

import sys
from ambiance import Atmosphere
from stack_functions import cell_model, stack_model, mass_flow_stack
from compressor_performance import compressor_performance_model
from compressor_mass import compressor_mass_model
from humidifier import humidifier_model
from heat_exchanger import heat_exchanger_model


def size_system(power_fc_sys, volt_req, h_cr, mach_cr, oversizing, beta, comp_bool, n_stacks_series):
    """
    Size a fuel cell system.

    :param n_stacks_series: Number of stacks in series
    :param power_fc_sys: Effective power to be delivered by the fuel cell system (in W)
    :param volt_req: Voltage to be delivered by the fuel cell system (in V)
    :param h_cr: Cruise altitude in m
    :param mach_cr: Cruise Mach number
    :param oversizing: Oversizing factor, 1-oversizing_factor = j_operating_point / j_max_power
    :param beta: Compression factor for cathode, i.e. a factor of 2 means that the cathode inlet pressure is twice the
    ambient pressure at the given altitude
    :param comp_bool: Boolean for whether to include a compressor
    :return: results: List of sizing results
    """
    # Atmospheric conditions
    atm_cr = Atmosphere(h_cr)
    c_cr = atm_cr.speed_of_sound[0]  # speed of sound at cruise in m
    v_cr = mach_cr * c_cr  # cruise true airspeed in m/s
    p_cr = atm_cr.pressure[0]  # static pressure at cruise altitude in Pa
    p_cr_tot = p_cr * (1 + 0.4 / 2 * mach_cr ** 2) ** (1.4 / 0.4)  # total pressure at cruise in Pa
    t_cr = atm_cr.temperature[0]  # static temperature at cruise altitude in K
    t_cr_tot = t_cr * (1 + 0.4 / 2 * mach_cr ** 2)  # total temperature at cruise in K
    rho_cr = atm_cr.density[0]  # air density at cruise altitude in kg/m3
    mu_cr = atm_cr.dynamic_viscosity[0]  # dynamic viscosity at cruise altitude in Pa s

    # Other inputs
    cell_temp = 273.15 + 80  # operating temperature inside cell
    mu_f = 0.95  # fuel utilisation

    # Compressor outlet conditions
    if comp_bool:
        pres_cathode_in = beta * p_cr_tot  # assuming the flow slows down completely before the compressor, see
        #                                  compressor_performance.py
    else:
        pres_cathode_in = p_cr_tot

    # Cell model
    pres_h = Atmosphere(0).pressure[0]  # assume that the anode inlet pressure is equal to sea level air pressure
    volt_cell, power_dens_cell, eta_cell = cell_model(pres_cathode_in, pres_h, cell_temp, oversizing)

    # Compressor models
    if comp_bool:
        # Iterate until the fuel cell stacks produce enough power for propulsion AND to run their compressor
        power_req = 0  # initiated as 0 to get iteration going
        power_req_new = power_fc_sys  # initially, the stacks only need to produce the propulsive power
        while abs(power_req_new - power_req) > 1e-3:  # while not converged within tolerance
            power_req = power_req_new  # this is the power produced by stacks
            geom_comp, power_comp, rho_humid_in, m_dot_comp = compressor_performance_model(power_req, volt_cell, beta,
                                                                                           p_cr_tot, t_cr_tot, mu_cr)
            power_req_new = power_fc_sys + power_comp  # (new) compressor power has been determined, add this to
            # propulsive
            #                                             power
        m_comp = compressor_mass_model(geom_comp, power_comp)  # determine compressor mass
    else:
        # no compressor
        m_comp = 0
        power_comp = 0
        power_req_new = power_fc_sys
        m_dot_comp = mass_flow_stack(power_req_new, volt_cell)  # mass flow of air for cathode in kg/s
        rho_humid_in = rho_cr  # humidifier inlet air density in kg/m3

    # Remaining BOP models
    m_humid = humidifier_model(m_dot_comp, rho_humid_in)  # mass of humidifier
    m_hx = heat_exchanger_model(power_req_new, volt_cell, cell_temp, mu_f, v_cr, mach_cr, p_cr_tot, t_cr_tot, rho_cr,
                                mu_cr)
    # mass of heat exchanger

    # Stack model
    m_stacks = stack_model(n_stacks_series, volt_req, volt_cell, power_req_new, power_dens_cell)  # mass of stack(s)

    # Sum up to find mass of FC system (all masses in kg)
    m_sys = m_stacks + m_comp + m_humid + m_hx
    print("Stack(s): {} kg, Compressor: {} kg, Humidifier: {} kg, Heat Exchanger: {} kg"
          .format(m_stacks, m_comp, m_humid, m_hx))
    # print("Power density of system in kW/kg: ", round(power_fc_sys/1000/m_sys, 3))
    # print("Stack prop output power: {} kW, Comp power: {} kW".format(power_fc_sys, power_comp))

    # Determine FC system efficiency
    eta_fcsys = eta_cell * power_fc_sys / (power_comp + power_fc_sys) * mu_f
    print("Cell efficiency: {}, Output efficiency: {}".format(eta_cell, eta_fcsys))

    # Make list of values to return
    results = [m_sys]

    return results


if __name__ == "__main__":
    # Run directly, for debugging
    power_fc_sys = 1e6
    volt_req = 500  # voltage to be produced by a fuel cell system
    h_cr = 5000
    mach_cr = 0.4
    oversizing = 0.2
    beta = 1.1
    comp_bool = 1
    n_stacks_series = 2

    sizing_results = size_system(power_fc_sys, volt_req, h_cr, mach_cr, oversizing, beta, comp_bool, n_stacks_series)

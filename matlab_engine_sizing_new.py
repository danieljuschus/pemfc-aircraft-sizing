# Sizing script for a fuel cell (FC) unit, used by the Class2WeightEstimation module of the Initiator (the Matlab
# aircraft sizing tool). The Initiator calls this when estimating the weight of the TurboGen engine of the serial hybrid
# configuration in getEngineWeight.m.
#
# This script can also be used directly, for debugging, sensitivity analyses or validation.
#
# One FC engine consists of a FC unit and an electric motor. One FC unit can consist of multiple stacks, but only has
# one set of BOP components. One stack consists of multiple cells.
#
# Important assumptions:
# - All engines produce same power (i.e. are identical)
# - Cruise is the sizing constraint for the fuel cell power loading (this needs to be ensured through the right
#   combinations of xi and phi in the missionSettings in the xml aircraft input files)
#
# Main inputs:
# 0: Propulsive power output of FC unit in W (the stacks need to produce more to run the BOP)
# 1: Output voltage in V (constant for now)
#
# Cruise condition inputs:
# 2: Cruise altitude in m
# 3: Cruise Mach number
#
# Inputs used only for sensitivity studies, validation and debugging:
# 4: Cell oversizing factor
# 5: Compressor pressure ratio (if equal to 0, then compressor is sized to provide sea level pressure)
# 6: Compressor boolean (False = no compressor)
#
# Outputs:
# 0: Mass of FC unit in kg
# 1: Cell efficiency

import sys
import streamlit as st
import pandas as pd
import plotly.express as px
from ambiance import Atmosphere
from stack_functions import cell_model, stack_model, mass_flow_stack
from compressor_performance import compressor_performance_model
from compressor_mass import compressor_mass_model
from humidifier import humidifier_model
from heat_exchanger import heat_exchanger_model
from scipy.io import savemat
import matplotlib.pyplot as plt


# Depending on number of inputs, determine if run by Matlab, by sensitivity analyses or directly (for debugging)
if len(sys.argv) == 1:
    # Run directly, for debugging
    power_fc_unit = 1e6
    volt_req = 500  # voltage to be produced by a fuel cell unit
    h_cr = 5000
    mach_cr = 0.4
    oversizing = 0.2
    beta = 0
    comp_bool = 1
elif len(sys.argv) > 4:
    # Run for sensitivity studies from Python script, so additional inputs are defined
    power_fc_unit = float(sys.argv[1])
    volt_req = 500
    h_cr = float(sys.argv[2])
    mach_cr = float(sys.argv[3])
    oversizing = float(sys.argv[4])
    beta = float(sys.argv[5])
    comp_bool = float(sys.argv[6])
elif len(sys.argv) == 4:
    # Run from Matlab
    power_fc_unit = float(sys.argv[1])
    volt_req = 500  # voltage to be produced by a fuel cell unit
    h_cr = float(sys.argv[2])
    mach_cr = float(sys.argv[3])
    oversizing = 0.2
    comp_bool = 1
else:
    raise TypeError("Wrong number of inputs for FC engine sizing python script.")

st.sidebar.title("User inputs")
st.sidebar.number_input("Required output power from system in MW", key="power", 
                value=power_fc_unit/1e6, step=0.1, min_value=0.1, max_value=2.)
st.sidebar.number_input("Cruise altitude in km", key="altitude", 
                value=h_cr/1e3, step=0.1, min_value=2., max_value=10.)

power_fc_unit = float(st.session_state.power)*1e6
h_cr = float(st.session_state.altitude)*1e3

# Atmospheric conditions
atm_cr = Atmosphere(h_cr)
c_cr = atm_cr.speed_of_sound[0]  # speed of sound at cruise in m
v_cr = mach_cr*c_cr  # cruise true airspeed in m/s
p_cr = atm_cr.pressure[0]  # static pressure at cruise altitude in Pa
p_cr_tot = p_cr*(1+0.4/2*mach_cr**2)**(1.4/0.4)  # total pressure at cruise in Pa
t_cr = atm_cr.temperature[0]  # static temperature at cruise altitude in K
t_cr_tot = t_cr * (1 + 0.4 / 2 * mach_cr ** 2)  # total temperature at cruise in K
rho_cr = atm_cr.density[0]  # air density at cruise altitude in kg/m3
mu_cr = atm_cr.dynamic_viscosity[0]  # dynamic viscosity at cruise altitude in Pa s

# Configuration of FC unit
n_stacks_series = 2  # number of stacks in series

# Other inputs
if len(sys.argv) == 4:
    # if run from Matlab, beta needs to be determined
    beta = Atmosphere(0).pressure[0]/p_cr_tot   # compression ratio at cruise level (pressure back to sea level)
elif not beta:
    # if beta is defined as 0, either from debugging or sensitivity study, redefine beta as usual
    # otherwise, if beta already has another value, then don't do anything
    beta = Atmosphere(0).pressure[0]/p_cr_tot  # compression ratio at cruise level (pressure back to sea level)
cell_temp = 273.15+80  # operating temperature inside cell
mu_f = 0.95  # fuel utilisation

# Compressor outlet conditions
if comp_bool:
    pres_cathode_in = beta*p_cr_tot  # assuming the flow slows down completely before the compressor, see
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
    power_req_new = power_fc_unit  # initially, the stacks only need to produce the propulsive power
    while abs(power_req_new-power_req) > 1e-3:  # while not converged within tolerance
        power_req = power_req_new  # this is the power produced by stacks
        geom_comp, power_comp, rho_humid_in, m_dot_comp = compressor_performance_model(power_req, volt_cell, beta,
                                                                                       p_cr_tot, t_cr_tot, mu_cr)
        power_req_new = power_fc_unit + power_comp  # (new) compressor power has been determined, add this to propulsive
        #                                             power
    m_comp = compressor_mass_model(geom_comp, power_comp)  # determine compressor mass
else:
    # no compressor
    m_comp = 0
    power_comp = 0
    power_req_new = power_fc_unit
    m_dot_comp = mass_flow_stack(power_req_new, volt_cell)  # mass flow of air for cathode in kg/s
    rho_humid_in = rho_cr  # humidifier inlet air density in kg/m3

# Remaining BOP models
m_humid = humidifier_model(m_dot_comp, rho_humid_in)  # mass of humidifier
m_hx = heat_exchanger_model(power_req_new, volt_cell, cell_temp, mu_f, v_cr, mach_cr, p_cr_tot, t_cr_tot, rho_cr, mu_cr)
# mass of heat exchanger

# Stack model
m_stacks = stack_model(n_stacks_series, volt_req, volt_cell, power_req_new, power_dens_cell)  # mass of stack(s)

# Sum up to find mass of FC unit (all masses in kg)
m_unit = m_stacks + m_comp + m_humid + m_hx
print("Stack(s): {} kg, Compressor: {} kg, Humidifier: {} kg, Heat Exchanger: {} kg"
      .format(m_stacks, m_comp, m_humid, m_hx))
# print("Power density of unit in kW/kg: ", round(power_fc_unit/1000/m_unit, 3))
# print("Stack prop output power: {} kW, Comp power: {} kW".format(power_fc_unit, power_comp))

# Determine FC unit efficiency
eta_fcsys = eta_cell*power_fc_unit/(power_comp+power_fc_unit)*mu_f
#print("Cell efficiency: {}, Output efficiency: {}".format(eta_cell, eta_fcsys))

# Write mass of engine and cell efficiency to output file
if len(sys.argv) == 4:
    # for Matlab
    savemat("External/FuelCell/output_from_python.mat", {"engineMass": m_unit, "engineEfficiency": eta_fcsys})
elif len(sys.argv) > 4:
    # for other purposes
    import hickle
    hickle.dump([m_stacks, m_comp, m_humid, m_hx, eta_fcsys], "sensitivity_output.hkl")
    
col1, col2, col3 = st.columns(3)
with col1:
    st.header("Numerical results")
    st.write("System mass {} kg".format(round(m_unit,2)))
    st.write("Gravimetric power density of system: {} kW/kg".format(round(power_fc_unit/1e3/m_unit,2)))
    
with col2: 
    mass_df = pd.DataFrame({"Component": ["Stack(s)", "Compressor", "Humidifier", "Heat exchanger"],
                           "Mass": [m_stacks, m_comp, m_humid, m_hx]})
    fig = px.pie(mass_df, values="Mass", names="Component")
    st.write(fig)
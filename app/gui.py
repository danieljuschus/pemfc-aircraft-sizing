# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 16:34:41 2022

@author: jusc_da
"""

import streamlit as st
import pandas as pd
import plotly.express as px
from ambiance import Atmosphere
from models.stack_functions import cell_model, stack_model, mass_flow_stack
from models.compressor_performance import compressor_performance_model
from models.compressor_mass import compressor_mass_model
from models.humidifier import humidifier_model
from models.heat_exchanger import heat_exchanger_model


def size_fc_system(h_cr, mach_cr, comp_bool, oversizing, volt_req, power_fc_unit):
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

    beta = Atmosphere(0).pressure[0]/p_cr_tot   # compression ratio at cruise level (pressure back to sea level)

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
    #print("Stack(s): {} kg, Compressor: {} kg, Humidifier: {} kg, Heat Exchanger: {} kg"
    #      .format(m_stacks, m_comp, m_humid, m_hx))
    # print("Power density of unit in kW/kg: ", round(power_fc_unit/1000/m_unit, 3))
    # print("Stack prop output power: {} kW, Comp power: {} kW".format(power_fc_unit, power_comp))

    # Determine FC unit efficiency
    eta_fcsys = eta_cell*power_fc_unit/(power_comp+power_fc_unit)*mu_f
    #print("Cell efficiency: {}, Output efficiency: {}".format(eta_cell, eta_fcsys))
    
    return m_unit, m_stacks, m_comp, m_humid, m_hx, eta_fcsys
    

# m = size_fc_system(5000, 0.4, 1, 0.2, 400, 1e6)

st.sidebar.title("User inputs")
form = st.sidebar.form(key="input")
form.number_input("Cruise altitude in km", key="altitude", 
                value=5., step=0.1, min_value=2., max_value=10.)
form.number_input("Mach number", key="mach", 
                value=0.4, step=0.05, min_value=0.05, max_value=0.65)
comp_bool = form.checkbox("Compressor", value=1)
form.number_input("Current oversizing factor", key="oversizing", 
                value=0.2, step=0.05, min_value=0., max_value=0.5)
form.number_input("Required output voltage in V", key="voltage", 
                value=400., step=50., min_value=100., max_value=850.)
form.number_input("Required output power from system in MW", key="power", 
                value=1., step=0.1, min_value=0.1, max_value=2.)
form.form_submit_button("Run")


h_cr = st.session_state.altitude*1e3
mach_cr = st.session_state.mach
oversizing = st.session_state.oversizing
voltage = st.session_state.voltage
power_fc_unit = st.session_state.power*1e6

m_unit, m_stacks, m_comp, m_humid, m_hx, eta_fcsys=  \
size_fc_system(h_cr, mach_cr, comp_bool, oversizing, voltage, power_fc_unit)

col1, col2, col3 = st.columns(3)
with col1:
    st.header("Numerical results")
    st.write("System mass: {} kg".format(round(m_unit,2)))
    st.write("Gravimetric power density of system: {} kW/kg".format(round(power_fc_unit/1e3/m_unit,2)))
    st.write("System efficiency: {}".format(round(eta_fcsys,2)))
    
with col2: 
    st.header("Component mass breakdown")
    mass_df = pd.DataFrame({"Component": ["Stack(s)", "Compressor", "Humidifier", "Heat exchanger"],
                           "Mass": [m_stacks, m_comp, m_humid, m_hx]})
    fig = px.pie(mass_df, values="Mass", names="Component")
    st.write(fig)
    
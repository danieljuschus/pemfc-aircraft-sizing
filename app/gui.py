# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 16:34:41 2022

@author: jusc_da
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import pyvista as pv
from ambiance import Atmosphere
from models.stack_functions import cell_model, stack_model, mass_flow_stack
from models.compressor_performance import compressor_performance_model
from models.compressor_mass import compressor_mass_model
from models.humidifier import humidifier_model
from models.heat_exchanger import heat_exchanger_model

import cadquery as cq
from cadquery import exporters, importers

import numpy as np


def size_fc_system(h_cr, mach_cr, comp_bool, oversizing, volt_req, power_fc_unit, beta, n_stacks_series):
    figs = []
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
    #n_stacks_series = 2  # number of stacks in series

    #beta = Atmosphere(0).pressure[0]/p_cr_tot   # compression ratio at cruise level (pressure back to sea level)

    cell_temp = 273.15+80  # operating temperature inside cell
    mu_f = 0.95  # fuel utilisation

    # Compressor outlet conditions
    # if comp_bool:
    #     pres_cathode_in = beta*p_cr_tot  # assuming the flow slows down completely before the compressor, see
    #     #                                  compressor_performance.py
    # else:
    #     pres_cathode_in = p_cr_tot
   
    pres_cathode_in = p_cr_tot*beta

    # Cell model
    pres_h = Atmosphere(0).pressure[0]  # assume that the anode inlet pressure is equal to sea level air pressure
    volt_cell, power_dens_cell, eta_cell, fig = cell_model(pres_cathode_in, pres_h, cell_temp, oversizing)
    figs.append(fig)

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
    m_hx, dim_hx = heat_exchanger_model(power_req_new, volt_cell, cell_temp, mu_f, v_cr, mach_cr, p_cr_tot, t_cr_tot, rho_cr, mu_cr)
    # mass of heat exchanger

    # Stack model
    m_stacks, dim_stack, res_stack = stack_model(n_stacks_series, volt_req, volt_cell, power_req_new, power_dens_cell)  # mass of stack(s)

    #res_stack{"n_cells": }

    # Sum up to find mass of FC unit (all masses in kg)
    m_unit = m_stacks + m_comp + m_humid + m_hx
    #print("Stack(s): {} kg, Compressor: {} kg, Humidifier: {} kg, Heat Exchanger: {} kg"
    #      .format(m_stacks, m_comp, m_humid, m_hx))
    # print("Power density of unit in kW/kg: ", round(power_fc_unit/1000/m_unit, 3))
    # print("Stack prop output power: {} kW, Comp power: {} kW".format(power_fc_unit, power_comp))

    # Determine FC unit efficiency
    eta_fcsys = eta_cell*power_fc_unit/(power_comp+power_fc_unit)*mu_f
    #print("Cell efficiency: {}, Output efficiency: {}".format(eta_cell, eta_fcsys))
    
    # Hydrogen comsumption
    mdot_h2 = 1.05e-8 * (power_comp+power_fc_unit)/volt_cell
    
    return m_unit, m_stacks, m_comp, m_humid, m_hx, eta_fcsys, mdot_h2, power_comp, figs, dim_stack, n_stacks_series, dim_hx, res_stack

    
def generate_system_geometry():
    n_stack_series = 2
    gap = 0.1
    exporters.export(cq.Workplane("front").box(1.5-gap, 2, 0.5).translate((-0.75, 0, 0.5+gap)), 
                     "media/stack_1.stl") 
    exporters.export(cq.Workplane("front").box(1.5-gap, 2, 0.5).translate((0.75, 0, 0.5+gap)),
                     "media/stack_2.stl")
    
    exporters.export(cq.Workplane("front").box(1.5, 1.5, 0.5).translate((0.3,0,2*(0.5+gap))),
                     "media/humidifier.stl")

    exporters.export(cq.Workplane("front").box(3, 1.5, 0.5), 
                     "media/hx.stl")
    
    exporters.export(cq.Workplane("left").circle(0.25).extrude(0.75).translate((-0.5,0,2*(0.5+gap))), 
                     "media/compressor.stl")
    
    exporters.export(cq.Workplane("left").circle(0.75).extrude(0.4).translate((-1.5,0,0.5+gap)), 
                     "media/em.stl")
    
    pts = np.genfromtxt("media/airfoil.dat")*3.5
    #exporters.export(cq.Workplane("top").circle(0.1).extrude(0.1).translate((0,1,0)), "media/nacelle.stl")
    #exporters.export(cq.Workplane("top").circle(0.1).extrude(0.1).translate((0,1,0)), "media/nacelle.stl")
    
    #exporters.export(cq.Workplane("top").circle(0.1).revolve(360, (0,1,0.5+gap), (1,1,0.5+gap)), "media/nacelle.stl")
    exporters.export(cq.Workplane("top").spline(pts).close().revolve(360, (-1,-1,0), (1,-1,0)).translate((-2,0,-0.5)), "media/nacelle.stl")
    #exporters.export(cq.Workplane("top").spline(pts).close().extrude(0.1), "media/nacelle.stl")
    #exporters.export(cq.Workplane("top").spline(pts).close().extrude(0.1).translate((0,0,0.5+gap)), "media/nacelle.stl")
    
    #prop = importers.importStep("media/Prop2.stp")
    #print(type(prop))
    #exporters.export(prop, "media/prop.stl")

    
#     res = (
#     cq.Assembly(boxes[0])
#     .add(boxes[1])
# )
    
    #res =(cq.Assembly().add(boxes[0]).add(boxes[1]))
    #exporters.export(res.toCompound(), "media/subsys{}.stl".format(i))


h_init = 5000.
mach_init = 0.4

st.set_page_config(layout="wide")
st.sidebar.title("User inputs")
form = st.sidebar.form(key="input")
form.number_input("Cruise altitude in km", key="altitude", 
                value=h_init/1e3, step=0.1, min_value=2., max_value=10.)
form.number_input("Mach number", key="mach", 
                value=mach_init, step=0.05, min_value=0.05, max_value=0.65)

h_cr = st.session_state.altitude*1e3

mach_cr = st.session_state.mach

atm_cr = Atmosphere(h_cr)
atm_init = Atmosphere(h_init)

p_cr_tot = atm_cr.pressure[0]*(1+0.4/2*mach_cr**2)**(1.4/0.4)
p_cr_tot_init = atm_init.pressure[0]*(1+0.4/2*mach_init**2)**(1.4/0.4)
beta_sl = 101325/p_cr_tot
beta_sl_init = 101325/p_cr_tot_init

col1, col2, col3 = form.columns(3)
with col1:
    comp_bool = st.checkbox("Compressor", value=1)
with col2:
    comp_sl_bool = st.checkbox("Provide sea level pressure to the compressor",
                               disabled=not comp_bool, value=1)
with col3: 
    st.number_input("Compressor pressure ratio", 
                    key="beta",value=beta_sl_init, step=0.05, 
                    min_value=0., max_value=3.5, 
                    disabled=not comp_bool or (comp_bool and comp_sl_bool))
    
    if not comp_bool:
        beta = 1
    elif comp_sl_bool:
        beta = beta_sl
    else:
        beta = st.session_state.beta
        
    st.write("Current pressure ratio: {}".format(round(beta,2)))
    


form.number_input("Current oversizing factor", key="oversizing", 
                value=0.2, step=0.05, min_value=0., max_value=0.5)
form.number_input("Number of stacks in series", key="n_stacks_series",
                  value=2, min_value=1, max_value=10)
form.number_input("Required output voltage in V", key="voltage", 
                value=400., step=50., min_value=100., max_value=850.)
form.number_input("Required output power from system in MW", key="power", 
                value=1., step=0.1, min_value=0.1, max_value=2.)
form.form_submit_button("Run")

oversizing = st.session_state.oversizing
n_stacks_series = st.session_state.n_stacks_series
voltage = st.session_state.voltage
power_fc_unit = st.session_state.power*1e6

m_unit, m_stacks, m_comp, m_humid, m_hx, eta_fcsys, mdot_h2, power_comp, figs,\
    dim_stack, n_stacks_series, dim_hx, res_stack =  \
size_fc_system(h_cr, mach_cr, comp_bool, oversizing, voltage, power_fc_unit, beta, n_stacks_series)

generate_system_geometry()    

col1, col2 = st.columns(2)
with col1:
    st.header("Numerical results")

    # hydrogen mass flow
    df = pd.DataFrame(data={'Value': [str(round(m_unit,1)) + " kg", 
                                      str(round(power_fc_unit/1e3/m_unit,1)) + " kW/kg",
                                      str(round(eta_fcsys*100,1)) + " %",
                                      str(round(mdot_h2*1000,1)) + " g/s",
                                      str(round(power_comp/1e3,1)) + " kW"]}, 
                      index=["System mass", 
                             "Gravimetric power density of system",
                             "System efficiency",
                             "Hydrogen consumption",
                             "Compressor power"],)
    st.table(df)
    
    df = pd.DataFrame(data={"Value": [round(dim_stack[0],1),
                                      round(dim_stack[2],1),
                                      round(res_stack[0],1),
                                      round(res_stack[1],1)]},
                      index=["Width/height in m", 
                             "Length in m", 
                             "Number of cells",
                             "Area per cell in m²"])
    
    st.warning('Numerical values below are just placeholders', icon="⚠️")
    
    with st.expander("Stack"):
        st.table(df)

    df = pd.DataFrame(data={"Value": [0,0]},
                      index=["Width/height in m", 
                             "Length in m"])

    with st.expander("Humidifier"):
        st.table(df)

    with st.expander("Heat exchanger"):
        st.table(df)
    
    with st.expander("Compressor"):
        st.table(df)

    st.header("System geometry")
    
    df = pd.DataFrame(data={"Color": ["Red", "Green", "Blue", "Yellow", "Orange"]},
                      index=["Heat exchanger", "Stack", "Humidifier", "Compressor", "EM"])
    with st.expander("Legend"): 
        st.warning('Legend in 3D-model viewer not yet working, replaced with table for now', icon="⚠️")
        st.table(df)
        
    plotter = pv.Plotter(
         window_size=[400,400]) 
    plotter.background_color = "white"
    
    reader = pv.STLReader("media/hx.stl")
    
    # ## Read data and send to plotter
    mesh = reader.read()
    plotter.add_mesh(mesh, color="red", label="Heat exchanger")

    reader = pv.STLReader("media/stack_1.stl")
   
    # ## Read data and send to plotter
    mesh = reader.read()
    plotter.add_mesh(mesh, color="green", label="Stack")
    
    reader = pv.STLReader("media/stack_2.stl")
   
    # ## Read data and send to plotter
    mesh = reader.read()
    plotter.add_mesh(mesh, color="green", label="Stack")
    
    reader = pv.STLReader("media/humidifier.stl")
   
    # ## Read data and send to plotter
    mesh = reader.read()
    plotter.add_mesh(mesh, color="blue", label="Stack")
    
    reader = pv.STLReader("media/compressor.stl")
   
    # ## Read data and send to plotter
    mesh = reader.read()
    plotter.add_mesh(mesh, color="yellow", label="Stack")
    
    reader = pv.STLReader("media/nacelle.stl")
   
    # ## Read data and send to plotter
    mesh = reader.read()
    plotter.add_mesh(mesh, label="Nacelle",opacity=0.5)
    
    reader = pv.STLReader("media/em.stl")
   
    # ## Read data and send to plotter
    mesh = reader.read()
    plotter.add_mesh(mesh, label="EM",color="orange")
    
    # reader = pv.STLReader("media/prop.stl")
   
    # # ## Read data and send to plotter
    # mesh = reader.read()
    # plotter.add_mesh(mesh)

    ## Export to an external pythreejs
    model_html = "model.html"
    other = plotter.export_html(model_html, backend='pythreejs')

    ## Read the exported model
    with open(model_html,'r') as file: 
        model = file.read()

    # ## Show in webpage
    st.warning('Not the actual geometry, just an example', icon="⚠️")
    st.components.v1.html(model,height=500)

    st.header("Compressor geometry")
    
    reader = pv.STLReader("media/comp.stl")
    plotter = pv.Plotter(
        window_size=[400,400]) 
    plotter.background_color = "white"
    #plotter.add_ruler()

    ## Read data and send to plotter
    mesh = reader.read()
    plotter.add_mesh(mesh)

    ## Export to an external pythreejs
    model_html = "model.html"
    other = plotter.export_html(model_html, backend='pythreejs')

    ## Read the exported model
    with open(model_html,'r') as file: 
        model = file.read()

    ## Show in webpage
    st.components.v1.html(model,height=500)
    
with col2: 
    st.header("Component mass breakdown")
    mass_df = pd.DataFrame({"Component": ["Stack(s)", "Compressor", "Humidifier", "Heat exchanger"],
                           "Mass": [m_stacks, m_comp, m_humid, m_hx]})
    fig = px.pie(mass_df, values="Mass", names="Component", color="Component",
                 color_discrete_map={'Stack(s)':'lightcyan',
                                 'Compressor':'cyan',
                                 'Humidifier':'royalblue',
                                 'Heat exchanger':'darkblue'})
    # make colors constant
    st.plotly_chart(fig, use_container_width=True)
    st.header("Polarization curve")
    figs[0].update_layout(height=800)
    st.plotly_chart(figs[0], use_container_width=True)
    
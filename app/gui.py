# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 16:34:41 2022

@author: jusc_da
"""

# GUI script
# first cd into app, then run in terminal with command "streamlit run gui.py"

from models.stack_functions import cell_model, stack_model, mass_flow_stack
from models.compressor_performance import compressor_performance_model
from models.compressor_mass import compressor_mass_model
from models.humidifier import humidifier_model
from models.heat_exchanger import heat_exchanger_model
from main import size_system

import cadquery as cq
from cadquery import exporters, importers
import streamlit as st
import pandas as pd
import plotly.express as px
import pyvista as pv
from stpyvista import stpyvista
from ambiance import Atmosphere
import numpy as np


def generate_system_geometry():
    """
    Generates the geometry for the GUI.

    Note: currently fixed dummy dimensions are used, that's why this function has no inputs yet.

    The resulting geometries are saved individually as stl files in the media directory.
    """
    n_stack_series = 2
    gap = 0.1
    exporters.export(cq.Workplane("front").box(1.5 - gap, 2, 0.5).translate((-0.75, 0, 0.5 + gap)),
                     "media/stack_1.stl")
    exporters.export(cq.Workplane("front").box(1.5 - gap, 2, 0.5).translate((0.75, 0, 0.5 + gap)),
                     "media/stack_2.stl")

    exporters.export(cq.Workplane("front").box(1.5, 1.5, 0.5).translate((0.3, 0, 2 * (0.5 + gap))),
                     "media/humidifier.stl")

    exporters.export(cq.Workplane("front").box(3, 1.5, 0.5),
                     "media/hx.stl")

    exporters.export(cq.Workplane("left").circle(0.25).extrude(0.75).translate((-0.5, 0, 2 * (0.5 + gap))),
                     "media/compressor.stl")

    exporters.export(cq.Workplane("left").circle(0.75).extrude(0.4).translate((-1.5, 0, 0.5 + gap)),
                     "media/em.stl")

    pts = np.genfromtxt("media/airfoil.dat") * 3.5
    exporters.export(
        cq.Workplane("top").spline(pts).close().revolve(360, (-1, -1, 0), (1, -1, 0)).translate((-2, 0, -0.5)),
        "media/nacelle.stl")

    # prop = importers.importStep("media/Prop2.stp")
    # exporters.export(prop, "media/prop.stl")


# ## Start of the GUI script

# # Sidebar with inputs
h_init = 5000.
mach_init = 0.4

st.set_page_config(layout="wide")
st.sidebar.title("User inputs")
form = st.sidebar.form(key="input")
form.number_input("Cruise altitude in km", key="altitude",
                  value=h_init / 1e3, step=0.1, min_value=2., max_value=10.)
form.number_input("Mach number", key="mach",
                  value=mach_init, step=0.05, min_value=0.05, max_value=0.65)

h_cr = st.session_state.altitude * 1e3

mach_cr = st.session_state.mach

atm_cr = Atmosphere(h_cr)
atm_init = Atmosphere(h_init)

p_cr_tot = atm_cr.pressure[0] * (1 + 0.4 / 2 * mach_cr ** 2) ** (1.4 / 0.4)
p_cr_tot_init = atm_init.pressure[0] * (1 + 0.4 / 2 * mach_init ** 2) ** (1.4 / 0.4)
beta_sl = 101325 / p_cr_tot
beta_sl_init = 101325 / p_cr_tot_init

col1, col2, col3 = form.columns(3)
with col1:
    comp_bool = st.checkbox("Compressor", value=1)
with col2:
    comp_sl_bool = st.checkbox("Provide sea level pressure to the compressor",
                               disabled=not comp_bool, value=1)
with col3:
    st.number_input("Compressor pressure ratio",
                    key="beta", value=beta_sl_init, step=0.05,
                    min_value=0., max_value=3.5,
                    disabled=not comp_bool or (comp_bool and comp_sl_bool))

    if not comp_bool:
        beta = 1
    elif comp_sl_bool:
        beta = beta_sl
    else:
        beta = st.session_state.beta

    st.write("Current pressure ratio: {}".format(round(beta, 2)))

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
volt_req = st.session_state.voltage
power_fc_sys = st.session_state.power * 1e6

# # Run sizing to generate results
m_sys, m_stacks, m_comp, m_humid, m_hx, eta_fcsys, mdot_h2, power_comp, figs, \
    dim_stack, n_stacks_series, dim_hx, res_stack = \
    size_system(power_fc_sys, volt_req, h_cr, mach_cr, oversizing, beta, comp_bool, n_stacks_series)

generate_system_geometry()

# # Display results
col1, col2 = st.columns(2)
with col1:
    st.header("Numerical results")

    df = pd.DataFrame(data={'Value': [str(round(m_sys, 1)) + " kg",
                                      str(round(power_fc_sys / 1e3 / m_sys, 1)) + " kW/kg",
                                      str(round(eta_fcsys * 100, 1)) + " %",
                                      str(round(mdot_h2 * 1000, 1)) + " g/s",
                                      str(round(power_comp / 1e3, 1)) + " kW"]},
                      index=["System mass",
                             "Gravimetric power density of system",
                             "System efficiency",
                             "Hydrogen consumption",
                             "Compressor power"], )
    st.table(df)

    df = pd.DataFrame(data={"Value": [round(dim_stack[0], 1),
                                      round(dim_stack[2], 1),
                                      round(res_stack[0], 1),
                                      round(res_stack[1], 1)]},
                      index=["Width/height in m",
                             "Length in m",
                             "Number of cells",
                             "Area per cell in m²"])

    st.warning('Numerical values below are just placeholders', icon="⚠️")

    with st.expander("Stack"):
        st.table(df)

    df = pd.DataFrame(data={"Value": [0, 0]},
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

    plotter = pv.Plotter(window_size=[400, 400])
    plotter.background_color = "#0e1117"

    reader = pv.STLReader("media/hx.stl")
    mesh = reader.read()
    plotter.add_mesh(mesh, color="red", label="Heat exchanger")

    reader = pv.STLReader("media/stack_1.stl")
    mesh = reader.read()
    plotter.add_mesh(mesh, color="green", label="Stack")

    reader = pv.STLReader("media/stack_2.stl")
    mesh = reader.read()
    plotter.add_mesh(mesh, color="green", label="Stack")

    reader = pv.STLReader("media/humidifier.stl")
    mesh = reader.read()
    plotter.add_mesh(mesh, color="blue", label="Stack")

    reader = pv.STLReader("media/compressor.stl")
    mesh = reader.read()
    plotter.add_mesh(mesh, color="yellow", label="Stack")

    reader = pv.STLReader("media/nacelle.stl")
    mesh = reader.read()
    plotter.add_mesh(mesh, label="Nacelle", opacity=0.5)

    reader = pv.STLReader("media/em.stl")
    mesh = reader.read()
    plotter.add_mesh(mesh, label="EM", color="orange")
    plotter.view_isometric()

    # plotter.show()

    # reader = pv.STLReader("media/prop.stl")
    # mesh = reader.read()
    # plotter.add_mesh(mesh)

    # ## Show in webpage
    st.warning('Not the actual geometry, just an example', icon="⚠️")
    stpyvista(plotter, key="engine")

    st.header("Compressor geometry")

    reader = pv.STLReader("media/comp.stl")
    plotter = pv.Plotter(
        window_size=[400, 400])
    plotter.background_color = "#0e1117"

    ## Read data and send to plotter
    mesh = reader.read()
    plotter.add_mesh(mesh)
    plotter.view_isometric()

    # ## Show in webpage
    stpyvista(plotter, key="compressor")

with col2:
    st.header("Component mass breakdown")
    mass_df = pd.DataFrame({"Component": ["Stack(s)", "Compressor", "Humidifier", "Heat exchanger"],
                            "Mass": [m_stacks, m_comp, m_humid, m_hx]})
    fig = px.pie(mass_df, values="Mass", names="Component", color="Component",
                 color_discrete_map={'Stack(s)': 'lightcyan',
                                     'Compressor': 'cyan',
                                     'Humidifier': 'royalblue',
                                     'Heat exchanger': 'darkblue'})
    # make colors constant
    st.plotly_chart(fig, use_container_width=True)
    st.header("Polarization curve")
    figs[0].update_layout(height=800)
    st.plotly_chart(figs[0], use_container_width=True)

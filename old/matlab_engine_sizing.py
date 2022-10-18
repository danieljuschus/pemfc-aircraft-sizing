# Sizing script for a fuel cell unit, used by the Class2WeightEstimation module of the Initiator. The Initiator
# calls this when estimating the weight of the TurboGen engine of the serial hybrid configuration in getEngineWeight.m.
#
# To save time, no new engine type was created in the Initiator. Instead, code was injected at several places ("if
# HydrogenPropulsionType == "FC"") to replace the GT/TurboGen of a serial hybrid configuration with a fuel cell unit.
#
# Due to the way that the hybrid code works, the TurboGen engine (fuel cell unit in this case) and the ElectroProp (EM1
# + prop) are implemented as separate Parts of the Aircraft. However, for the purpose of this thesis, the fuel cell
# engine will be a single physical component on the wing, housing both.
#
# Nomenclature used in this script:
# - One FC stack consists of multiple cells
# - One FC unit consists either of one stack and its BOP, or multiple stacks in series and ONE set of BOP components -
#   the latter is useful to simplify the structural design of the stacks, only the effect on the weight is taken into
#   account
# - One FC engine consists of a FC unit and the electric motor and propeller - unit and engine are sometimes used
#   interchangeably in this script, but not in the thesis report
#
# Important assumptions:
# - all engines produce same power (i.e. are identical)
# - one engine can have multiple stacks in series, but it only has one of each BOP component. all stacks are identical.
# - cruise is the sizing constraint for the fuel cell power loading (this needs to be ensured through the right
#   combinations of xi and phi in the missionSettings in the xml aircraft input files)
#
# Inputs from Matlab:
# 1: power per FC unit in W - this is the "propulsive" power output that the Initiator requires from the FC unit (the 
#                             stacks need to produce more to also run the BOP)
# 2: number of engines (i.e. number of FC units in parallel)
# 3: cruise altitude in m
# 4: cruise mach number
#
# Outputs:
# 1: mass of FC unit in kg
# 2: cell efficiency

import sys
from ambiance import Atmosphere
from stack_functions import cell_operating_point, mass_flow_stack, stack_mass
from compressor_performance import compressor_performance_model
from compressor_mass import compressor_mass_model
from humidifier import humidifier_model
from heat_exchanger import heat_exchanger_model
from scipy.io import savemat

# Inputs from Matlab
power_fc_unit = float(sys.argv[1])  # power output of FC unit in W
n_engines = float(sys.argv[2])  # number of engines in parallel
h_cr = float(sys.argv[3])  # cruise altitude in m
mach_cr = float(sys.argv[4])  # cruise mach number

# Atmospheric/cruise conditions
atm_cr = Atmosphere(h_cr)
c_cr = atm_cr.speed_of_sound[0]  # cruise speed of sound in m
v_cr = mach_cr*c_cr
p_cr = atm_cr.pressure[0]  # static pressure at cruise altitude
p_cr_tot = p_cr*(1+0.4/2*mach_cr**2)**(1.4/0.4)  # total pressure at cruise in Pa
t_cr = atm_cr.temperature[0]  # static temperature at cruise altitude in K
t_cr_tot = t_cr * (1 + 0.4 / 2 * mach_cr ** 2)  # total temperature at cruise in K
rho_cr = atm_cr.density[0]  # air density at cruise altitude in kg/m3
mu_cr = atm_cr.dynamic_viscosity[0]  # dynamic viscosity at cruise altitude in Pa s

# Configuration of FC unit
n_stacks_series = 2  # number of stacks in series

# Other constant inputs
volt_req = 500  # voltage to be produced by a fuel cell unit
comp_rpm = 10000  # compressor RPM
beta = 2  # compression ratio at cruise level - roughly back to sea level
cell_temp = 273.15+80  # operating temperature inside cell
mu_f = 0.95  # fuel utilisation - assume recirculation, i.e. this doesn't increase required amount of hydrogen

# Compressor outlet conditions
pres_comp_out = p_cr*beta  # assuming the flow slows down completely before the compressor, see
#                            compressor_performance.py

# Find operating point and number of cells
pres_h = pres_comp_out  # assume that the anode inlet pressure is the same as at the cathode TODO: is this reasonable?
volt_cell, power_dens_cell = cell_operating_point(pres_comp_out, pres_h, cell_temp)
eta_cell = volt_cell/1.229  # cell efficiency
n_cells = volt_req/volt_cell  # number of cells in FC unit

# Compressor sizing
# Iterate until the fuel cell stacks produce enough power for propulsion AND to run their compressor
power_req = 0  # initiated as 0 to get iteration going
power_req_new = power_fc_unit  # initially, the stacks only need to produce the propulsive power
while abs(power_req_new-power_req) > 1e-3:
    power_req = power_req_new
    m_dot_stack = mass_flow_stack(power_req/n_stacks_series, volt_cell)
    geom_comp, power_comp, rho_humid_in = compressor_performance_model(beta, m_dot_stack * n_stacks_series, p_cr_tot,
                                                                       t_cr_tot, mu_cr, comp_rpm)
    power_req_new = power_fc_unit + power_comp
power_req = power_req
m_comp = compressor_mass_model(geom_comp, power_comp)
    
# Size other BOP components and find stack mass
m_humid = humidifier_model(m_dot_stack * n_stacks_series, rho_humid_in)  # mass of single humidifier in kg
m_hx = heat_exchanger_model(power_req, volt_cell, cell_temp, mu_f, v_cr, mach_cr, p_cr, t_cr, rho_cr, mu_cr)
area_cell = power_req/n_cells/power_dens_cell  # area of single cell
m_stack = stack_mass(area_cell, n_cells/n_stacks_series)  # mass of single stack in kg

# Sum up to find size of FC unit
m_unit = m_stack*n_stacks_series + m_comp + m_humid + m_hx
#print("Power density of unit in kW/kg: ", round(power_fc_unit/1000/m_unit, 3))

# Write mass of engine and cell efficiency to output file for Matlab
savemat("External/FuelCell/output_from_python.mat", {"engineMass": m_unit, "engineEfficiency": eta_cell})
#with open("output_to_matlab.txt", "w") as file:
#    file.write(str(m_unit) + "\n" + str(eta_cell))

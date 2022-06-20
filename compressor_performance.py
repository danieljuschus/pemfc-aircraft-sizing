import pyromat as pm
import warnings
import os
from math import sin, cos, acos, tan, atan, sqrt, pi, radians, degrees, log, log10, exp
from stack_functions import mass_flow_stack

# use standard SI units
pm.config["unit_energy"], pm.config["unit_molar"], pm.config["unit_pressure"] = "J", "mol", "Pa"
# create air ideal gas mixture
air = pm.get("ig.air")
# warnings.simplefilter("always")

r_air = pm.units.const_Ru/(air.mw()/1000000)
# Wrong order of magnitude for weight from air.mw: 28964 g/kmol --> Corrected Value by dividing with 1e6: 0,0289 g/mol


def visc(t):
    """
    Determine dynamic viscosity of air from temperature.

    :param t: Temperature in K
    :return: Dynamic viscosity in kg/m/s
    """
    return 1.458e-6 * t**1.5 / (t + 110.4)


def p_s_theo(s_in, t_in):
    """
    Determine pressure from entropy and temperature. Not included in Pyromat yet.

    :param s_in: Entropy in J/K
    :param t_in: Temperature in K
    :return: Pressure in Pa
    """
    return exp((air.s(t_in, 1)[0]-s_in)/r_air)


def compressor_performance_model(power_req, volt_cell, beta, p1, t1, mu1, n=10000):
    """
    Determine compressor performance by performing sizing process. Method as presented in Gambini & Vellini.
    Overall isentropic efficiency assumed constant. See report.

    Assumption: flow slows down to almost zero velocity at compressor inlet --> total conditions of freestream become
    static input

    :param power_req: Electrical power to be produced by stack(s) in W
    :param volt_cell: Cell voltage in V
    :param beta: Compression ratio
    :param p1: Static inlet pressure in Pa (i.e. total pressure of freestream)
    :param t1: Static inlet temperature in K (i.e. total temperature of freestream)
    :param mu1: Dynamic viscosity in kg/m/s
    :param n: RPM
    :return:
        - geom: List of resulting dimensions of compressor, in m (d1h, d2, d2s, d3, la, b1, b2, b3)
        - p_comp: Power required to run compressor, in W
        - rho3: Compressor outlet air density in kg/m^3
        - m: mass flow through compressor in kg/s
    """
    # # Inputs
    # Mass flow
    m = mass_flow_stack(power_req, volt_cell)

    # Some preliminary things
    gamma = 1.4
    p3 = beta*p1  # outlet pressure
    del_h_is = air.cp(t1)[0] * t1 * (beta**((gamma-1)/gamma)-1)
    # r_air = pm.units.const_Ru/air.mw()

    # Assumed input parameters
    psi = 0.9  # work coefficient
    dt = 0.7  # rotor tip diameter ratio
    dh = dt * 0.7  # rotor hub diameter ratio
    alpha2 = radians(60)  # rotor outlet flow angle
    eta_is = 0.75  # isentropic efficiency
    eta_r = 0.95  # rotor efficiency
    alpha1 = 0  # inlet flow angle

    # # Kinematics
    # Directly from assumptions
    psi_is = eta_is*psi  # isentropic work coefficient
    u2 = sqrt(del_h_is/psi_is)  # blade speed at rotor outlet
    u1 = u2*dt  # blade speed at rotor inlet
    c2u = u2*psi  # assuming c1u = 0
    c2m = c2u/tan(alpha2)  # rotor outlet meridional velocity
    d2 = 60/pi*u2/n  # rotor outlet diameter
    rho1 = air.d(t1, p1)[0]  # rotor inlet density
    v1 = m/rho1  # rotor inlet volumetric flow rate
    b1 = d2/2 * (dt - dh)  # blade height at rotor inlet
    c1m = 2*v1 / (pi * d2 * (dt+dh) * b1)  # rotor inlet meridional velocity
    c1 = c1m * cos(alpha1)
    phi = c1m/u2  # flow coefficient
    xi = 1/(phi*tan(alpha2)) * (psi + phi*dt*tan(alpha1))  # rotor meridional velocity ratio
    # r = 1 - psi/2 + phi**2/(2*psi) * ((1-xi**2) + (tan(alpha1))**2 * (1-dt**2)) - phi * dt * tan(alpha1)
    beta1 = atan(dt/phi - tan(alpha1))  # rotor inlet flow angle
    beta2 = atan(1/(phi*xi) * (1-psi) - dt/xi*tan(alpha1))  # rotor outlet flow angle
    w1 = u2 * phi * sqrt(1 + (tan(beta1))**2)
    w2 = u2 * xi * phi * sqrt(1 + (tan(beta2))**2)
    c2 = u2 * xi * phi * sqrt(1 + (tan(alpha2))**2)

    # # Thermodynamics
    s1 = air.s(t1, p1)[0]
    h1 = air.h(t1)[0]
    i1 = h1 + w1**2/2 - u1**2/2
    i2 = i1  # rothalpy conservation in rotor
    h2 = i2 - w2**2/2 + u2**2/2
    h2_is = h1 + eta_r*(h2-h1)
    # limit the rotor outlet enthalpy, otherwise things go wrong
    if h2_is > 500000:
        # if output file still exists from previous run, remove it - otherwise it looks like it actually converged
        if os.path.isfile("External/FuelCell/output_from_python.mat"):
            os.remove("External/FuelCell/output_from_python.mat")
        raise ValueError("Rotor outlet enthalpy too high. Reduce RPM!")
    t2_is = float(air.T_h(h2_is))
    t2 = float(air.T_h(h2))
    p2 = p_s_theo(s1, t2_is)
    rho2 = air.d(t2, p2)[0]
    cs2 = sqrt(gamma*r_air*t2)
    t3_is = float(air.T_s(s1, p3))
    h3_is = air.h(t3_is)[0]
    h3 = h1 + (h3_is - h1)/eta_is
    t3 = float(air.T_h(h3))
    rho3 = air.d(t3, p3)[0]

    # # Geometry
    b2 = m / (rho2 * c2m * pi * d2)  # rotor outlet blade height
    # dh = 0.35  # rotor hub diameter ratio
    # dt = sqrt(dh**2 + 4*m/(rho1 * c1m * pi * d2**2))  # rotor tip diameter ratio
    d1h, d1t = dh*d2, dt*d2  # rotor inlet hub and tip diameter
    d1m = (d1h + d1t)/2  # rotor inlet mean diameter

    # iteratively determine beta2b, sf, nbr
    tau_a = 0.05 * b2  # rotor exit clearance
    tb = 0.003 * d2  # blade thickness - unshrouded impeller
    tan_alpha1m = tan(alpha1) * d1t/d1m  # inlet mean absolute flow angle
    beta1m = atan(dt/phi * d1m/d1t - tan_alpha1m)  # inlet mean relative flow angle
    w1m = c1m / cos(beta1m)  # mean relative velocity
    sf_s = sin(radians(19) + 0.2*(pi/2 - beta2))
    sf = 0.1
    sf_new = 1
    nbr, beta2b = 1, 1
    while abs(sf-sf_new)/sf > 0.01:
        sf = sf_new
        beta2b = atan(1/(xi*phi) - tan(alpha2/sf))
        betam = (beta1m + beta2b)/2
        nbr = int((2*pi*cos(betam))/(0.4*log(1/dt)))
        sf_int = 1 - sqrt(cos(beta2b))/nbr**0.7
        sf_new = (sf + sf_int)/2
        # if d1m / d2 > (sf - sf_s) / (1 - sf_s):
        #    warnings.warn("Rotor mean diameter ratio too high")
        # print(sf_new)

    s1r = pi * (d1t + d1h) / (2 * nbr)  # rotor inlet blade pitch
    s2r = pi * d2 / nbr  # rotor outlet blade pitch
    o1, o2 = s1r*cos(beta1m), s2r*cos(beta2b)
    d_hyd = o1*b1/(o1+b1) + o2*b2/(o2+b2)
    la_alt = d2 * (0.014 + 0.023/dh + 1.58*(dt**2 - dh**2)*phi)
    la = (d2-d1t)/2 + b2
    l_mr = pi/8 * (2*la - b2 + d2 - d1t + b1)
    l_hyd = l_mr/cos((beta1m+beta2b)/2)

    if alpha2 > radians(72):
        alpha2s = radians(72) + (alpha2 - radians(72))/4
    else:
        alpha2s = radians(72)
    ma2 = c2/cs2
    d2s = d2 * (1 + (radians(90) - alpha2s)/(2*pi) + ma2**2/15)
    c2su = c2u*d2/d2s
    c2sm = c2su/tan(alpha2s)
    c2s = c2sm/cos(alpha2s)
    b2s = m/(rho2*c2sm*pi*d2s)  # TODO: check correct rho
    if b2s > b2:
        b2s = b2
    l_hyd_vaneless = (d2s-d2)/2
    d_hyd_vaneless = b2 + b2s
    d3 = d2 * (1.55 + (dt**2 - dh**2)*phi)
    b3 = b2s
    c3m = m/(rho3*pi*d3*b3)
    if c3m < c1:
        c3 = c1
        alpha3 = acos(c3m/c3)
    else:
        c3 = c3m  # TODO: check if stuff needs to be recalculated here
        alpha3 = 0
    c3u = c3 * sin(alpha3)
    # if 10<nbr<20:
    #     nbs = nbr - 1
    # else:
    #     warnings.warn("Number of rotor blades too high")
    nbs = nbr-1
    s2s = pi*d2s/nbs
    s3 = pi*d3/nbs
    o2s = s2s*cos(alpha2s)
    o3 = s3*cos(alpha3)
    l_hyd_vaned = (d3-d2s)/(2*cos((alpha2s+alpha3)/2))
    d_hyd_vaned = o2s*b2s/(o2s + b2s) + o3*b3/(o3 + b3)

    # applicability check
    theta_c = pi * (d3*cos(alpha3) - d2s*cos(alpha2s)) / (2*nbs*l_hyd_vaned)
    bl = pi * (d2s*c2su - d3*c3u) / (nbs*l_hyd_vaned*(c2s-c3))
    ar = d3*b3*cos(alpha3)/(d2s*b2s*cos(alpha2s))
    # if not 7<2*degrees(theta_c)<11 or not 0<bl<1/3 or not 1.4<ar<2.4:
    #    warnings.warn("Check applicability of kinematics and geometry")
    # print(2*degrees(theta_c), bl, ar)

    geom = [d1h, d2, d2s, d3, la, b1, b2, b3]
    p_comp = (h3-h1)*m

    return geom, p_comp, rho3, m


if __name__ == "__main__":
    geom, power_comp, rho3, m_dot = compressor_performance_model(1e5, 0.9, 1.5, 40000, 250, 1.6e-5, 1e4)

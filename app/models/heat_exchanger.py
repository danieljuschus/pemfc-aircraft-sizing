from math import sqrt
import numpy as np

# constants:
r_air = 287.05
c_p = 1003.5
gamma = 1.4


def heat_exchanger_model(power_stacks, volt_cell, cell_temp, mu_f, v_cr, m_cr, p_cr, t_cr, rho_cr, mu_cr):
    """ Determine mass and dimensions of heat exchanger. Based on Kožulović 2020.

    Drag/ram effect not yet added. See report

    :param power_stacks: Power output of the stacks in W
    :param volt_cell: Cell voltage in V
    :param cell_temp: Cell temperature in K
    :param mu_f: Fuel utilisation
    :param v_cr: Cruise TAS in m/s
    :param m_cr: Cruise mach number
    :param p_cr: Cruise static pressure in Pa
    :param t_cr: Cruise static temperature in K
    :param rho_cr: Cruise density in kg/m^3
    :param mu_cr: Cruise dynamic viscosity in Pa s
    :return: Mass of heat exchanger in kg"""
    # design constants
    rho_al = 2710
    t_al = 0.2e-3
    lambda_al = 220

    q_all = power_stacks * (1.482 / (volt_cell * mu_f) - 1)  # heat release of engine in W - HHV value - Dicks
    del_t_all = cell_temp - t_cr

    v1_v0 = 0.15
    v_1 = v1_v0*v_cr
    f_t_split = 0.6
    del_t_fluid = del_t_all*f_t_split
    del_t_loc = del_t_all - del_t_fluid
    h_pass = 3e-3  # ?
    peri = (2 + 2*sqrt(2))*h_pass  # perimeter of tube
    displ = 0.776
    v_b = v_1/displ

    pt0 = p_cr*(1+(gamma-1)/2*m_cr**2)**(gamma/(gamma-1))
    tt0 = t_cr*(1+(gamma-1)/2*m_cr**2)
    t1 = tt0 - v_1**2/(2*c_p)
    rho_1 = 0.96*pt0/(r_air*t1 + 0.5*v_1**2)
    rho_11 = 0.96*p_cr/(0.4/1.4*(v_cr**2-v_1**2)/2 + p_cr/rho_cr)
    rho_b = rho_1
    a_pass_tube = h_pass**2
    m_dot_tube = rho_b * a_pass_tube * v_b
    q_tube = m_dot_tube*c_p*del_t_fluid
    n_tube = q_all/q_tube  # total number of tubes

    pr_h2o = 4.5
    rho_h2o = 1000
    c_p_h2o = 4182
    nu_h2o = 6.49
    lambda_h2o = 0.63
    height_h2o = 0.2e-3
    width_h2o = 1.6e-3
    d_h_h2o = 2*height_h2o*width_h2o/(height_h2o+width_h2o)
    alpha_h2o = nu_h2o*lambda_h2o/d_h_h2o
    pr_air = 0.71
    d_h_air = 4*h_pass**2/peri
    re_d_air = v_b*rho_b*d_h_air/mu_cr
    nu_air = 0.022*sqrt(pr_air)*re_d_air**0.8  # assuming always turbulent
    st_air = nu_air/(pr_air*re_d_air)
    alpha_air = st_air*v_b*rho_b*c_p
    alpha = (1/alpha_h2o + t_al/lambda_al + 1/alpha_air)**-1
    l_tube = 1.2 * q_tube / (alpha * peri * del_t_loc)

    # Resulting masses
    m_al = n_tube * (peri/2+h_pass) * t_al * l_tube * rho_al  # mass of aluminium
    m_h2o = n_tube * h_pass * height_h2o * l_tube * rho_h2o  # mass of coolant

    # Resulting dimensions - not used yet
    w_hx = (int(sqrt(n_tube))+1)*h_pass
    h_hx = int(sqrt(n_tube))*h_pass/displ

    return 1.2*(m_al + 2*m_h2o), [(np.ceil(np.sqrt(n_tube))+1)*h_pass,
                                  (np.ceil(np.sqrt(n_tube))+1)*h_pass,
                                  l_tube]


if __name__ == "__main__":
    t_cr = 216.65
    rho_cr = 0.363918
    mu_cr = 0.0000143226
    v_cr = 337/3.6  # cruise TAS - m/s
    m_cr = 0.8
    p_cr = 22632.1
    del_t_all = 110  # total temperature difference - K
    m_hx = heat_exchanger_model(11.46e6, 110, v_cr, m_cr, p_cr, t_cr, rho_cr, mu_cr)

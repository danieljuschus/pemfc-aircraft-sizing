from math import sqrt, ceil


def humidifier_model(mass_flow_air, density_air):
    """Determine mass of humidifier. Based on Huizing et al. 2008.

    Internal drag/pressure losses not yet accounted for. See report.

    :param mass_flow_air: Inlet mass flow of air through dry side (from compressor)
    :param density_air: Air density of inlet mass flow through through dry side (from compressor)
    :return: Mass of humidifier in kg"""
    dwa = 0.3e-4  # diffusivity of water in air in m**2/s - dependence on temperature ignored
    w = 5e-3  # width of channels in m
    h = 1.5e-3  # height of channels in m
    volum_air = mass_flow_air/density_air  # volumetric flow

    l = 4.5*h**2/dwa  # length of channels in m
    n = 2*volum_air/(1.5*w*h)  # total number of channels

    t_al = 0.2e-3
    t_mem = 0.2e-3
    rho_al = 2710
    rho_mem = 970  # UHMWPE

    m = ceil(sqrt(n))  # number of plates = number of channels per plate -> roughly square cross section
    return 1.2*l*m*(m+1) * (h*t_al*rho_al + w*t_mem*rho_mem)


if __name__ == "__main__":
    mass_flow_air = 0.8
    density_air = 1

    m_humid = humidifier_model(mass_flow_air, density_air)

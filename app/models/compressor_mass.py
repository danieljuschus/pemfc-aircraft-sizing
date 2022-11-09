import cadquery as cq
from cadquery import exporters

rho_ss = 8e3  # kg/m3 - SS304L


def compressor_mass_model(geom, power):
    """
    Determine compressor weight from geometry and power required. For details see report.

    :param geom: List of dimensions of compressor, in m
    :param power: Power consumption of the compressor (i.e. power the the electric motor has to provide to produce the
    required compression ratio) in W
    :return: mass of compressor (compressor body + electric motor + bearing assembply) in kg
    """
    # # unpack and scale up dimensions (the latter is required for Cadquery to work properl)
    d1h, d2, d2s, d3, la, b1, b2, b3 = [g*100 for g in geom]

    # # compressor body
    # backplate
    d4 = 1.1*d3
    t_mat = 0.01*d4
    t_bp = t_mat
    backplate = cq.Workplane("XY").circle(d4/2).extrude(t_bp)

    # impeller
    t_imp_b = t_mat
    imp = cq.Workplane("XZ").center(0, t_bp).line(0, la+t_imp_b).line(d1h/2, 0)\
        .ellipseArc((d2-d1h)/2, la, 180, 270)\
        .line(0, -t_imp_b).line(-d2/2, 0).close().revolve()

    # casing
    d_scr = 1.3*(d3-d2s)/2
    b_e = b3*0.1
    cas_i = cq.Workplane("XZ").center((d2s + d3)/4, t_bp+b_e+d_scr/2).circle(d_scr/2-t_mat)\
        .revolve(360,[-(d2s + d3)/4, 0], [-(d2s + d3)/4, 1])
    # cas_ic = cas_i.copyWorkplane(cq.Workplane("XY", origin=(0, 0, d_scr/4)))\
    #    .split(keepTop=True)
    cas_o = cq.Workplane("XZ").center((d2s + d3)/4, t_bp+b_e+d_scr/2).circle(d_scr/2)\
        .revolve(360,[-(d2s + d3)/4, 0], [-(d2s + d3)/4, 1])
    cas_h = cas_o.cut(cas_i)
    cas_ai = cq.Workplane("XZ").center(d1h/2 + b1, t_bp+t_imp_b+la)\
        .line((d2-d1h)/2, 0).line(0, -la+b2)\
        .ellipseArc((d2-d1h)/2-b1, la-b2, 180, 270, sense=-1).close()\
        .revolve(360, [-d1h/2 - b1, 0], [-d1h/2 - b1, 1])
    cas_ai_2 = cas_ai.cut(cas_o)
    out = cq.Workplane("XZ").center((d2s + d3)/4, t_bp+b_e+d_scr/2).circle(d_scr/2)\
        .extrude(d4/2)
    out_h = out.faces(">Y or <Y").shell(-t_mat)
    cas_ao = cq.Workplane("XZ").center((d2s + d3)/4+d_scr/2, t_bp+b_e+d_scr/2)\
        .line(0,-b_e-d_scr/2).line(-d_scr/2*0.1, 0).line(0,b_e+d_scr/2).close()\
        .revolve(360, [-(d2s + d3)/4-d_scr/2, 0], [-(d2s + d3)/4-d_scr/2, 1])
    cas_ao_c = cas_ao.cut(out)
    out_hc = out_h.cut(cas_i)
    cas_res3 = cas_h.copyWorkplane(cq.Workplane("XY", origin=(0, 0, d_scr/4)))\
        .split(keepTop=True)
    cas_res3_c = cas_res3.cut(out)

    # all parts of compressor body
    parts = [backplate, imp, cas_ai_2, cas_ao_c, out_hc, cas_res3_c]
    #exporters.export(cas_res3_c, 'app/media/comp.json', tolerance=0.01, angularTolerance=0.1,
    #                 exportType=exporters.ExportTypes.TJS)
    
    res = backplate
    for part in parts:
        res = res.union(part)
    
    #res = backplate.union(imp).union(cas_res3_c)
    exporters.export(res, "media/comp.stl")

    # sum up the material volumes of all parts, scale down and mulitply with material density
    m_cu = sum([part.objects[0].Volume() for part in parts]) / 100**3 * rho_ss

    # # EM
    m_em = power/2200  # based on statistical EM power density

    # # Total, including bearing assembly
    frac_bearings = 1.15  # from celeroton ct-17-700
    m_tot = frac_bearings*(m_cu+m_em)

    return m_tot


if __name__ == "__main__":
    d1h = 0.062
    d2 = 0.192
    d2s = 0.205
    d3 = 0.307
    la = 0.054
    b1 = 0.017
    b2 = 0.007
    b3 = 0.007

    geom = [d1h, d2, d2s, d3, la, b1, b2, b3]
    m_tot = compressor_mass_model(geom, 5000)

# BOP mass distribution pie charts - Fig 3.16
import matplotlib.pyplot as plt
import os
import hickle

os.system("/home/daniel/anaconda3/envs/thesis_code/bin/python3.8 ../matlab_engine_sizing_new.py 1e6 5000 0.3 0.2 0")
res = hickle.load("sensitivity_output.hkl")[1:-1]
# m_s, m_c, m_h, m_hx, eta = hickle.load("sensitivity_output.hkl")

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']

plt.pie(res,
        autopct="%.1f %%", pctdistance=1.25, radius=0.7, colors=colors)
plt.legend(["Compressor", "Humidifier", "Heat exchanger"], loc='upper center', bbox_to_anchor=(0.5, 0.1))
plt.tight_layout()
plt.savefig("figures/pie_sys_own.pdf", transparent=True, bbox_inches='tight', pad_inches=0.1)
plt.show()
plt.close()

plt.pie([46.6, 23.6, 100-46.6-23.6],
        autopct="%.1f %%", pctdistance=1.25, radius=0.7, colors=[colors[0], colors[2], colors[3]])
plt.legend(["Compressors", "Heat exchanger", "Other"], loc='upper center', bbox_to_anchor=(0.5, 0.1))
plt.tight_layout()
plt.savefig("figures/pie_sys_bradley.pdf", transparent=True, bbox_inches='tight', pad_inches=0.1)
plt.show()
plt.close()

# Propulsion system pie chart - not used in report
import matplotlib.pyplot as plt

res = [4008, 3424, 490, 2793, 201, 1351]
labels = ["FC systems", "Electric motors", "Nacelles", "Battery", "Cables", "Tank"]

plt.pie(res, labels=labels, radius=0.7)
# plt.legend(labels, loc='upper center', bbox_to_anchor=(0.5, 0.1))
plt.tight_layout()
plt.savefig("figures/pie_prop_sys.pdf", transparent=True, bbox_inches='tight', pad_inches=0.1)
plt.show()
# plt.close()

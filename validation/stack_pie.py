# Stack mass distribution pie charts - Fig 3.15
from stack_functions import cell_model, stack_model
import matplotlib.pyplot as plt

volt, pdens, _ = cell_model(101325, 101325, 350, 0.2)
n_cells = 72
m = stack_model(1, n_cells*volt, volt, 1e6, pdens)
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']

plt.pie([263.6, 119.6, 5, 32.9], autopct="%.1f %%", pctdistance=1.25, radius=0.7, colors=colors)
plt.legend(["Bipolar plates", "Endplates", "Bolts", "MEA"], loc='upper center', bbox_to_anchor=(0.5, 0.1))
plt.tight_layout()
# plt.savefig("figures/pie_own.pdf", transparent=True, bbox_inches='tight', pad_inches=0.1)
plt.show()
plt.close()

# fig2, ax2 = plt.subplots()
# ax2.pie([78.7, 11.7, 9.6], labels=["Bipolar plates", "Endplates & bolts", "MEA"], autopct="%.1f %%")
# ax2.axis("equal")
#
# fig3, ax3 = plt.subplots()
# ax3.pie([53, 19, 17, 3.9, 7.1], labels=["Bipolar plates", "Endplates", "Bolts", "MEA & gaskets", "Other"], autopct="%.1f %%")
# ax3.axis("equal")

plt.show()

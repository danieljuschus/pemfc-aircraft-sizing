# Effect of system voltage and number of stacks in system on stack mass - Fig. 4.11
from stack_functions import stack_model
import matplotlib.pyplot as plt
import numpy as np

volts = np.linspace(1e2, 1e3)

ms_1 = [stack_model(1, vi, 0.9, 1e6, 6000) for vi in volts]
ms_2 = [stack_model(2, vi, 0.9, 1e6, 6000) for vi in volts]

plt.plot(volts, ms_1, "k-", label="Single stack")
plt.plot(volts, ms_2, "k--", label="Two stacks in series")
plt.xlabel("System voltage in V")
plt.ylabel("Stack mass in kg")
plt.grid()
plt.legend()
plt.savefig("figures/sp_voltages.pdf")
plt.show()

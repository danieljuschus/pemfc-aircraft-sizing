from stack_functions import cell_model, stack_model
import numpy as np
import matplotlib.pyplot as plt

pres = 54019.9
betas = np.linspace(1.5, 2.5, 10)
powers = np.linspace(1045000, 1105000, len(betas))
ms = []
for i in range(len(betas)):
    volt, pdens, _ = cell_model(betas[i]*pres, 1e5, 350, 0.2)
    m = stack_model(2, 500, volt, powers[i], pdens)
    ms.append(m)
    # print(betai, volt, pdens)

plt.plot(betas, ms)
plt.show()

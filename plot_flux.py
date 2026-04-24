import matplotlib.pyplot as plt
import numpy as np


data = np.loadtxt('output/heat_flux_kappa_avg.txt')

# matplotlib.rcParams.update({'font.size': 24})
plt.figure(figsize=(12, 8))
plt.plot(data[:, 0], data[:, 3])
plt.xlabel('Шаг симуляции')
plt.ylabel('$\\kappa_{2D}$, Вт/К')
plt.title('Оценка электронной теплопроводности')
plt.grid()
plt.show()

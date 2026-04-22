import matplotlib.pyplot as plt
import numpy as np


flux = np.loadtxt('output/energy_flux_avg.txt')

# matplotlib.rcParams.update({'font.size': 24})
plt.figure(figsize=(12, 8))
plt.plot(flux)
plt.title('Поток энергии (не удельный)')
plt.grid()
plt.show()
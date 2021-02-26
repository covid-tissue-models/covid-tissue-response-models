import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.grid'] = True

ODEnames = ['Time','Virus','Viability','P','IFNe','STATP','IRF7','IRF7P','IFN']
ODE = np.genfromtxt('JordanOriginalODE_1.txt', skip_header=1, delimiter=',', names=ODEnames)
CC3Dnames = ['mcs', 'Time', 'IFN', 'STATP', 'IRF7', 'IRF7P', 'Viability', 'Virus', 'IFNe']
CC3D = np.genfromtxt('ifn_data.dat', skip_header=1, delimiter=',', names=CC3Dnames)

plot_names = ['STATP','IRF7','IRF7P','Virus','IFN','IFNe','Viability']
plot_colors = ['#0033FF','#FF6600','#666666','#6600FF','#CC0033','#993300','#339900']

fig, axs = plt.subplots(3,3)
counter = 0
for i in range(3):
    for j in range(3):
        if i == 2 and j == 0:
            axs[i,j].set_visible(False)
        elif i == 2 and j == 2:
            axs[i, j].set_visible(False)
        else:
            axs[i, j].plot(ODE['Time'], ODE[plot_names[counter]], color=plot_colors[counter], linestyle='dotted')
            axs[i, j].plot(CC3D['Time'], CC3D[plot_names[counter]], color=plot_colors[counter])
            axs[i, j].set_title(plot_names[counter])
            counter += 1
plt.tight_layout()
plt.savefig('validation.pdf',transparent=True)
plt.show()
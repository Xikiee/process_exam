import numpy as np
import matplotlib.pyplot as plt
import control as ctrl 
import sympy as sp

#### plotting an equation
def plotting_equation(t, eq):
    plt.plot(t,eq, label = 'equation')
    plt.xlabel('Time (s)')
    plt.ylabel('y(t)')
    plt.grid(which = 'both', linewidth = 0.5)
    plt.legend()
    plt.show()

t = np.linspace(0,50,10000)
eq = 22.8*0.4907*((1/12)*np.exp(-6*t) - (1/8)*np.exp(-4*t) + 1/24)
plotting_equation(t,eq)


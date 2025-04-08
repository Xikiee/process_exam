import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import control as ctrl

# #partial fraction decomposition

# s = sp.symbols('s', real = True)

# num = 5 * sp.exp(-0.5 * s)
# num = 1
# den = s * (s**2 + 10*s +24)
# eq = num/den
# partial_fractions = sp.apart(eq)
# print(partial_fractions)


t = np.linspace(0,50,10000)
eq = 22.8*0.4907*((1/12)*np.exp(-6*t) - (1/8)*np.exp(-4*t) + 1/24)
plt.plot(t,eq)
plt.show()


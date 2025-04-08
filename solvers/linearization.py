import numpy as np
import matplotlib.pyplot as plt
import control as ctrl 
import sympy as sp

#### linearization of a system
x = sp.symbols('x', real= True)
eq = 4*x**0.5
x0 = 1



import sympy as sp
def  linearization_of_system(eq,x0):

# Define the symbol
x = sp.Symbol('x')
xs = sp.Symbol('x_s')  # Operating point

# Define the function
F = 4 * sp.sqrt(x)

# First order Taylor series expansion around x_s
F_taylor = F.subs(x, xs) + sp.diff(F, x).subs(x, xs) * (x - xs)

# Display the linearized form
sp.pprint(F_taylor)

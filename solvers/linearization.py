import numpy as np
import matplotlib.pyplot as plt
import control as ctrl 
import sympy as sp

#### linearization of a system
def  linearization_of_system(eq):
    x = sp.Symbol('x')
    xs = sp.Symbol('x_s') 

    F_taylor = eq.subs(x, xs) + sp.diff(eq, x).subs(x, xs) * (x - xs)

    return sp.pprint(F_taylor)

#### how to use
# x = sp.Symbol('x')
# xs = sp.Symbol('x_s')  # Operating point
# eq = 4 * sp.sqrt(x)
# linearization_of_system(eq)
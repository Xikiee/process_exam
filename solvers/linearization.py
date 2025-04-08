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

#### linearization of a system with multiple variables

def linearization_system_with_multiple_variables(eq, X_o, Y_o):
    X, Y = sp.symbols('X Y')
    eq = eq.subs({X: X_o, Y: Y_o}) + sp.diff(eq,X).subs({X: X_o, Y: Y_o}) * (X - X_o) + sp.diff(eq,Y).subs({X: X_o, Y: Y_o}) * (Y -Y_o)
    return sp.pprint(eq)


# #### how to use 
# X = sp.symbols('X')
# Y = sp.symbols('Y')
# X_o, Y_o = sp.symbols('X_o Y_o') # can use this or act values
# # X_o = 50
# # Y_o = 1.23
# eq = (1.23*X)/(Y+2)
# linearization_system_with_multiple_variables(eq,X_o,Y_o)

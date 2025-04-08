import numpy as np
import matplotlib.pyplot as plt
import control as ctrl 
import sympy as sp


######## partial fraction decomposition
def partial_fraction_decomposition(num, den):
    s = sp.symbols('s', real = True)
    eq = num / den
    partial_fraction_decomposition = sp.simplify(sp.apart(eq))
    return print(partial_fraction_decomposition)

##### inverse laplace transform
def inverse_laplace_transform(eq,s,t):
    output = sp.inverse_laplace_transform(eq,s,t)
    return print(output)


#### splitting the system into imaginary and real parts 
def split_system(num,den):
    s = sp.symbols('s', real = True)
    eq = num / den 
    real_part = sp.simplify(sp.re(eq))
    imag_part = sp.simplify(sp.im(eq))
    return print(f'Real part: {real_part}, Imaginary part: {imag_part}')



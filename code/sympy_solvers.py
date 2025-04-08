import numpy as np
import matplotlib.pyplot as plt
import control as ctrl 
import sympy as sp


######## partial fraction decomposition
def partial_fraction_decomposition(num, den):
    s = sp.symbols('s', real = True)
    eq = num / den
    partial_fraction_decomposition = sp.simplify(sp.apart(eq))
    return print(f"The decomposition is: {partial_fraction_decomposition}")

##### inverse laplace transform
def inverse_laplace_transform(eq,s):
    t = sp.symbols('t',real= True)
    output = sp.inverse_laplace_transform(eq,s,t)
    return print(output)


#### splitting the system into imaginary and real parts 
def split_system(num,den):
    s = sp.symbols('s', real = True)
    eq = num / den 
    real_part = sp.simplify(sp.re(eq))
    imag_part = sp.simplify(sp.im(eq))
    return print(f'Real part: {real_part}, Imaginary part: {imag_part}')


#### finding magnitude 
def find_magnitude(eq):
    s = sp.symbols('s', real = True)
    real_part = sp.simplify(sp.re(eq))
    imag_part = sp.simplify(sp.im(eq))
    magnitude = sp.sqrt(real_part**2 + imag_part**2)
    return print(f'Magnitude: {magnitude}')


#### finding the poles/ zeros of the system
def find_poles_and_zeros(num,den):
    s =sp.symbols('s', real = True)
    num = sp.poly(num,s)
    den = sp.poly(den,s)
    zeros = num.all_roots()
    poles = den.all_roots()
    return print(f'Zeros: {zeros}, Poles: {poles}')


#### sample runs
s =sp.symbols('s', real = True)
num = s+1
den = (s+2)*(s+3)
# partial_fraction_decomposition(num,den)
inverse_laplace_transform(num/den,s)
# find_poles_and_zeros(num,den)


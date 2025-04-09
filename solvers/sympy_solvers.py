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
    w = sp.symbols('w', real = True)
    new_eq = eq.subs(s,w*sp.I)

    real_part = sp.simplify(sp.re(new_eq))
    imag_part = sp.simplify(sp.im(new_eq))
    magnitude = sp.simplify(sp.sqrt(real_part**2 + imag_part**2))
    return print(f'Magnitude: {magnitude}')


#### finding the poles/ zeros of the system
def find_poles_and_zeros(num,den):
    s =sp.symbols('s', real = True)
    num = sp.poly(num,s)
    den = sp.poly(den,s)
    zeros = num.all_roots()
    poles = den.all_roots()
    return print(f'Zeros: {zeros}, Poles: {poles}')


def find_phase_margin(num, den):
    s = sp.symbols('s', real=True)
    eq = num / den

    w = sp.symbols('w', real=True, positive=True)
    eq = eq.subs(s, w * sp.I)

    real_part = sp.re(eq)
    imag_part = sp.im(eq)
    magnitude = sp.sqrt(real_part**2 + imag_part**2)


    w_g_solutions = sp.solve(sp.Eq(magnitude, 1), w)
    w_g_solutions = [sol.evalf() for sol in w_g_solutions if sol.is_real and sol > 0]

    if not w_g_solutions:
        return print("No valid gain crossover frequency found.")

    w_g = w_g_solutions[-1]  

    phase_rad = sp.arg(eq.subs(w, w_g))  
    phase_deg = sp.deg(phase_rad)        

    phase_margin = 180 + phase_deg.evalf()

    return print(f"Cross-over frequency: {w_g:.3f} rad/s, Phase margin: {phase_margin:.2f} degrees")


s = sp.symbols('s', real = True)
num = 40
den = (s+1)**2
find_phase_margin(num,den)
# find_magnitude(num/den)


#### sample runs
# s =sp.symbols('s', real = True)
# num = s+1
# den = (s+2)*(s+3)
# partial_fraction_decomposition(num,den)
# inverse_laplace_transform(num/den,s)
# find_poles_and_zeros(num,den)





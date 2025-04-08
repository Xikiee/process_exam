import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import control as ctrl 

#partial fraction decomposition

# s = sp.symbols('s', real = True)

# num = 3*s + 5
# den = s*(s+6)*(s+4)
# eq = num/den
# partial_fractions = sp.apart(eq)
# print(f"decomposition = {partial_fractions}")



# #find the inverse laplace transform
# s = sp.symbols('s')
# t= sp.symbols('t')
# eq1 = (2*(s+2))/(s*(s+6)*(s+4))
# eq2 = (s+2)/(s*(s+6)*(s+4))
# solve = sp.inverse_laplace_transform((eq1+eq2),s,t)
# print(f"solve = {solve}")

#####################   ploting bode plots
# num_1 = [1,1]
# den_1 = [1,10,24]
# sys1 = ctrl.TransferFunction(num_1,den_1)

# num_2 = [1,2]
# den_2 = [1,10,24]
# sys2 = ctrl.TransferFunction(num_2,den_2)

# sys = sys1 + sys2
# omega = np.logspace(-1, 2, 500)  # From 0.1 to 100 rad/s
# # Plot Bode plot
# ctrl.bode(sys, omega, dB=True)  # `dB=True` for magnitude in decibels
# plt.show()



##### seperating the equation into real and imaginary parts 
# w = sp.symbols('s', real = True)
# eq = (w*sp.I+2)/((w*sp.I)**2 + 10*(w*sp.I) + 24)
# real = sp.re(eq)
# imag = sp.im(eq)

# print(f" real  ={real}")
# print(f"imaginary = {imag}")




#################### find the inverse laplace transform
# s= sp.symbols('s')
# t = sp.symbols('t')

# eq = 40/(s**3*(s+4)) - (40*sp.exp(-10*s))/(s**3*(s+4))
# solve = sp.inverse_laplace_transform(eq,s,t)
# print(f'solve = {solve}')


################ finding the partial fraction decomposition 
s= sp.symbols('s')

eq = (40 * sp.exp(-10*s))/(s*(s**2 + 2*s))
partial_fractions = sp.apart(eq)
print(f'decomposition = {partial_fractions}')
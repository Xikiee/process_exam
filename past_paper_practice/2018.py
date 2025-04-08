import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import control as ctrl

######## partial fraction decomposition

# s, m= sp.symbols('s m')
# eq = ((-s+1))/((s+1)*(0.2*s + 1)*s)
# part_frac = sp.apart(eq)
# print(f'partial fraction = {part_frac}')



#############  plotting bode plots 

# num = [-1,1]
# den = [0.2, 1.2, 1]
# sys = ctrl.TransferFunction(num,den)
# gm, pm,wg,wp = ctrl.margin(sys)             # to find the different margins
# omega = np.logspace(-1,2,500, base=10)
# ctrl.bode(sys, omega, dB=True)
# plt.show()
# print(f"gain margin ={gm}")


# ### finding the imaginary and real parts of the equation

# w = sp.symbols('w',real = True)
# e1  = (-w*sp.I + 1)/((w*sp.I + 1)* (0.2*w*sp.I + 1))
# eq = sp.simplify(e1)
# imag = sp.simplify(sp.im(eq))
# real = sp.re(eq) 
# print(f'real = {real}')
# print(f'imaginary = {imag}')

### finding the step response of the system 
# s= sp.symbols('s')
# h_1 = (-s +1)/((s+1)*(0.2*s+1)*s)
# h_2 = (1/(2*s+1)*s)
# h = h_1 + h_2 
# ## findint the transfer function of the system 
# t = sp.symbols('t')
# output = sp.inverse_laplace_transform(h,s,t)
# print(f'output= {output}')


## finding the poles and zeros of the system
# num_1 = [-1,1]
# den_1 = [0.2,1.2,1]

# num_2 = [1]
# den_2 = [2,1]

# h_1 = ctrl.TransferFunction(num_1, den_1)
# h_2 = ctrl.TransferFunction(num_2,den_2)
# sys = h_1 + h_2
# ctrl.pzmap(sys, Plot = True)
# plt.show()

### partial fraction decomposition 
# s = sp.symbols('s',real = True)
# eq = (2.25 *0.462) /((0.25*s**2 + 1.25*s + 1)*s)
# part_frac = sp.apart(eq)
# print(f'partial fraction = {part_frac}')

#### finding the magnitude
s = sp.symbols('s')
eq = ((s+3)*sp.exp(-s))/(s**2 + 2*s + 2)
imag = sp.simplify(sp.im(eq))
real = sp.simplify(sp.re(eq))

mag = sp.sqrt(real**2 + imag**2)
print(f'magnitude = {mag}')


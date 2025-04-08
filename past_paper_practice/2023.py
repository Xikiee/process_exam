import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import control as ctrl 

# s = sp.symbols('s')
# eq = (-s +1)/((2*s + 1)*(0.5*s + 1)*s)
# part_frac = sp.apart(eq)
# print(part_frac)

## plotting the response
# t = np.linspace(0,50,10000)
# h1 = 10*(-2*np.exp(-t/2) + np.exp(-t/2) + 1)
# h2 = 1- np.exp(-t)
# plt.plot(t,h1, label = 'h1')
# plt.plot(t,h2,label = 'h2')
# plt.show()

### making bode plot h1 
# num = [-10,10]
# den = [1, 2.5, 1]
# sys = ctrl.TransferFunction(num,den)
# gm, pm,wg,wp = ctrl.margin(sys)             # to find the different margins
# omega = np.logspace(-1,2,500, base=10)
# ctrl.bode(sys, omega, dB=True)
# plt.show()
# print(f"gain margin ={gm}")


###3 making bode plots for h2 
# num = [1]
# den = [1,1]
# sys = ctrl.TransferFunction(num,den)
# gm, pm,wg,wp = ctrl.margin(sys)             # to find the different margins
# omega = np.logspace(-1,2,500, base=10)
# ctrl.bode(sys, omega, dB=True)
# plt.show()
# print(f"gain margin ={gm}")

###### plotting the step response of the whole system 
# num1 = [-10,10]
# den1 = [1,2.5,1]
# sys1 = ctrl.TransferFunction(num1,den1)

# num2 =[1]
# den2 = [1,1]
# sys2= ctrl.TransferFunction(num2,den2)

# sys = sys1 + sys2 
# t = np.linspace(0,50,10000)
# t,y = ctrl.step_response(sys,t)
# # plt.plot(t,y)
# # plt.show()

# ## plotting the bode plot for the same system 
# omega = np.logspace(-1,2,500, base=10)
# ctrl.bode(sys, omega, dB=True)
# plt.show()



### finding the poles and zeros of the system
num1 = [-10,10]
den1 = [1,2.5,1]
sys1 = ctrl.TransferFunction(num1,den1)

num2 =[1]
den2 = [1,1]
sys2= ctrl.TransferFunction(num2,den2)

sys = sys1 + sys2 

ctrl.pzmap(sys, Plot = True)
plt.show()

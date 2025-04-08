import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import control as ctrl 

#### h1 system step response plotting

def multi_system_step_response(num1, den1):
    sys = ctrl.TransferFunction(num1,den1)

    t = np.linspace(0,50,10000)
    t,y = ctrl.step_response(sys,t)

    plt.plot(t,y, label = 'step response')
    plt.xlabel('Time (s)')
    plt.ylabel('y(t)')
    plt.grid(which = 'both', linewidth = 0.5)
    plt.legend()
    plt.show()


#### h1 + h2 system ste presonse plotting
def multi_system_step_response(num1, den1, num2, den2):
    sys1 = ctrl.TransferFunction(num1,den1)
    sys2 = ctrl.TransferFunction(num2,den2)
    sys = sys1 + sys2

    t = np.linspace(0,50,10000)
    t,y = ctrl.step_response(sys,t)

    plt.plot(t,y, label = 'step response')
    plt.xlabel('Time (s)')
    plt.ylabel('y(t)')
    plt.grid(which = 'both', linewidth = 0.5)
    plt.legend()
    plt.show()

# num1 = [-10,10]
# den1 = [1,2.5,1]

# num2 =[1]
# den2 = [1,1]
# multi_system_step_response(num1,den1,num2,den2) ## sample plotting

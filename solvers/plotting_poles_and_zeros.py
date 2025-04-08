import numpy as np
import matplotlib.pyplot as plt
import control as ctrl 
import sympy as sp

### finding the poles and zeros of the system
def plotting_poles_and_zeros_multi_sys(num1, den1):
    sys = ctrl.TransferFunction(num1,den1)
    ctrl.pzmap(sys,Plot = True)
    plt.show()

#### plotting poles and zeros of a multi system
def plotting_poles_and_zeros_multi_sys(num1, den1, num2, den2):
    sys1 = ctrl.TransferFunction(num1,den1)
    sys2 = ctrl.TransferFunction(num2,den2)

    sys = sys1 + sys2 
    ctrl.pzmap(sys,Plot = True)
    plt.show()

# num1 = [-10,10]
# den1 = [1,2.5,1]

# num2 =[1]
# den2 = [1,1]
# plotting_poles_and_zeros_multi_sys(num1,den1,num2,den2) ## sample plotting



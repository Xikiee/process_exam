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

#### plotting the root locus of the system
def update_root_locus(num, den, K=1):
    system = ctrl.TransferFunction(num,den)

    plt.figure(figsize=(6, 6))
    
    # Compute closed-loop poles for given K
    closed_loop = ctrl.feedback(K * system)  
    poles = ctrl.poles(closed_loop)  
    
    # Plot root locus
    ctrl.root_locus(system, grid=True)
    
    # Highlight the poles for the given K
    plt.scatter(np.real(poles), np.imag(poles), color='red', s=100, label=f'Poles at K={K}')
    plt.show()

# Interactive slider for gain K
# interact(update_root_locus, K=(0.01, 300, 0.5))

num = [18,6]
den = [3,22,19]
update_root_locus(num, den)

# num1 = [-10,10]
# den1 = [1,2.5,1]

# num2 =[1]
# den2 = [1,1]
# plotting_poles_and_zeros_multi_sys(num1,den1,num2,den2) ## sample plotting




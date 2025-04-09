import numpy as np
import matplotlib.pyplot as plt
import control as ctrl 
import sympy as sp

#### plotting bode plots
def plotting_bode_plot(num, den):
    sys = ctrl.TransferFunction(num,den)
    gm, pm,wg,wp = ctrl.margin(sys)             # to find the different margins
    omega = np.logspace(-1,2,500, base=10)
    ctrl.bode(sys, omega, dB=True)
    plt.xlim(0.1,100)
    plt.show()
    return print(f"gain margin ={gm}, phase margin = {pm}, gain frequency = {wg}, phase frequency = {wp}")

num = [40]
den = [1,2,1]
plotting_bode_plot(num,den) ####### sample plotting


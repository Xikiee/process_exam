import numpy as np
import matplotlib.pyplot as plt
import control as ctrl 
import sympy as sp

#### plotting bode plots
def bode_plot(num, den):
    sys = ctrl.TransferFunction(num,den)
    gm, pm,wg,wp = ctrl.margin(sys)           
    omega = np.logspace(-1,2,500, base=10)
    ctrl.bode(sys, omega, dB=True)
    plt.xlim(0.1,100)
    plt.show()
    return print(f"gain margin ={gm}, phase margin = {pm}, gain frequency = {wg}, phase frequency = {wp}")

# num = [40]
# den = [1,2,1]
# bode_plot(num,den) ####### sample plotting


#### bode plot of a system with time delay
def bode_plot_with_delay(num, den, delay):
    num_d, den_d = ctrl.pade(delay,8)
    d_fun = ctrl.TransferFunction(num_d,den_d)

    sys = ctrl.TransferFunction(num, den)
    sys_d = sys * d_fun
    gm, pm, wg, wp = ctrl.margin(sys_d)
    omega = np.logspace(-1,2,500, base=10)

    ctrl.bode(sys_d, omega, dB=True )
    plt.xlim(0.1,100)
    plt.show()

    return print(f"gain margin = {gm}, phase margin = {pm}, gain frequency = {wg}, phase frequency = {wp}")


# bode_plot_with_delay([1,4],[1,2,4], 1) #### sample plotting


#### plotting bode plots
def bode_plot_multi_sys(num1, den1, num2, den2):
    sys1 = ctrl.TransferFunction(num1,den2)
    sys2 = ctrl.TransferFunction(num2, den2)
    sys = sys1 + sys2
    gm, pm,wg,wp = ctrl.margin(sys)           
    omega = np.logspace(-1,2,500, base=10)
    ctrl.bode(sys, omega, dB=True)
    plt.xlim(0.1,100)
    plt.show()
    return print(f"gain margin ={gm}, phase margin = {pm}, gain frequency = {wg}, phase frequency = {wp}")

# num1 = [40]
# den1 = [1,2,1]
# num2 = [1]
# den2 = [1,2,1]
# bode_plot_multi_sys(num1,den1,num2, den2) ####### sample plotting


#### bode plot of a system with time delay
def bode_plot_with_delay_multi_sys(num1, den1, num2, den2, delay):
    num_d, den_d = ctrl.pade(delay,8)
    d_fun = ctrl.TransferFunction(num_d,den_d)

    sys1 = ctrl.TransferFunction(num1, den1)
    sys2 = ctrl.TransferFunction(num2, den2)
    
    sys_d = (sys1 + sys2) * d_fun
    gm, pm, wg, wp = ctrl.margin(sys_d)
    omega = np.logspace(-1,2,500, base=10)

    ctrl.bode(sys_d, omega, dB=True )
    plt.xlim(0.1,100)
    plt.show()

    return print(f"gain margin = {gm}, phase margin = {pm}, gain frequency = {wg}, phase frequency = {wp}")


# num1 = [40]
# den1 = [1,2,1]
# num2 = [1]
# den2 = [1,2,1]
# delay = 1
# bode_plot_with_delay_multi_sys(num1,den1,num2, den2,delay) ####### sample plotting


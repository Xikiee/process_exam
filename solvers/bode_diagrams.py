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


#### bode plot of a system with time delay
delay = 1
num, den = ctrl.pade(delay, 8)
d_fun = ctrl.TransferFunction(num, den)
# Define transfer function
h1 = ctrl.TransferFunction([1, 4], [1, 2, 4])
h2 = h1*d_fun
# Compute gain margin, phase margin, and crossover frequencies
gm, pm, wg, wp = ctrl.margin(h2)


gm_db = 20 * np.log10(gm)
print(f"Gain Margin (GM): {gm:.5f}")
print(f"Gain Margin db = {gm_db:.5f}")
print(f"Phase Margin (PM): {pm:.5f} degrees")
print(f"Gain Crossover Frequency (Wg): {wg:.5f} rad/s")
print(f"Phase Crossover Frequency (Wp): {wp:.5f} rad/s")

# Plot Bode plot
mag, phase, omega = ctrl.bode_plot(h2, dB=True, Hz=True, margins=True)

# Show plot
plt.show()
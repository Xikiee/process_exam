import numpy as np
import matplotlib.pyplot as plt
import control as ctrl 
import sympy as sp

#### linearization of a system
def  linearization_of_system(eq):
    x = sp.Symbol('x')
    xs = sp.Symbol('x_s') 

    F_taylor = eq.subs(x, xs) + sp.diff(eq, x).subs(x, xs) * (x - xs)

    return sp.pprint(F_taylor)

### how to use
x = sp.Symbol('x')
xs = sp.Symbol('x_s')  # Operating point
eq = (0.05*x)/(2+3*x)
# linearization_of_system(eq)

#### plotting bode plots
def bode_plot(num, den):
    sys = ctrl.TransferFunction(num,den)
    gm, pm,wp,wg = ctrl.margin(sys)           
    omega = np.logspace(-2,3,500, base=10)
    ctrl.bode(sys, omega, dB=True)
    plt.xlim(0.01,500)
    plt.show()
    return print(f"gain margin ={gm}, phase margin = {pm}, phase frequency = {wg}, gain frequency = {wp}")


# bode_plot([1,3],[1,6,5])
s= sp.symbols('s')
eq = (s+1)*(s+50)
# print(sp.expand(eq))


num = [1,3]
den = [1,6,5,0]
# bode_plot(num,den) ####### sample plotting


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


num = s+3
den = s**3+6*s**2 + 5*s
# find_phase_margin(num,den)




def bode_plot(num, den):
    sys = ctrl.TransferFunction(num,den)
    gm, pm,wp,wg = ctrl.margin(sys)           
    omega = np.logspace(-2,3,500, base=10)
    ctrl.bode(sys, omega, dB=True)
    plt.xlim(0.01,500)
    plt.show()
    return print(f"gain margin ={gm}, phase margin = {pm}, phase frequency = {wg}, gain frequency = {wp}")


den = (s*(s+1)*(s+5))
# print(sp.expand(den))

# bode_plot([1,3],[1,6,5,0])


##### inverse laplace transform
def inverse_laplace_transform(eq,s):
    t = sp.symbols('t',real= True)
    output = sp.inverse_laplace_transform(eq,s,t)
    return print(output)


k = sp.symbols('k', real = True)

# num = (k*(s+3))
# den = (s*(s**3+6*s**2+5*s+k*(s+3)))
# inverse_laplace_transform(eq,s)

######## partial fraction decomposition
def partial_fraction_decomposition(num, den):
    s = sp.symbols('s', real = True)
    eq = num / den
    partial_fraction_decomposition = sp.apart(eq)
    return print(f"The decomposition is: {partial_fraction_decomposition}")


# num = s+3
# den = s*(s**3 + 6*s**2 + 6*s + 3)

# -(s + 1)*(s + 5)/(s**3 + 6*s**2 + 6*s + 3)
# num = -(s + 1)*(s + 5)
# den = (s**3 + 6*s**2 + 6*s + 3)



# partial_fraction_decomposition(num,den)

eq= (s+3)
# de = (s*(s+1)*(s+5))
# partial_fraction_decomposition(eq, de)


def system_step_response(num1, den1):
    sys = ctrl.TransferFunction(num1,den1)

    t = np.linspace(0,50,10000)
    t,y = ctrl.step_response(sys,t)

    plt.plot(t,y, label = 'step response')
    plt.xlabel('Time (s)')
    plt.ylabel('y(t)')
    plt.grid(which = 'both', linewidth = 0.5)
    plt.legend()
    plt.show()

num = [0.7, 0.55*3]
den = [1, 6, 5.55, 3*0.7]
system_step_response(num,den)
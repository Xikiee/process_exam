import numpy as np
import matplotlib.pyplot as plt
import control as ctrl 
import sympy as sp

#### feedback system
def feedback(c_s, g_s):
    """
    This function calculates the feedback system of a given transfer function. 
    The feedback system is defined as G(s) / (1 + G(S)H(s)), where G(s) is the transfer function and H(s) is the feedback system. 
    """

    c_s = ctrl.TransferFunction(c_s[0], c_s[1])
    g_s = ctrl.TransferFunction(g_s[0], g_s[1])

    feedback_system = c_s/ (1 + c_s*g_s)
    return print(f"feedback system :{feedback_system}") 


#### how to use 
c_s = [[40],[1,2,1]]
g_s = [[1],[1,1]]
feedback(c_s,g_s)

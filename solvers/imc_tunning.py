import numpy as np
import control as ct
import pandas as pd
import sympy as sp
from sympy.abc import s
from imc_tuning_table import *

def match_to_table(G, controller_type='PI', epsilon=None):
    """
    Match a transfer function to the IMC tuning table and return the tuned controller.
    
    Args:
        G (str or TransferFunction): Transfer function (e.g., "2/(3*s + 1)" or ct.TransferFunction).
        controller_type (str): 'PI' or 'PID'.
        epsilon (float): IMC filter time constant. If None, defaults to τ/3.
    
    Returns:
        dict: Matched model, extracted parameters (k, τ, etc.), and controller TF.
    """
    # Convert string input to symbolic expression (if needed)
    if isinstance(G, str):
        G_sym = sp.sympify(G.replace('^', '**'))
    elif isinstance(G, ct.TransferFunction):
        num = G.num[0][0]
        den = G.den[0][0]
        G_sym = sp.Poly(num, s) / sp.Poly(den, s)
    else:
        raise ValueError("Input must be a transfer function (string or control.TransferFunction).")
    
    # Simplify and factor the transfer function
    G_sym = sp.factor(G_sym)
    k, tau, zeta, beta, epsilon_sym = sp.symbols('k tau zeta beta epsilon', real=True, positive=True)
    
    # Match G_sym to the closest model in the IMC table
    matched_model = None
    params = {}
    
    # Extract numerator and denominator coefficients
    num = sp.Poly(sp.numer(G_sym), s).all_coeffs()
    den = sp.Poly(sp.denom(G_sym), s).all_coeffs()
    
    # Model A: k/(τ*s + 1)
    if len(num) == 1 and len(den) == 2 and den[1] != 0:
        matched_model = 'A'
        params['k'] = float(num[0] / den[1])
        params['tau'] = float(den[0] / den[1])
    
    # Model B: k/((τ₁*s + 1)(τ₂*s + 1)) (simplified for τ₁ = τ₂)
    elif len(num) == 1 and len(den) == 3 and den[2] != 0:
        matched_model = 'B'
        params['k'] = float(num[0] / den[2])
        params['tau1'] = params['tau2'] = float(np.sqrt(den[0] / den[2]))
    
    # Model C: k/(τ²*s² + 2ζτ*s + 1)
    elif len(num) == 1 and len(den) == 3 and den[2] != 0:
        matched_model = 'C'
        params['k'] = float(num[0] / den[2])
        params['tau'] = float(np.sqrt(den[0] / den[2]))
        params['zeta'] = float(den[1] / (2 * params['tau'] * den[2]))
    
    else:
        raise ValueError("Transfer function does not match any model in the IMC table.")
    
    # Default epsilon (IMC filter time constant)
    if epsilon is None:
        epsilon = params.get('tau', 1) / 3  # Conservative choice
    
    # Get tuning rules from the IMC table
    df_imc_pid = pd.DataFrame(imc_pid_table)
    model_row = df_imc_pid[df_imc_pid['Model'] == matched_model].iloc[0]
    
    # Substitute parameters into the tuning formulas and ensure numeric evaluation
    if controller_type == 'PI':
        Kc_expr = sp.sympify(model_row['k_c'].replace('ϵ', 'epsilon'))
        tau_I_expr = sp.sympify(model_row['τ_I'].replace('ϵ', 'epsilon'))
        
        # Substitute and evaluate
        Kc = float(Kc_expr.subs({
            'k': params['k'],
            'tau': params.get('tau', 1),
            'epsilon': epsilon
        }).evalf())
        tau_I = float(tau_I_expr.subs({
            'tau': params.get('tau', 1),
            'epsilon': epsilon
        }).evalf())
        tau_D = 0
    
    elif controller_type == 'PID':
        Kc_expr = sp.sympify(model_row['k_c'].replace('ϵ', 'epsilon'))
        tau_I_expr = sp.sympify(model_row['τ_I'].replace('ϵ', 'epsilon'))
        tau_D_expr = sp.sympify(model_row['τ_D'].replace('ϵ', 'epsilon'))
        
        # Substitute and evaluate
        Kc = float(Kc_expr.subs({
            'k': params['k'],
            'tau': params.get('tau', 1),
            'epsilon': epsilon
        }).evalf())
        tau_I = float(tau_I_expr.subs({
            'tau': params.get('tau', 1),
            'epsilon': epsilon
        }).evalf())
        tau_D = float(tau_D_expr.subs({
            'tau': params.get('tau', 1),
            'epsilon': epsilon
        }).evalf())
    
    else:
        raise ValueError("Controller type must be 'PI' or 'PID'.")
    
    # Create controller transfer function
    if controller_type == 'PI':
        C = ct.TransferFunction([Kc * tau_I, Kc], [tau_I, 0])
    elif controller_type == 'PID':
        C = ct.TransferFunction([Kc * tau_D * tau_I, Kc * tau_I, Kc], [tau_I, 0])
    
    return {
        'Matched Model': matched_model,
        'Parameters': params,
        'Kc': Kc,
        'tau_I': tau_I,
        'tau_D': tau_D,
        'epsilon': epsilon,
        'Controller': C
    }

# Example Usage
if __name__ == "__main__":
    # Input: Transfer function and controller type
    G = "2/(3*s + 1)"  # or ct.TransferFunction([2], [3, 1])
    controller_type = 'PI'
    
    # Match to table and get controller
    result = match_to_table(G, controller_type)
    print(f"Matched Model: {result['Matched Model']}")
    print(f"Parameters: {result['Parameters']}")
    print(f"Kc: {result['Kc']:.2f}")
    print(f"tau_I: {result['tau_I']:.2f}")
    print(f"Controller TF: {result['Controller']}")
    
    # Simulate closed-loop response
    sys_cl = ct.feedback(result['Controller'] * ct.TransferFunction([2], [3, 1]), 1)
    t, y = ct.step_response(sys_cl)
    import matplotlib.pyplot as plt
    plt.plot(t, y)
    plt.xlabel('Time')
    plt.ylabel('Output')
    plt.title('Closed-Loop Step Response')
    plt.grid()
    plt.show()
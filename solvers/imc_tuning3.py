from sympy import symbols, simplify, Eq
from sympy.abc import s

def generalize_tf(tf_expr):
    """
    Generalize a given transfer function using standard IMC-PID notation.
    Currently works for first- and second-order transfer functions.
    
    Parameters:
    - tf_expr: sympy expression representing the transfer function (function of 's')
    
    Returns:
    - Generalized symbolic expression
    """
    # Define symbols
    k, tau1, tau2, beta, zeta = symbols('k tau1 tau2 beta zeta', real=True, positive=True)
    
    # Try matching different known forms
    # 1st Order: k / (tau1*s + 1)
    tf1 = k / (tau1*s + 1)
    
    # 1st Order + Zero: k*(-beta*s + 1)/(tau1*s + 1)
    tf2 = k * (-beta*s + 1) / (tau1*s + 1)
    
    # 2nd Order: k / (tau1**2 * s**2 + 2*zeta*tau1*s + 1)
    tf3 = k / (tau1**2 * s**2 + 2*zeta*tau1*s + 1)
    
    # 2nd Order + Zero: k*(-beta*s + 1)/(tau1**2 * s**2 + 2*zeta*tau1*s + 1)
    tf4 = k * (-beta*s + 1) / (tau1**2 * s**2 + 2*zeta*tau1*s + 1)
    
    # Try matching
    models = [tf1, tf2, tf3, tf4]
    for model in models:
        if simplify(tf_expr / model).is_constant():
            return simplify(tf_expr / model) * model

    return "No standard form match found."

# Example usage:
from sympy import Rational
numerical_tf = Rational(2) / (3*s + 1)
general_form = generalize_tf(numerical_tf)
print("Generalized form:", general_form)

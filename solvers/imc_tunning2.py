import numpy as np
import control as ct
import pandas as pd
import sympy as sp
from sympy.abc import s
from imc_tuning_table import *
import re

def identify_model_and_calculate_params(eq, epsilon=1.0):
    """
    Identifies the model type and calculates PID parameters based on IMC tuning rules
    
    Args:
        eq: The process model equation as a string (e.g., "2/(3s+1)")
        epsilon: The IMC filter time constant (default=1.0)
        
    Returns:
        A dictionary containing:
        - model_type: The identified model (A-N)
        - parameters: Extracted model parameters
        - pid_params: Dictionary with k_c, τ_I, τ_D, τ_F
        - comments: Any special notes about the model
    """
    eq = eq.replace(" ", "").replace("²", "2").lower()
    
    # Helper function to extract parameters
    def get_params(match, names):
        return {name: float(match.group(i+1)) for i, name in enumerate(names)}
    
    result = {
        'model_type': None,
        'parameters': {},
        'pid_params': {'k_c': None, 'τ_I': None, 'τ_D': None, 'τ_F': None},
        'comments': None
    }
    
    # Model A: k/(τs + 1)
    match = re.fullmatch(r"([\d\.+-]+)/\(([\d\.+-]+)s\+([\d\.+-]+)\)", eq)
    if match and float(match.group(3)) == 1:
        params = get_params(match, ['k', 'τ'])
        result.update({
            'model_type': "A",
            'parameters': params,
            'pid_params': {
                'k_c': params['τ'] / (params['k'] * epsilon),
                'τ_I': params['τ'],
                'τ_D': None,
                'τ_F': None
            },
            'comments': "First-order process"
        })
        return result
    
    # Model B: k/((τ₁s + 1)(τ₂s + 1))
    match = re.fullmatch(r"([\d\.+-]+)/\(\(([\d\.+-]+)s\+1\)\(([\d\.+-]+)s\+1\)\)", eq)
    if match:
        params = get_params(match, ['k', 'τ₁', 'τ₂'])
        result.update({
            'model_type': "B",
            'parameters': params,
            'pid_params': {
                'k_c': (params['τ₁'] + params['τ₂']) / (params['k'] * epsilon),
                'τ_I': params['τ₁'] + params['τ₂'],
                'τ_D': (params['τ₁'] * params['τ₂']) / (params['τ₁'] + params['τ₂']),
                'τ_F': None
            },
            'comments': "Second-order process (two time constants)"
        })
        return result
    
    # Model C: k/(τ2s2 + 2ζτs + 1)
    match = re.fullmatch(r"([\d\.+-]+)/\(([\d\.+-]+)s2\+([\d\.+-]+)s\+1\)", eq)
    if match:
        params = get_params(match, ['k', 'τ2', '2ζτ'])
        params['τ'] = (params['τ2'])**0.5
        params['ζ'] = params['2ζτ'] / (2 * params['τ'])
        result.update({
            'model_type': "C",
            'parameters': params,
            'pid_params': {
                'k_c': params['2ζτ'] / (params['k'] * epsilon),
                'τ_I': 2 * params['ζ'] * params['τ'],
                'τ_D': params['τ'] / (2 * params['ζ']),
                'τ_F': None
            },
            'comments': "Second-order underdamped process"
        })
        return result
    
    # Model D: k(-βs + 1)/(τs + 1)
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(([\d\.+-]+)s\+([\d\.+-]+)\)", eq)
    if match and float(match.group(4)) == 1:
        params = get_params(match, ['k', 'β', 'τ'])
        result.update({
            'model_type': "D",
            'parameters': params,
            'pid_params': {
                'k_c': params['τ'] / (params['k'] * (epsilon + params['β'])),
                'τ_I': params['τ'],
                'τ_D': None,
                'τ_F': epsilon
            },
            'comments': "First-order with RHP zero (β > 0)"
        })
        return result
    
    # Model E: k(-βs + 1)/((τs + 1)(βs + 1))
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(\(([\d\.+-]+)s\+1\)\(([\d\.+-]+)s\+1\)\)", eq)
    if match and float(match.group(3)) == float(match.group(4)):
        params = get_params(match, ['k', 'β', 'τ'])
        result.update({
            'model_type': "E",
            'parameters': params,
            'pid_params': {
                'k_c': params['τ'] / (params['k'] * (2*epsilon + params['β'])),
                'τ_I': params['τ'],
                'τ_D': None,
                'τ_F': 2*epsilon + params['β']
            },
            'comments': "First-order with RHP zero and matching pole (β > 0)"
        })
        return result
    
    # Model F: k(-βs + 1)/(τ2s2 + 2ζτs + 1)
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(([\d\.+-]+)s2\+([\d\.+-]+)s\+1\)", eq)
    if match:
        params = get_params(match, ['k', 'β', 'τ2', '2ζτ'])
        params['τ'] = (params['τ2'])**0.5
        params['ζ'] = params['2ζτ'] / (2 * params['τ'])
        result.update({
            'model_type': "F",
            'parameters': params,
            'pid_params': {
                'k_c': params['2ζτ'] / (params['k'] * (epsilon + params['β'])),
                'τ_I': 2 * params['ζ'] * params['τ'],
                'τ_D': params['τ'] / (2 * params['ζ']),
                'τ_F': epsilon
            },
            'comments': "Second-order underdamped with RHP zero (β > 0)"
        })
        return result
    
    # Model G: k(-βs + 1)/((τ2s2 + 2ζτs + 1)(βs + 1))
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(\(([\d\.+-]+)s2\+([\d\.+-]+)s\+1\)\(([\d\.+-]+)s\+1\)\)", eq)
    if match and float(match.group(5)) == float(match.group(2)):
        params = get_params(match, ['k', 'β', 'τ2', '2ζτ'])
        params['τ'] = (params['τ2'])**0.5
        params['ζ'] = params['2ζτ'] / (2 * params['τ'])
        result.update({
            'model_type': "G",
            'parameters': params,
            'pid_params': {
                'k_c': params['2ζτ'] / (params['k'] * (2*epsilon + params['β'])),
                'τ_I': 2 * params['ζ'] * params['τ'],
                'τ_D': params['τ'] / (2 * params['ζ']),
                'τ_F': 2*epsilon + params['β']
            },
            'comments': "Second-order underdamped with RHP zero and matching pole (β > 0)"
        })
        return result
    
    # Model H: k/s
    match = re.fullmatch(r"([\d\.+-]+)/s", eq)
    if match:
        params = get_params(match, ['k'])
        result.update({
            'model_type': "H",
            'parameters': params,
            'pid_params': {
                'k_c': 1 / (params['k'] * epsilon),
                'τ_I': None,
                'τ_D': None,
                'τ_F': None
            },
            'comments': "Integrator process"
        })
        return result
    
    # Model I: k(2s + 1)/s
    match = re.fullmatch(r"([\d\.+-]+)\(([\d\.+-]+)s\+1\)/s", eq)
    if match and float(match.group(2)) == 2:
        params = get_params(match, ['k'])
        result.update({
            'model_type': "I",
            'parameters': params,
            'pid_params': {
                'k_c': 2 / (params['k'] * epsilon),
                'τ_I': 2,
                'τ_D': 2*epsilon,
                'τ_F': None
            },
            'comments': "Integrator with lead"
        })
        return result
    
    # Model J: k/(s(τs + 1))
    match = re.fullmatch(r"([\d\.+-]+)/\(s\(([\d\.+-]+)s\+([\d\.+-]+)\)\)", eq)
    if match and float(match.group(3)) == 1:
        params = get_params(match, ['k', 'τ'])
        result.update({
            'model_type': "J",
            'parameters': params,
            'pid_params': {
                'k_c': 1 / (params['k'] * epsilon),
                'τ_I': None,
                'τ_D': None,
                'τ_F': params['τ']
            },
            'comments': "Integrator with first-order lag"
        })
        return result
    
    # Model K: k(2ϵs + 1)/(s(τs + 1)(2ϵs + 1))
    match = re.fullmatch(r"([\d\.+-]+)\(([\d\.+-]+)s\+1\)/\(s\(([\d\.+-]+)s\+1\)\(([\d\.+-]+)s\+1\)\)", eq)
    if match and float(match.group(2)) == 2*epsilon and float(match.group(4)) == 2*epsilon:
        params = get_params(match, ['k', 'τ'])
        result.update({
            'model_type': "K",
            'parameters': params,
            'pid_params': {
                'k_c': (2*epsilon + params['τ']) / (params['k'] * epsilon**2),
                'τ_I': 2*epsilon + params['τ'],
                'τ_D': (2*epsilon * params['τ']) / (2*epsilon + params['τ']),
                'τ_F': None
            },
            'comments': "Complex integrator process"
        })
        return result
    
    # Model L: k(-βs + 1)/s
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/s", eq)
    if match:
        params = get_params(match, ['k', 'β'])
        result.update({
            'model_type': "L",
            'parameters': params,
            'pid_params': {
                'k_c': 1 / (params['k'] * (epsilon + params['β'])),
                'τ_I': None,
                'τ_D': None,
                'τ_F': epsilon
            },
            'comments': "Integrator with RHP zero (β > 0)"
        })
        return result
    
    # Model M: k(-βs + 1)/(ϵs + 1)
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(([\d\.+-]+)s\+([\d\.+-]+)\)", eq)
    if match and float(match.group(3)) == epsilon and float(match.group(4)) == 1:
        params = get_params(match, ['k', 'β'])
        result.update({
            'model_type': "M",
            'parameters': params,
            'pid_params': {
                'k_c': 1 / (params['k'] * (2*epsilon + params['β'])),
                'τ_I': 1,
                'τ_D': None,
                'τ_F': 2*epsilon + params['β']
            },
            'comments': "First-order with RHP zero and ϵ filter (β > 0)"
        })
        return result
    
    # Model N: k(-βs + 1)/((ϵs + 1)(2s + 1))
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(\(([\d\.+-]+)s\+1\)\(([\d\.+-]+)s\+1\)\)", eq)
    if match and float(match.group(3)) == epsilon and float(match.group(4)) == 2:
        params = get_params(match, ['k', 'β'])
        result.update({
            'model_type': "N",
            'parameters': params,
            'pid_params': {
                'k_c': 1 / (params['k'] * (2*epsilon + params['β'])),
                'τ_I': 2*epsilon,
                'τ_D': 2*epsilon,
                'τ_F': 2*epsilon + params['β']
            },
            'comments': "Complex dynamics with RHP zero (β > 0)"
        })
        return result
    
    return result

#### how to use
################ input (s+1) as (1s+1)
################ input (s^2+1) as (1s2+1)

# eq_1 = sp.simplify*()

eq = "5/((1s+1)(2s+1))" #### insert your equation here as a string 

result = identify_model_and_calculate_params(eq, epsilon=1.0)
model = result['model_type']
if model:
    print(f"Model {model} with parameters: {result['parameters']}")
    print(f"  k_c: {result['pid_params']['k_c']:.4f}")
    print(f"  τ_I: {result['pid_params']['τ_I'] if result['pid_params']['τ_I'] is not None else '-'}")
    print(f"  τ_D: {result['pid_params']['τ_D'] if result['pid_params']['τ_D'] is not None else '-'}")
    print(f"  τ_F: {result['pid_params']['τ_F'] if result['pid_params']['τ_F'] is not None else '-'}")
else:
    print("No matching model found")
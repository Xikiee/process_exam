import numpy as np
import control as ct
import pandas as pd
import sympy as sp
from sympy.abc import s
from imc_tuning_table import *
import re

def identify_model_and_calculate_kc(eq, epsilon=1.0):
    """
    Identifies the model type and calculates k_c based on the IMC tuning rules
    
    Args:
        eq: The process model equation as a string (e.g., "2/(3s+1)")
        epsilon: The IMC filter time constant (default=1.0)
        
    Returns:
        A tuple containing (model_type, parameters, k_c_value)
        or (None, {}, None) if no match found
    """
    eq = eq.replace(" ", "").replace("²", "2").lower()
    
    # Helper function to extract parameters
    def get_params(match, names):
        return {name: float(match.group(i+1)) for i, name in enumerate(names)}
    
    # Model A: k/(τs + 1)
    match = re.fullmatch(r"([\d\.+-]+)/\(([\d\.+-]+)s\+([\d\.+-]+)\)", eq)
    if match and float(match.group(3)) == 1:
        params = get_params(match, ['k', 'τ'])
        k_c = params['τ'] / (params['k'] * epsilon)
        return "A", params, k_c
    
    # Model B: k/((τ₁s + 1)(τ₂s + 1))
    match = re.fullmatch(r"([\d\.+-]+)/\(\(([\d\.+-]+)s\+1\)\(([\d\.+-]+)s\+1\)\)", eq)
    if match:
        params = get_params(match, ['k', 'τ₁', 'τ₂'])
        k_c = (params['τ₁'] + params['τ₂']) / (params['k'] * epsilon)
        return "B", params, k_c
    
    # Model C: k/(τ2s2 + 2ζτs + 1)
    match = re.fullmatch(r"([\d\.+-]+)/\(([\d\.+-]+)s2\+([\d\.+-]+)s\+1\)", eq)
    if match:
        params = get_params(match, ['k', 'τ2', '2ζτ'])
        params['τ'] = (params['τ2'])**0.5
        params['ζ'] = params['2ζτ'] / (2 * params['τ'])
        k_c = params['2ζτ'] / (params['k'] * epsilon)
        return "C", params, k_c
    
    # Model D: k(-βs + 1)/(τs + 1)
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(([\d\.+-]+)s\+([\d\.+-]+)\)", eq)
    if match and float(match.group(4)) == 1:
        params = get_params(match, ['k', 'β', 'τ'])
        k_c = params['τ'] / (params['k'] * (epsilon + params['β']))
        return "D", params, k_c
    
    # Model E: k(-βs + 1)/((τs + 1)(βs + 1))
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(\(([\d\.+-]+)s\+1\)\(([\d\.+-]+)s\+1\)\)", eq)
    if match and float(match.group(3)) == float(match.group(4)):
        params = get_params(match, ['k', 'β', 'τ'])
        k_c = params['τ'] / (params['k'] * (2*epsilon + params['β']))
        return "E", params, k_c
    
    # Model F: k(-βs + 1)/(τ2s2 + 2ζτs + 1)
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(([\d\.+-]+)s2\+([\d\.+-]+)s\+1\)", eq)
    if match:
        params = get_params(match, ['k', 'β', 'τ2', '2ζτ'])
        params['τ'] = (params['τ2'])**0.5
        params['ζ'] = params['2ζτ'] / (2 * params['τ'])
        k_c = params['2ζτ'] / (params['k'] * (epsilon + params['β']))
        return "F", params, k_c
    
    # Model G: k(-βs + 1)/((τ2s2 + 2ζτs + 1)(βs + 1))
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(\(([\d\.+-]+)s2\+([\d\.+-]+)s\+1\)\(([\d\.+-]+)s\+1\)\)", eq)
    if match and float(match.group(5)) == float(match.group(2)):
        params = get_params(match, ['k', 'β', 'τ2', '2ζτ'])
        params['τ'] = (params['τ2'])**0.5
        params['ζ'] = params['2ζτ'] / (2 * params['τ'])
        k_c = params['2ζτ'] / (params['k'] * (2*epsilon + params['β']))
        return "G", params, k_c
    
    # Model H: k/s
    match = re.fullmatch(r"([\d\.+-]+)/s", eq)
    if match:
        params = get_params(match, ['k'])
        k_c = 1 / (params['k'] * epsilon)
        return "H", params, k_c
    
    # Model I: k(2s + 1)/s
    match = re.fullmatch(r"([\d\.+-]+)\(([\d\.+-]+)s\+1\)/s", eq)
    if match and float(match.group(2)) == 2:
        params = get_params(match, ['k'])
        k_c = 2 / (params['k'] * epsilon)
        return "I", params, k_c
    
    # Model J: k/(s(τs + 1))
    match = re.fullmatch(r"([\d\.+-]+)/\(s\(([\d\.+-]+)s\+([\d\.+-]+)\)\)", eq)
    if match and float(match.group(3)) == 1:
        params = get_params(match, ['k', 'τ'])
        k_c = 1 / (params['k'] * epsilon)
        return "J", params, k_c
    
    # Model K: k(2ϵs + 1)/(s(τs + 1)(2ϵs + 1))
    match = re.fullmatch(r"([\d\.+-]+)\(([\d\.+-]+)s\+1\)/\(s\(([\d\.+-]+)s\+1\)\(([\d\.+-]+)s\+1\)\)", eq)
    if match and float(match.group(2)) == 2*epsilon and float(match.group(4)) == 2*epsilon:
        params = get_params(match, ['k', 'τ'])
        k_c = (2*epsilon + params['τ']) / (params['k'] * epsilon**2)
        return "K", params, k_c
    
    # Model L: k(-βs + 1)/s
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/s", eq)
    if match:
        params = get_params(match, ['k', 'β'])
        k_c = 1 / (params['k'] * (epsilon + params['β']))
        return "L", params, k_c
    
    # Model M: k(-βs + 1)/(ϵs + 1)
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(([\d\.+-]+)s\+([\d\.+-]+)\)", eq)
    if match and float(match.group(3)) == epsilon and float(match.group(4)) == 1:
        params = get_params(match, ['k', 'β'])
        k_c = 1 / (params['k'] * (2*epsilon + params['β']))
        return "M", params, k_c
    
    # Model N: k(-βs + 1)/((ϵs + 1)(2s + 1))
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(\(([\d\.+-]+)s\+1\)\(([\d\.+-]+)s\+1\)\)", eq)
    if match and float(match.group(3)) == epsilon and float(match.group(4)) == 2:
        params = get_params(match, ['k', 'β'])
        k_c = 1 / (params['k'] * (2*epsilon + params['β']))
        return "N", params, k_c
    
    return None, {}, None


#### calculation

eq = "2/(3s+1)" #### insert your equation here as a string

model, params, k_c = identify_model_and_calculate_kc(eq, epsilon=1.0) 
if model:
    print(f"Model {model} with parameters: {params}")
    print(f"Calculated k_c value: {k_c}")
else:
    print("No matching model found")
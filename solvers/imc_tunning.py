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
    """
    eq = eq.replace(" ", "").lower()
    
    # Model A: k/(τs+1)
    match = re.fullmatch(r"([\d\.+-]+)/\(([\d\.+-]+)s\+([\d\.+-]+)\)", eq)
    if match and float(match.group(3)) == 1:
        k = float(match.group(1))
        τ = float(match.group(2))
        k_c = τ / (k * epsilon)
        return "A", {"k": k, "τ": τ}, k_c
    
    # Model B: k/((τ₁s+1)(τ₂s+1))
    match = re.fullmatch(r"([\d\.+-]+)/\(\(([\d\.+-]+)s\+1\)\(([\d\.+-]+)s\+1\)\)", eq)
    if match:
        k = float(match.group(1))
        τ1 = float(match.group(2))
        τ2 = float(match.group(3))
        k_c = (τ1 + τ2) / (k * epsilon)
        return "B", {"k": k, "τ₁": τ1, "τ₂": τ2}, k_c
    
    # Model C: k/(τ²s² + 2ζτs + 1)
    match = re.fullmatch(r"([\d\.+-]+)/\(([\d\.+-]+)s²\+([\d\.+-]+)s\+1\)", eq)
    if match:
        k = float(match.group(1))
        τ_sq = float(match.group(2))
        two_zeta_tau = float(match.group(3))
        τ = (τ_sq)**0.5
        ζ = two_zeta_tau / (2 * τ)
        k_c = two_zeta_tau / (k * epsilon)
        return "C", {"k": k, "τ": τ, "ζ": ζ}, k_c
    
    # Model D: k(-βs + 1)/(τs + 1)
    match = re.fullmatch(r"([\d\.+-]+)\(-([\d\.+-]+)s\+1\)/\(([\d\.+-]+)s\+([\d\.+-]+)\)", eq)
    if match:
        k = float(match.group(1))
        β = float(match.group(2))
        τ = float(match.group(3))
        k_c = τ / (k * (epsilon + β))
        return "D", {"k": k, "β": β, "τ": τ}, k_c
    
    # Add patterns for other models as needed...
    
    return None, {}, None

# Test with your equation
model, params, k_c = identify_model_and_calculate_kc("2/(3s+1)", epsilon=1.0)
if model:
    print(f"Model {model} with parameters: {params}")
    print(f"Calculated k_c value: {k_c}")
else:
    print("No matching model found")
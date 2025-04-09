import sympy as sp
import pandas as pd
from imc_tuning_table import imc_pid_table

# Convert the dictionary to a DataFrame
df_imc_pid = pd.DataFrame(imc_pid_table)

def match_to_table(transfer_func, controller_type="PI", epsilon_value=None):
    # Define consistent symbols
    s = sp.Symbol('s')
    k, tau, zeta, epsilon, beta, tau1, tau2 = sp.symbols('k tau zeta epsilon beta tau1 tau2')

    # Simplify and expand transfer function
    num_expr, den_expr = sp.fraction(sp.simplify(transfer_func))
    input_num = list(reversed(sp.Poly(num_expr, s).all_coeffs()))
    input_den = list(reversed(sp.Poly(den_expr, s).all_coeffs()))

    matched_rows = []

    # Locals dictionary for sympify to avoid 'taus' issue
    locals_dict = {
        's': s, 'k': k, 'tau': tau, 'zeta': zeta, 'epsilon': epsilon,
        'beta': beta, 'tau1': tau1, 'tau2': tau2
    }

    for i, row in df_imc_pid.iterrows():
        try:
            model_expr = row["Process Model (g)"]
            model_expr_cleaned = (
                model_expr.replace("ϵ", "epsilon")
                          .replace("ζ", "zeta")
                          .replace("τ", "tau")
                          .replace("β", "beta")
                          .replace("τ₁", "tau1")
                          .replace("τ₂", "tau2")
            )

            # Use proper context
            model_expr_sympy = sp.sympify(model_expr_cleaned, locals=locals_dict)
            model_num, model_den = sp.fraction(sp.simplify(model_expr_sympy))

            model_num_coeffs = list(reversed(sp.Poly(model_num, s).all_coeffs()))
            model_den_coeffs = list(reversed(sp.Poly(model_den, s).all_coeffs()))

            params = list(model_expr_sympy.free_symbols - {s})
            equations = []

            # Match coefficients
            for model_c, input_c in zip(model_num_coeffs, input_num):
                equations.append(sp.Eq(model_c, input_c))
            for model_c, input_c in zip(model_den_coeffs, input_den):
                equations.append(sp.Eq(model_c, input_c))

            solution = sp.solve(equations, params, dict=True)

            if solution:
                final_solution = {}
                for param, val in solution[0].items():
                    final_solution[str(param)] = float(val.evalf()) if val.free_symbols == set() else str(val)
                matched_rows.append((i, final_solution))

        except Exception as e:
            continue

    if not matched_rows:
        return "No match found in the table."

    # Controller formulas
    controller_columns = {
        "P": ["k_c"],
        "PI": ["k_c", "τ_I"],
        "PID": ["k_c", "τ_I", "τ_D"],
        "PIDF": ["k_c", "τ_I", "τ_D", "τ_F"]
    }.get(controller_type.upper(), [])

    results = []
    for row_idx, params in matched_rows:
        row = df_imc_pid.iloc[row_idx]
        result = {
            "Model": row["Model"],
            "Process Model": row["Process Model (g)"],
            "Parameters": params
        }

        # Compute k_c numerically
        if epsilon_value is not None and "k_c" in row:
            try:
                k_c_expr = row["k_c"].replace("ϵ", "epsilon").replace("τ", "tau").replace("k", "k")
                k_c_sym = sp.sympify(k_c_expr, locals=locals_dict)
                subs_dict = {sp.Symbol(k): v for k, v in params.items() if isinstance(v, (int, float))}
                subs_dict[epsilon] = epsilon_value
                result["k_c_numeric"] = float(k_c_sym.subs(subs_dict).evalf())
            except Exception:
                result["k_c_numeric"] = "Could not evaluate"

        for col in controller_columns:
            result[col] = row[col]

        results.append(result)

    return results







#### example use
s = sp.Symbol('s')
G = 2 / (3*s + 1)
print(match_to_table(G, "PI", epsilon_value=1))








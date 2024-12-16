import numpy as np

from ortools.linear_solver import pywraplp
from sklearn.linear_model import LinearRegression


def do_regression(lmps: list[list[float]], Y: np.ndarray):
    """
    Perform linear regression on lineage mutation profiles against observed data,
    adjusting the model iteratively if the sum of coefficients exceeds 1 to find the best model
    based on regression score that satisfies this constraint.

    Args:
        lmps: List of lineage mutation profiles.
        Y: Observed mutation frequencies.

    Returns:
        tuple: Tuple containing the matrix of mutation profiles and the best regression coefficients.
    """
    X = np.array(lmps).T
    reg = LinearRegression(fit_intercept=False, positive=True).fit(X, Y)
    best_score = reg.score(X, Y)
    best_reg = reg

    if sum(reg.coef_) > 1:
        for i in range(len(X)):
            if Y[i] == 0:
                new_X = np.concatenate((X[:i], X[i + 1 :]))
                new_Y = np.concatenate((Y[:i], Y[i + 1 :]))
                new_reg = LinearRegression(fit_intercept=False, positive=True).fit(
                    new_X, new_Y
                )
                new_score = new_reg.score(new_X, new_Y)
                if sum(new_reg.coef_) <= 1 and new_score > best_score:
                    best_reg = new_reg
                    best_score = new_score

    return X, [coef for coef in best_reg.coef_]


def do_regression_linear(lmps, Y, muts):
    """
    Perform a linear programming approach to minimize the error between predicted and observed mutation frequencies.

    Args:
        lmps (list of list of floats): Lineage mutation profiles.
        Y (list of floats): Observed mutation frequencies.
        muts (list of str): List of mutation identifiers.

    Returns:
        tuple: Matrix of mutation profiles, solution values for lineages, and a dictionary of mutation differences.
    """

    solver = pywraplp.Solver.CreateSolver("GLOP")
    num_lins = len(lmps)
    num_muts = len(lmps[0])
    lins = [solver.NumVar(0, 1, f"x_{i}") for i in range(num_lins)]
    t = [solver.NumVar(0, solver.infinity(), f"t_{i}") for i in range(num_muts)]

    # Set up constraints
    for i, mut_freq in enumerate(Y):
        constraint_less = solver.Constraint(-solver.infinity(), mut_freq)
        constraint_more = solver.Constraint(mut_freq, solver.infinity())
        for j in range(num_lins):
            constraint_less.SetCoefficient(lins[j], lmps[j][i])
            constraint_more.SetCoefficient(lins[j], lmps[j][i])
        constraint_less.SetCoefficient(t[i], -1)
        constraint_more.SetCoefficient(t[i], 1)

    # Objective: minimize the sum of t, which represents the residuals
    objective = solver.Objective()
    for ti in t:
        objective.SetCoefficient(ti, 1)
    objective.SetMinimization()

    solver.Solve()

    X = np.array(lmps).T
    mut_diffs = {
        mut: Y[i] - sum(lins[j].solution_value() * lmps[j][i] for j in range(num_lins))
        for i, mut in enumerate(muts)
    }

    return X, [lin.solution_value() for lin in lins], mut_diffs

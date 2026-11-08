# SEIQR COVID-19 model: solve the equations and plot the five compartments
# (this is the code of Appendix C of the project report)
# Case 1 (R0 > 1): beta = 4e-5.   For Case 2 (R0 < 1) change beta to 1e-5.
# Run on your computer (python seiqr_charts.py) or on Google Colab.

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

A, mu, beta = 100, 0.02, 4e-5
alpha, gamma, sigma, theta = 0.1, 0.05, 0.08, 0.07
k1, k2, k3 = alpha + mu + gamma, sigma + mu, theta + mu

def model(t, y):
    S, E, I, Q, R = y
    inf = beta * S * (E + I)
    return [A - mu * S - inf,
            inf - k1 * E,
            alpha * E - k2 * I,
            gamma * E + sigma * I - k3 * Q,
            theta * Q - mu * R]

t = np.linspace(0, 365, 3651)
sol = solve_ivp(model, [0, 365], [4990, 10, 0, 0, 0], t_eval=t)

for y, name in zip(sol.y, ["S", "E", "I", "Q", "R"]):
    plt.figure()
    plt.plot(t, y)
    plt.xlabel("time (days)")
    plt.ylabel(name)
plt.show()

"""Notations:
S: Susceptible
I: Infected
R: Reconvered
Q: Quarantined

alpha: R -> S
beta: S -> I
gamma: I -> R
phi: S -> Q

I_cutoff: threshold after which S -> Q becomes significant
R_hist[i]: recovered cases until i-th day
"""

import numpy as np
import matplotlib.pyplot as plt


def fear_function(I, dIdt, phi):
    """This function decides likely hood of S -> Q

    Parameters
    ----------
    I : array like
        Infected
    dIdt : array like
        rate of change in Infected
    phi : float
        controlling parameter

    Returns
    -------
    array like
    """
    # if I < I_cutoff: 0 else phi * ((I - I_cutoff)) ** 2
    return np.piecewise(
        I, [I < I_cutoff, I >= I_cutoff], [0, lambda I: phi * ((I - I_cutoff)) ** 2]
    )


def dydt(t, y, R_hist, params):
    """Create the system of eqaution

    Parameters
    ----------
    t : float
        time
    y : array like
        [S, I, R, Q]
    R_hist : array like
        R in the the previous days
    params : array like
        list of parameters, [alpha, beta, gamma, phi]

    Returns
    -------
    array like
        list of differential equation
    """
    [S, I, R, Q] = y
    [alpha, beta, gamma, phi] = params

    dIdt = beta * I * S / N - gamma * I
    dQdt = S * fear_function(I, dIdt, phi) / N

    # waning is considered only after 6 months
    if len(R_hist) > 180:
        dRdt = gamma * I - alpha * (R_hist[-179] - R_hist[-180])
        dSdt = (
            -beta * I * S / N
            - S * fear_function(I, dIdt, phi) / N
            + alpha * (R_hist[-179] - R_hist[-180])
        )

    else:
        dRdt = gamma * I
        dSdt = -beta * I * S / N - S * fear_function(I, dIdt, phi) / N

    return np.array([dSdt, dIdt, dRdt, dQdt])


def custom_rk4(t, y0, t0, I_cutoff, params):
    """Custom rk4 integrator

    Parameters
    ----------
    t : float
        time
    y0 : array like
        initial values of [S, I, R, Q]
    t0 : float
        initial time
    I_cutoff : int
        threshold after which S -> Q becomes significant
    params : array like
        list of parameters, [alpha, beta, gamma, phi]

    Returns
    -------
    list
        list of [S, I, R, Q] where each of which are themselves np.array
    """
    h = 1e-2
    n = int((t - t0) / h)
    y = y0
    R_hist = []
    y_hist = []
    for i in range(1, n + 1):

        [S, I, R, Q] = y

        # store only every 100th step
        if i % (int(1 / h)) == 1:
            R_hist.append(R)
            y_hist.append(y)

        if I < I_cutoff:
            y = np.array([S + Q, I, R, 0])

        k1 = h * dydt(t0, y, R_hist, params)
        k2 = h * dydt(t0 + 0.5 * h, y + 0.5 * k1, R_hist, params)
        k3 = h * dydt(t0 + 0.5 * h, y + 0.5 * k2, R_hist, params)
        k4 = h * dydt(t0 + h, y + k3, R_hist, params)
        y = y + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

        t0 = t0 + h
        y = [y[i] for i in range(4)]

    return np.transpose((np.array(y_hist)))


if __name__ == "__main__":

    alpha = 0.5
    beta = 0.08
    gamma = 0.02
    phi = 1e-3
    N = 1.3526e9
    I_cutoff = 1e4

    S_init, I_init, R_init, Q_init = N - 2006, 2000, 6, 0
    y_init = [S_init, I_init, R_init, Q_init]
    params = [alpha, beta, gamma, phi]

    days = 400

    [S, I, R, Q] = custom_rk4(days, y_init, 0, I_cutoff, params)

    plt.figure()
    tot = S + I + R + Q
    plt.plot(np.arange(days), I, "b", label="Infectives")
    # plt.plot(np.arange(days), tot, label="total")
    # plt.plot(np.arange(days), R, "g", label="Recovered")
    # plt.plot(np.arange(days), S, "r", label="Susceptible")
    # plt.plot(np.arange(days), Q, "black", label="Quarantined")
    plt.legend(loc="best")
    plt.xlabel("time (days)")
    plt.grid()
    plt.show()

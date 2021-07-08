import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from model import Model, integrator
from scipy.optimize import least_squares

model = Model()
beta = 0.08
gamma = 0.02
phi = 0.001
N = model.N
I_cutoff = model.I_cutoff

"""
The columns are :
['Date', 'Date_YMD', 'Daily Confirmed', 'Total Confirmed', 'Daily Recovered', 'Total Recovered',
'Daily Deceased','Total Deceased']

April 15, 2020 - Feb 9, 2021 : First wave
Feb 10, 2021 - Now: Second Wave


15 April 2020 is row index 76
10 Feb 2020 is row index 377
"""
data = pd.read_csv("https://api.covid19india.org/csv/latest/case_time_series.csv")

# First wave
params_1 = [beta, gamma, phi]
data_first = np.transpose(np.array(data.values[76:376, 2:], dtype=np.int))
R_1 = data_first[3] + data_first[5]
I_1 = data_first[1] - R_1
S_1 = N - data_first[1]
Q_1 = 0


def residue_1(params):
    """function whose value to be minimised for first wave

    Parameters
    ----------
    params : array like
        [beta, gamma, phi]

    Returns
    -------
    array like
        difference between the data and fit function
    """
    I_computed = integrator(
        model.dydt_fear, len(S_1), [S_1[0], I_1[0], R_1[0], Q_1], 0, I_cutoff, params,
    )[1]
    return I_1 - I_computed


opt_1 = least_squares(residue_1, params_1)
print(opt_1.x)
plt.plot(
    integrator(
        model.dydt_fear, len(S_1), [S_1[0], I_1[0], R_1[0], Q_1], 0, I_cutoff, opt_1.x,
    )[1],
    label="first wave fit",
)
plt.plot(I_1, label="first wave data")
plt.text(
    250, 9e5, f"beta_1 = {opt_1.x[0]};\ngamma_1 = {opt_1.x[1]};\nphi_1 = {opt_1.x[2]};"
)

# Second wave
params_2 = [opt_1.x[0], opt_1.x[1]]
model.phi = opt_1.x[2]  # must set this
data_second = np.transpose(np.array(data.values[377:, 2:], dtype=np.int))
R_2 = data_second[3] + data_second[5]
I_2 = data_second[1] - R_2
S_2 = N - data_second[1]
Q_2 = 0


def residue_2(params):
    """function whose value to be minimised for second wave

    Parameters
    ----------
    params : array like
        [beta, gamma]

    Returns
    -------
    array like
        difference between the data and fit function
    """
    I_computed = integrator(
        model.dydt_fixed_fear,
        len(S_2),
        [S_2[0], I_2[0], R_2[0], Q_2],
        0,
        I_cutoff,
        params,
    )[1]
    return I_2 - I_computed


opt_2 = least_squares(residue_2, params_2)
print(opt_2.x)
plt.plot(
    integrator(
        model.dydt_fixed_fear,
        len(S_2),
        [S_2[0], I_2[0], R_2[0], Q_2],
        0,
        I_cutoff,
        opt_2.x,
    )[1],
    label="second wave fit",
)

plt.plot(I_2, label="second wave data")
plt.text(10, 3e6, f"beta_2 = {opt_2.x[0]};\ngamma_2 = {opt_2.x[1]};")
plt.legend()
plt.show()

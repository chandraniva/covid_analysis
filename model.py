import numpy as np


class Model:
    """Notations:
    S: Susceptible
    I: Infected
    R: Removed
    Q: Quarantined

    alpha: R -> S
    beta: S -> I
    gamma: I -> R
    phi: S -> Q

    I_cutoff: threshold after which S -> Q becomes significant
    R_hist[i]: recovered cases until i-th day
    """

    def __init__(self):
        self.N = 1.3526e9
        self.I_cutoff = 1e4

    def fear_function(self, I, dIdt, phi):
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
            I,
            [I < self.I_cutoff, I >= self.I_cutoff],
            [0, lambda I: phi * ((I - self.I_cutoff)) ** 2],
        )

    def dydt_fear(self, t, y, R_hist, params):
        """Create the system of eqaution with fear taken into account

        Parameters
        ----------
        t : float
            time
        y : array like
            [S, I, R, Q]
        R_hist : array like
            R in the the previous days
        params : array like
            list of parameters, [beta, gamma, phi]

        Returns
        -------
        array like
            list of differential equation
        """
        [S, I, R, Q] = y
        [beta, gamma, phi] = params

        dSdt, dRdt = 0.0, 0.0
        dIdt = beta * I * S / self.N - gamma * I

        dRdt = gamma * I
        dSdt = -beta * I * S / self.N - S * self.fear_function(I, dIdt, phi) / self.N

        dQdt = S * self.fear_function(I, dIdt, phi) / self.N

        return np.array([dSdt, dIdt, dRdt, dQdt])

    def dydt_fixed_fear(self, t, y, R_hist, params):
        """Create the system of eqaution with fear taken into account

        Parameters
        ----------
        t : float
            time
        y : array like
            [S, I, R, Q]
        R_hist : array like
            R in the the previous days
        params : array like
            list of parameters, [beta, gamma, phi]

        Returns
        -------
        array like
            list of differential equation
        """
        if self.phi is None:
            raise "phi is not defined. assign value to phi before calling this function"
        [S, I, R, Q] = y
        [beta, gamma] = params

        dSdt, dRdt = 0.0, 0.0
        dIdt = beta * I * S / self.N - gamma * I

        dRdt = gamma * I
        dSdt = (
            -beta * I * S / self.N - S * self.fear_function(I, dIdt, self.phi) / self.N
        )

        dQdt = S * self.fear_function(I, dIdt, self.phi) / self.N

        return np.array([dSdt, dIdt, dRdt, dQdt])

    def dydt_fear_waning(self, t, y, R_hist, params):
        """Create the system of eqaution with fear and waning taken into account

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

        dIdt = beta * I * S / self.N - gamma * I
        dQdt = S * self.fear_function(I, dIdt, phi) / self.N

        # waning is considered only after 6 months
        if len(R_hist) > 180:
            dRdt = gamma * I - alpha * (R_hist[-179] - R_hist[-180])
            dSdt = (
                -beta * I * S / self.N
                - S * self.fear_function(I, dIdt, phi) / self.N
                + alpha * (R_hist[-179] - R_hist[-180])
            )

        else:
            dRdt = gamma * I
            dSdt = (
                -beta * I * S / self.N - S * self.fear_function(I, dIdt, phi) / self.N
            )

        return np.array([dSdt, dIdt, dRdt, dQdt])


def integrator(func, t, y0, t0, I_cutoff, params):
    """Custom rk4 integrator

    Parameters
    ----------
    func : function
        function which is to be integrated
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
    h = 1e-1
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

        k1 = h * func(t0, y, R_hist, params)
        k2 = h * func(t0 + 0.5 * h, y + 0.5 * k1, R_hist, params)
        k3 = h * func(t0 + 0.5 * h, y + 0.5 * k2, R_hist, params)
        k4 = h * func(t0 + h, y + k3, R_hist, params)
        y = y + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

        t0 = t0 + h
        y = [y[i] for i in range(4)]

    return np.transpose((np.array(y_hist)))

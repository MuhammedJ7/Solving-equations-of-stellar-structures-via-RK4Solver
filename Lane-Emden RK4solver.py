import matplotlib.pyplot as plt
import numpy as np

class Polytrope:

    """
    A polytrope of index n.
    Attributes:
    -----------
    n : float
        The polytropic index.
    xi : numpy.ndarray
        An array containing the values of the independent variable xi.
    theta : numpy.ndarray
        An array containing the values of the dependent variable theta.
    dtheta_dxi : numpy.ndarray
        An array containing the values of the derivative of theta with respect to xi.

    """
    TOL = 1.e-12
    H0 = 1.e-2
    def __init__(self, n):
        self.n = n
        self.xi = []
        self.theta = []
        self.dtheta_dxi = []
        self._integrate()
#defining the class type 

    def _integrate(self):

        """Integrate the Lane-Emden system."""

        q = np.zeros(2, dtype=np.float64)
        xi = 0.0
        h = Polytrope.H0
        q[0] = 1.0
        q[1] = 0.0

        while h > Polytrope.TOL:

            k1 = self._lane_emden_equation(xi, q)

            k2 = self._lane_emden_equation(xi + 0.5 * h, q + 0.5 * h * k1)

            k3 = self._lane_emden_equation(xi + 0.5 * h, q + 0.5 * h * k2)

            k4 = self._lane_emden_equation(xi + h, q + h * k3)

            q += (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

            xi += h

            R_est = xi - q[0] / q[1]


            if xi + h > R_est:

                h = -q[0] / q[1]

            self.xi.append(xi)
            self.theta.append(q[0])
            self.dtheta_dxi.append(q[1])

 
        self.xi = np.array(self.xi)
        self.theta = np.array(self.theta)
        self.dtheta_dxi = np.array(self.dtheta_dxi)

    def _lane_emden_equation(self, xi, q):
        """The right-hand side of the Lane-Emden system, q' = f."""

        f = np.zeros_like(q)
        f[0] = q[1]
        if xi == 0.0:
            f[1] = (2.0 / 3.0) - q[0] ** self.n

        else:

            f[1] = -2.0 * q[1] / xi - q[0] ** self.n

        return f

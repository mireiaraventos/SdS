"""
two_waves.py
==================
Perturbed slab model generating two islands described in Mireia Raventós' TPIV 2 project.
"""
from pyoculus.fields.toroidal_bfield import ToroidalBfield
import numpy as np

class SlabBfieldTwoPert(ToroidalBfield):
    """Slab magnetic field with two perturbations"""


    def __init__(self, m, n, delta, iota_prime, x0):
        """
        Set up the problem
            
        Arguments:
            m (int): Poloidal mode number.
            n (int): Toroidal mode number.
            delta (float): Perturbation amplitude.
            iota_prime (float): Magnetic shear.
            x0 (float): Resonant surface position.
        """

        super().__init__() #inheritance: call ToroidalBfield.__init__()

        self.m = m
        self.n = n
        self.delta = delta
        self.iota_prime = iota_prime
        self.x0 = x0
    

    def B(self, coords, *args):
        x = coords[0]
        y = coords[1]
        z = coords[2]

        alpha = self.m*y - self.n*z
        epsilon = self.delta*x*(1-x)
        epsilon_prime = self.delta*(1-2*x)

        Bx = self.m*epsilon*np.sin(alpha)
        By = self.iota_prime*(x-self.x0) + epsilon_prime*np.cos(alpha)
        Bz = 1.0

        return np.array([Bx, By, Bz], dtype=np.float64)

    
    def dBdX(self, coords, *args):
        x = coords[0]
        y = coords[1]
        z = coords[2]

        alpha = self.m*y - self.n*z
        epsilon = self.delta*x*(1-x)
        epsilon_prime = self.delta*(1-2*x)
        epsilon_double_prime = -2*self.delta
        
        Bx = self.m*epsilon*np.sin(alpha)
        By = self.iota_prime*(x-self.x0) + epsilon_prime*np.cos(alpha)
        Bz = 1.0

        dBxdx = self.m*epsilon_prime*np.sin(alpha)
        dBxdy = self.m**2*epsilon*np.cos(alpha)
        dBxdz = -self.m*self.n*epsilon*np.cos(alpha)

        dBydx = self.iota_prime + epsilon_double_prime*np.cos(alpha)
        dBydy = -self.m*epsilon_prime*np.sin(alpha)
        dBydz = self.n*epsilon_prime*np.sin(alpha)

        dBu = np.zeros([3, 3], dtype=np.float64)

        dBu[0, 0] = dBxdx
        dBu[0, 1] = dBxdy
        dBu[0, 2] = dBxdz

        dBu[1, 0] = dBydx
        dBu[1, 1] = dBydy
        dBu[1, 2] = dBydz

        Bu = np.array([Bx, By, Bz], dtype=np.float64)

        return Bu, dBu


    def A(self, coords, *args):
        pass

    def convert_coords(self, incoords):
        return np.array(
            [
                incoords[0],
                np.mod(incoords[1], 2.0 * np.pi),   # y is periodic
                np.mod(incoords[2], 2.0 * np.pi),   # z is periodic
            ],
            dtype=np.float64,
        )
"""
two_waves.py
==================
Perturbed slab model generating two islands described in Mireia Raventós' TPIV 2 project.
"""
from pyoculus.fields.toroidal_bfield import ToroidalBfield
import numpy as np

class SlabBfieldTwoPert(ToroidalBfield):
    """Slab magnetic field with two perturbations"""


    def __init__(self, m1, m2, n1, n2, delta1, delta2, iota_prime, x0):
        """
        Set up the problem
            
        Arguments:
            m1 (int): First poloidal mode number.
            m2 (int): Second poloidal mode number.
            n1 (int): First toroidal mode number.
            n2 (int): Second toroidal mode number.
            delta1 (float): First perturbation amplitude.
            delta2 (float): Second perturbation amplitude.
            iota_prime (float): Magnetic shear.
            x0 (float): Resonant surface position.
        """

        super().__init__() #inheritance: call ToroidalBfield.__init__()

        self.m1 = m1
        self.m2 = m2
        self.n1 = n1
        self.n2 = n2
        self.delta1 = delta1
        self.delta2 = delta2
        self.iota_prime = iota_prime
        self.x0 = x0
    

    def B(self, coords, *args):
        x = coords[0]
        y = coords[1]
        z = coords[2]

        alpha1 = self.m1*y - self.n1*z
        alpha2 = self.m2*y - self.n2*z
        epsilon1 = self.delta1*x*(2-x)/4
        epsilon2 = self.delta2*x*(2-x)/4
        epsilon_prime1 = self.delta1*(1-x)/2
        epsilon_prime2 = self.delta2*(1-x)/2

        Bx = self.m1*epsilon1*np.sin(alpha1) + self.m2*epsilon2*np.sin(alpha2)
        By = self.iota_prime*(x-self.x0) + epsilon_prime1*np.cos(alpha1) + epsilon_prime2*np.cos(alpha2)
        Bz = 1.0

        return np.array([Bx, By, Bz], dtype=np.float64)

    
    def dBdX(self, coords, *args):
        x = coords[0]
        y = coords[1]
        z = coords[2]

        alpha1 = self.m1*y - self.n1*z
        alpha2 = self.m2*y - self.n2*z
        epsilon1 = self.delta1*x*(2-x)/4
        epsilon2 = self.delta2*x*(2-x)/4
        epsilon_prime1 = self.delta1*(1-x)/2
        epsilon_prime2 = self.delta2*(1-x)/2
        epsilon_double_prime1 = -self.delta1/2
        epsilon_double_prime2 = -self.delta2/2

        Bx = self.m1*epsilon1*np.sin(alpha1) + self.m2*epsilon2*np.sin(alpha2)
        By = self.iota_prime*(x-self.x0) + epsilon_prime1*np.cos(alpha1) + epsilon_prime2*np.cos(alpha2)
        Bz = 1.0

        dBxdx = self.m1*epsilon_prime1*np.sin(alpha1) + self.m2*epsilon_prime2*np.sin(alpha2)
        dBxdy = self.m1**2*epsilon1*np.cos(alpha1) + self.m2**2*epsilon2*np.cos(alpha2)
        dBxdz = -self.m1*self.n1*epsilon1*np.cos(alpha1) - self.m2*self.n2*epsilon2*np.cos(alpha2)

        dBydx = self.iota_prime + epsilon_double_prime1*np.cos(alpha1) + epsilon_double_prime2*np.cos(alpha2)
        dBydy = -self.m1*epsilon_prime1*np.sin(alpha1) - self.m2*epsilon_prime2*np.sin(alpha2)
        dBydz = self.n1*epsilon_prime1*np.sin(alpha1) + self.n2*epsilon_prime2*np.sin(alpha2)

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
import numpy as np
import scipy.constants
import math

class Rolls:

    def __init__(self,r1,r2,q,nu1=0.3,nu2=0.3,e_1=21000,e_2=21000):
        self.r1 = r1    #半径r1
        self.r2 = r2    #半径r2
        self.nu1 = nu1  #ポアソン比1 デフォルト0.3
        self.nu2 = nu2  #ポアソン比2 デフォルト0.3
        self.e1 = e_1    #縦弾性係数1 デフォルト21000kgf/mm2
        self.e2 = e_2    #縦弾性係数2 デフォルト21000kgf/mm2
        self.q = q      #単位長さあたりの荷重　kgf/mm

    def _radius_ratio():
        return (self.r1+self.r2) / (self.r1*self.r2)

    def _poisson_e():
        return (1-self.nu1**2) / self.e1, (1-self.nu2**2) / self.e2

    def get_press_max():
        k = self._radius_ratio()
        po_e1, po_e2 = self._poisson_e()
        po2 = 1/np.pi * k * self.q / (po_e1+po_e2)
        return math.sqrt(po2)

    def get_touchwide():
        k_rev = 1/self._radius_ratio()
        po_e1, po_e2 = self._poisson_e()
        b2 = 4/np.pi * k_rev * (po_e1+po_e2) * self.q
        return math.sqrt(b2)

    def get_press_x(x):
        po = self.get_press_max()
        b = self.get_touchwide()
        return po*math.sqrt(1-x**2/b**2)

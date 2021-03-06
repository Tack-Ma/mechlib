import numpy as np
import scipy.constants
import math

class Rolls:

    def __init__(self,r1,r2,q,nu1=0.3,nu2=0.3,e_1=21000,e_2=21000):
        self.r1 = r1    #半径r1 mm
        self.r2 = r2    #半径r2 mm
        self.nu1 = nu1  #ポアソン比1 デフォルト0.3
        self.nu2 = nu2  #ポアソン比2 デフォルト0.3
        self.e1 = e_1    #縦弾性係数1 デフォルト21000kgf/mm2
        self.e2 = e_2    #縦弾性係数2 デフォルト21000kgf/mm2
        self.q = q      #単位長さあたりの荷重　kgf/mm

    def _radius_ratio(self):
        return (self.r1/self.r2+1) / self.r1

    def _poisson_e(self):
        return (1-self.nu1**2) / self.e1, (1-self.nu2**2) / self.e2

    def get_press_max(self):
        k = self._radius_ratio()
        po_e1, po_e2 = self._poisson_e()
        po2 = 1/np.pi * k * self.q / (po_e1+po_e2)
        return math.sqrt(po2)

    def get_touchwide(self):
        k_rev = 1/self._radius_ratio()
        po_e1, po_e2 = self._poisson_e()
        b2 = 4/np.pi * k_rev * (po_e1+po_e2) * self.q
        return 2*math.sqrt(b2)

    def get_press_x(self,x):
        po = self.get_press_max()
        b = self.get_touchwide()
        return po*math.sqrt(1-x**2/b**2)

    def get_press_sameenu(self):
        k = self._radius_ratio()
        return 0.418*math.sqrt(self.q * self.e1 *k)

    def get_touchwide_sameenu(self):
        k_rev = 1/self._radius_ratio()
        return 2*1.522*math.sqrt(self.q /self.e1 * k_rev)

class RollsByLoad(Rolls):
    """"""
    def __init__(self, r1, r2, length, load,
                 nu1=0.3,nu2=0.3,e_1=21000,e_2=21000):
        q = load/length
        super().__init__(r1,r2,q,nu1,nu2,e_1,e_2)

class RolltoPlane(Rolls):

    def __init__(self,ro,q,nu1=0.3,nu2=0.3,e_1=21000,e_2=21000):
        super().__init__(ro,np.inf,q,nu1,nu2,e_1,e_2)

class RolltoPlaneByLoad(RolltoPlane):

    def __init__(self,ro,length,load,nu1=0.3,nu2=0.3,e_1=21000,e_2=21000):
        q = load/length
        super().__init__(ro,q,nu1,nu2,e_1,e_2)

class RolltoConcave(Rolls):
    """docstring for RolltoConcave.Rolls"""
    def __init__(self, r1,r2,q,nu1=0.3,nu2=0.3,e_1=21000,e_2=21000):
        super().__init__(r1,-1*r2,q,nu1,nu2,e_1,e_2)

class RolltoConcaveByLoad(RolltoConcave):
    """docstring for RolltoConcaveByLoad.RolltoConcave"""
    def __init__(self, r1, r2, length, load,
                 nu1=0.3,nu2=0.3,e_1=21000,e_2=21000):
        q = load/length
        super().__init__(r1,r2,q,nu1,nu2,e_1,e_2)

"""sample"""
if __name__ == '__main__':
    inp1 = {'r1':25, 'r2':50, 'length':50, 'load':4000}
    hz1 = RollsByLoad(**inp1)
    press1 = hz1.get_press_max()
    l = [press1,hz1.r1,hz1.r2]
    print(l)

    inp2 = {'ro':25, 'length':50, 'load':4000}
    hz2 = RolltoPlaneByLoad(**inp2)
    press2 = hz2.get_press_max()
    l = [press2,hz2.r1,hz2.r2]
    print(l)

    inp3 = {'r1':25, 'r2':50, 'length':50, 'load':4000}
    hz3 = RolltoConcaveByLoad(**inp3)
    press3 = hz3.get_press_max()
    l = [press3,hz3.r1,hz3.r2]
    print(l)

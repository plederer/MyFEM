import numpy as np 

class H1FiniteElement():
    def __init__(self, order = 1):
        self.order = order
        self.ndofs = int((order+2) * (order+1)/2)

    def Evaluate(self, x, y):
        phi = np.zeros(self.ndofs)
        phi[0] = x
        phi[1] = y
        phi[2] = 1-x-y

        if self.order > 1:
            phi[4] = y - x * y - y**2 # (1-x-y) * y * 4
            phi[3] = x - x * y - x**2  # (1-x-y) * x * 4
            phi[5] = x*y
        return phi

    def DEvaluate(self, x,y):
        Dphi = np.zeros((self.ndofs, 2))
        Dphi[0,0] = 1
        Dphi[0,1] = 0
        Dphi[1,0] = 0
        Dphi[1,1] = 1
        Dphi[2,0] = -1
        Dphi[2,1] = -1

        if self.order > 1:
            Dphi[4,0] = -y
            Dphi[4,1] = -x - 2 * y + 1
            Dphi[3,0] = -y - 2 * x + 1
            Dphi[3,1] = -x
            Dphi[5,0] = y
            Dphi[5,1] = x
        return Dphi
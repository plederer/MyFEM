import numpy as np 
from .fem import *
from .mesh import *

class LagrangeSpace:
    def __init__(self, mesh, dirichlet = ".*"):
        self.order = 1
        self.mesh = mesh
        self.ndof = mesh.nv
        
        self.fem = H1FiniteElement()
        self.dirichlet = dirichlet
    
        self.FreeDofs = [1] * self.ndof
        self.CalcFreeDofs()

    def CalcFreeDofs(self):
        # surf_els = [e for e in self.mesh.Boundaries(self.dirichlet).Elements()]
        for e in self.mesh.GetSurfaceElements(self.dirichlet):
            for v in e:
                self.FreeDofs[v] = 0

class GridFunction:
    def __init__(self, space):
        self.space = space
        self.vec = np.zeros(self.space.ndof)
    
    # def Solve(self, blf, lf):
        
class Laplace():
    def __init__(self, space, coeff = 1):
        self.space = space
        self.mat = np.zeros((space.ndof,space.ndof))
        self.invmat = np.zeros((space.ndof,space.ndof))

        self.coeff = coeff

    def Assemble(self):
        print("in Assemble")
        # for e in self.space.mesh.els
        #...
        
    def CalcInverse(self):
        print("in CalcInverse")
        #...

class Source():
    def __init__(self, space, load = 1):
        self.space = space
        self.vec = np.zeros((space.ndof))
        self.load = load

    def Assemble(self):
        print("in Assemble")
        #...
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
        # print("in Assemble")
        for nr, el in enumerate(self.space.mesh.els):
            lndofs = self.space.fem.ndofs
            A_T = np.zeros((lndofs,lndofs))

            ir = IntegrationRule()

            trafo = MeshTrafo(nr, self.space.mesh)

            for ip in ir.points:
                omega = ip[1]
                x = ip[0][0]
                y = ip[0][1]
                
                Dphi_hat = self.space.fem.DEvaluate(x,y) @ trafo.Dinv
                A_T += trafo.J * omega * self.coeff * (Dphi_hat @ Dphi_hat.T) 

            for i in range(lndofs):
                for j in range(lndofs):
                    self.mat[el[i],el[j]] += A_T[i,j]
            
    def CalcInverse(self):
        dofs = []
        for i in range(self.space.ndof):
            if self.space.FreeDofs[i]:
                dofs.append(i)
        # print(dofs)
        submat = self.mat[np.ix_(dofs,dofs)]
        # print(submat)
        submat_inv = np.linalg.inv(submat)
        self.invmat[np.ix_(dofs,dofs)] = submat_inv


class Source():
    def __init__(self, space, load = 1):
        self.space = space
        self.vec = np.zeros((space.ndof))
        self.load = load

    def Assemble(self):
        for nr, el in enumerate(self.space.mesh.els):
            lndofs = self.space.fem.ndofs
            F_T = np.zeros(lndofs)

            ir = IntegrationRule()
            trafo = MeshTrafo(nr, self.space.mesh)

            for ip in ir.points:
                omega = ip[1]
                xhat = ip[0][0]
                yhat = ip[0][1]
                x,y = trafo.map(xhat,yhat)
                F_T += trafo.J * omega * self.load(x,y) * self.space.fem.Evaluate(xhat,yhat)
                
            for i in range(lndofs):
                    self.vec[el[i]] += F_T[i]

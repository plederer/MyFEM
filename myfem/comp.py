import numpy as np 
from .fem import *
from .mesh import *

class LagrangeSpace:
    def __init__(self, mesh, order = 1, dirichlet = ""):
        self.order = order
        self.mesh = mesh
        self.ndof = mesh.nv
        if order > 1:
            self.ndof += self.mesh.nedge
        
        self.fem = H1FiniteElement(order = self.order)
        self.dirichlet = dirichlet
    
        self.FreeDofs = [1] * self.ndof
        self.CalcFreeDofs()

    def GetDofs(self, elnr):
        dofs = list(self.mesh.els[elnr])
        
        if self.order > 1:
            for ed in self.mesh.edges[elnr]:
                dofs.append(ed+self.mesh.nv)

        return dofs

    def CalcFreeDofs(self):
        for e in self.mesh.GetSurfaceElements(self.dirichlet):
            for v in e.vertices:
                self.FreeDofs[v.nr] = 0
            if self.order > 1:
                self.FreeDofs[self.mesh.nv + e.elementnode.nr] = 0


class GridFunction:
    def __init__(self, space):
        self.space = space
        self.vec = np.zeros(self.space.ndof)
    
    def SetBND(self, coeff, bnd):
        for e in self.space.mesh.GetSurfaceElements(bnd):
            for v in e.vertices:
                p = self.space.mesh.pnts[v.nr]
                self.vec[v.nr] = coeff(p[0], p[1])
            if self.space.order > 1:
                meanp = [0,0]
                val = 0
                for v in e.vertices:
                    meanp[0] += self.space.mesh.pnts[v.nr][0] * 0.5
                    meanp[1] += self.space.mesh.pnts[v.nr][1] * 0.5
                    val += 0.5 * self.vec[v.nr]
                val -= coeff(meanp[0], meanp[1])

                self.vec[self.space.mesh.nv + e.elementnode.nr] = -val * 4


    def Evaluate(self, x, y, enr):
        trafo = MeshTrafo(enr, self.space.mesh)

        xhat,yhat = trafo.mapinv(x,y)
        uu = self.space.fem.Evaluate(xhat,yhat)
        u = 0
        
        dofs = self.space.GetDofs(enr)
        for vv in range(self.space.fem.ndofs):
            u += self.vec[dofs[vv]] * uu[vv]
            
        return u
    
    def DEvaluate(self, x, y, enr):
        trafo = MeshTrafo(enr, self.space.mesh)

        xhat, yhat = trafo.mapinv(x,y)
        Duu = self.space.fem.DEvaluate(xhat,yhat) @ trafo.Dinv
        Du = [0, 0]
        
        dofs = self.space.GetDofs(enr)
        for vv in range(self.space.fem.ndofs):
            Du[0] += self.vec[dofs[vv]] * Duu[vv,0]
            Du[1] += self.vec[dofs[vv]] * Duu[vv,1]

        return Du
        
class Laplace():
    def __init__(self, space, coeff = 1):
        self.space = space
        self.mat = np.zeros((space.ndof,space.ndof))
        self.invmat = np.zeros((space.ndof,space.ndof))

        self.coeff = coeff

    def Assemble(self):
        for nr, el in enumerate(self.space.mesh.els):
            lndofs = self.space.fem.ndofs
            A_T = np.zeros((lndofs,lndofs))

            ir = IntegrationRule()
            trafo = MeshTrafo(nr, self.space.mesh)
            
            for i, ip in enumerate(ir.points):
                omega = ir.weights[i] #ip[1]
                x = ip[0]
                y = ip[1]
                
                Dphi_hat = self.space.fem.DEvaluate(x,y) @ trafo.Dinv
                A_T += trafo.J * omega * self.coeff * (Dphi_hat @ Dphi_hat.T) 
            
            dofs = self.space.GetDofs(nr)
            for i in range(lndofs):
                for j in range(lndofs):
                    self.mat[dofs[i],dofs[j]] += A_T[i,j]
            print(A_T)
            
    def CalcInverse(self):
        dofs = []
        for i in range(self.space.ndof):
            if self.space.FreeDofs[i]:
                dofs.append(i)      
        submat = self.mat[np.ix_(dofs,dofs)]
        submat_inv = np.linalg.inv(submat)
        self.invmat[np.ix_(dofs,dofs)] = submat_inv

class Mass():
    def __init__(self, space, coeff = 1):
        self.space = space
        self.mat = np.zeros((space.ndof,space.ndof))
        self.invmat = np.zeros((space.ndof,space.ndof))

        self.coeff = coeff

    def Assemble(self):
        for nr, el in enumerate(self.space.mesh.els):
            lndofs = self.space.fem.ndofs
            A_T = np.zeros((lndofs,lndofs))

            ir = IntegrationRule()

            trafo = MeshTrafo(nr, self.space.mesh)

            for i, ip in enumerate(ir.points):
                omega = ir.weights[i] #ip[1]
                x = ip[0]
                y = ip[1]
                
                phi_hat = self.space.fem.Evaluate(x,y)
                A_T += trafo.J * omega * self.coeff * np.outer(phi_hat, phi_hat) 

            dofs = self.space.GetDofs(nr)
            for i in range(lndofs):
                for j in range(lndofs):
                    self.mat[dofs[i],dofs[j]] += A_T[i,j]
            
    def CalcInverse(self):
        dofs = []
        for i in range(self.space.ndof):
            if self.space.FreeDofs[i]:
                dofs.append(i)      
        submat = self.mat[np.ix_(dofs,dofs)]
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

            for i, ip in enumerate(ir.points):
                omega = ir.weights[i]
                xhat = ip[0]
                yhat = ip[1]
                x,y = trafo.map(xhat,yhat)
                F_T += trafo.J * omega * self.load(x,y) * self.space.fem.Evaluate(xhat,yhat)
                
            dofs = self.space.GetDofs(nr)
            for i in range(lndofs):
                    self.vec[dofs[i]] += F_T[i]

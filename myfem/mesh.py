import numpy as np 
from ngsolve import IntegrationRule as ngsIntegrationRule
from ngsolve import TRIG

class IntegrationRule:
    def __init__(self, order = 1):
        ir = ngsIntegrationRule(TRIG, order)
        self.weights = ir.weights #[0.5/3, 0.5/3, 0.5/3]
        self.points =  ir.points #[(1,0), (0,1),(0,0)]
        self.order = order

class Mesh:
    def __init__(self, ngsmesh):
        self.ngsmesh = ngsmesh
        self.nv = ngsmesh.nv
        self.nedge = ngsmesh.nedge

        self.pnts = [p.point for p in ngsmesh.vertices]
        self.els = [[v.nr for v in e.vertices] for e in ngsmesh.Elements()]
        self.edges = [[v.nr for v in e.edges] for e in ngsmesh.Elements()]

    def GetSurfaceElements(self, bnd):
        return [e for e in self.ngsmesh.Boundaries(bnd).Elements()]
    
    def GetSubDivision(self, nr, sd = 0):
        self.ev = self.els[nr]
        V0 = np.array(self.pnts[self.ev[2]])
        V1 = np.array(self.pnts[self.ev[0]])
        V2 = np.array(self.pnts[self.ev[1]])

        es = []
        if sd > 0:
            V01 = 0.5 * np.array([V1[0]+V0[0],V1[1]+V0[1]])
            V02 = 0.5 * np.array([V2[0]+V0[0],V2[1]+V0[1]])
            V12 = 0.5 * np.array([V2[0]+V1[0],V2[1]+V1[1]])

            e0 = [V0, V01, V02]
            e1 = [V01, V1, V12]
            e2 = [V02, V2, V12]
            e3 = [V01, V12, V02]
            es = [e0, e1, e2 , e3]
        else:
            es = [[V0, V1, V2]]
        return es
    
class MeshTrafo:
    def __init__(self, nr, mesh):
        self.ev = mesh.els[nr]
        
        V0 = mesh.pnts[self.ev[2]]
        V1 = mesh.pnts[self.ev[0]]
        V2 = mesh.pnts[self.ev[1]]
        
        # phi = D * x + b
        self.D = np.array([[V1[0]-V0[0],V2[0]-V0[0]],[V1[1]-V0[1],V2[1]-V0[1]]])
        self.b = np.array(V0)
        self.J = np.abs(self.D[0,0] * self.D[1,1] - self.D[1,0] * self.D[0,1])
        self.Dinv = 1/self.J * np.array([[self.D[1,1],-self.D[0,1]], [-self.D[1,0], self.D[0,0]]])
    
    def map(self, xhat, yhat):
        Xhat = np.array([xhat,yhat])
        return self.D @ Xhat + self.b
        
    def mapinv(self, x, y):
        X = np.array([x,y]) - self.b
        return self.Dinv @ X


def Integrate(f, mesh, order = 1):
    integral = 0
    for nr, el in enumerate(mesh.els):
        ir = IntegrationRule(order)
        trafo = MeshTrafo(nr, mesh)

        for i, ip in enumerate(ir.points):
            omega = ir.weights[i]
            
            xhat = ip[0]
            yhat = ip[1]
            x,y = trafo.map(xhat,yhat)
            
            integral += trafo.J * omega * f(x, y, nr)
    
    return integral


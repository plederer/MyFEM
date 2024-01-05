import numpy as np 
    
class IntegrationRule:
    def __init__(self, order = 1):
        # (x,y),omega
        # self.points = [[(1/3,1/3),1]]
        self.points = [[(1,0),1/3],
                       [(0,1),1/3],
                       [(0,0),1/3]]
        self.order = order

class Mesh:
    def __init__(self, ngsmesh):
        self.ngsmesh = ngsmesh
        self.nv = ngsmesh.nv
        self.pnts = [p.point for p in ngsmesh.vertices]
        self.els = [[v.nr for v in e.vertices] for e in ngsmesh.Elements()]

    def GetSurfaceElements(self, bnd):
        return [[v.nr for v in e.vertices] for e in self.ngsmesh.Boundaries(bnd).Elements()]
    
class MeshTrafo:
    def __init__(self, nr, mesh):
        self.ev = mesh.els[nr]
        
        V0 = mesh.pnts[self.ev[0]]
        V1 = mesh.pnts[self.ev[1]]
        V2 = mesh.pnts[self.ev[2]]
        
        # phi = D * x + b
        self.D = np.array([[V1[0]-V0[0],V2[0]-V0[0]],[V1[1]-V0[1],V2[1]-V0[1]]])
        self.b = np.array(V0)
        self.J = np.abs(self.D[0,0] * self.D[1,1] - self.D[1,0] * self.D[0,1])
        self.Dinv = 1/self.J * np.array([[self.D[1,1],-self.D[0,1]], [-self.D[1,0], self.D[0,0]]])
    
    def map(self, xhat, yhat):
        Xhat = np.array([xhat,yhat])
        return self.D @ Xhat + self.b


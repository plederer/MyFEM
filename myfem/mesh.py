

    
class IntegrationRule():
    def ___init__(self):
        self.points = np.array([1/3,1/3])
        self.weights = np.array([1])
        self.order = 1

class Mesh():
    def __init__(self, ngsmesh):
        self.ngsmesh = ngsmesh
        self.nv = ngsmesh.nv
        self.pnts = [p for p in ngsmesh.vertices]
        self.els = [[v.nr for v in e.vertices] for e in ngsmesh.Elements()]

    def GetSurfaceElements(self, bnd):
        return [[v.nr for v in e.vertices] for e in self.ngsmesh.Boundaries(bnd).Elements()]
    
class MeshTrafo:
    def __init__(self):
        self.mesh = None
        print("initi of Trafo")
    

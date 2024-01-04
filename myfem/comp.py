
class LagrangeSpace:
    def __init__(self, mesh, dirichlet = ".*"):
        self.order = 1
        self.mesh = mesh
        self.ndof = mesh.nv
        
        self.dirichlet = dirichlet

        self.FreeDofs = [1] * self.ndof
        self.CalcFreeDofs()

    
    def CalcFreeDofs(self):
        surf_els = [e for e in self.mesh.Boundaries(self.dirichlet).Elements()]
        for e in surf_els:
            for v in e.vertices:
                self.FreeDofs[v.nr] = 0

class GridFunction:
    def __init__(self, space):
        self.space = space
        self.vec = np.zeros(self.space.ndof)
    
    # def Solve(self, blf, lf):
        
    

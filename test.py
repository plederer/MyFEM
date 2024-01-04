from myfem import *
from netgen.occ import unit_square
from ngsolve import Mesh
# from netgen.meshing import *

mesh = Mesh(unit_square.GenerateMesh(maxh = 0.3))

# number of vertices

pnts = [p for p in mesh.vertices]
els = [e for e in mesh.Elements()]
surf_els = [e for e in mesh.Boundaries(".*").Elements()]

V = LagrangeSpace(mesh, dirichlet = ".*")

# print(V.FreeDofs)
uh = GridFunction(V)

# a = Laplace(V, alpha)
# a.Assemble()
# # invmat = a.CalcInverse(V.FreeDofs)

# f = Source(V, load)
# f.Assemble()

# uh.Solve(a, f)

# Draw(uh)


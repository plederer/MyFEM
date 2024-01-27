from netgen.geom2d import unit_square
from ngsolve import Mesh as ngsmesh
# from netgen.meshing import *

from netgen.occ import *
wp = WorkPlane()

rec = wp.LineTo(1,0).LineTo(0,1).LineTo(0,0).Face()

geo = OCCGeometry(rec, dim=2)

# geo = unit_square
mmesh = ngsmesh(geo.GenerateMesh(maxh = 0.2))
print(mmesh.nv)
from myfem import *
mesh = Mesh(mmesh)

# number of vertices

V = LagrangeSpace(mesh, order = 1) #, dirichlet = ".*")
# print(mmesh.nv)
# print(V.ndof)
print(V.FreeDofs)
# V.FreeDofs[5] = 1
uh = GridFunction(V)

# a = Laplace(V, coeff = 1)
a = Mass(V, coeff = 1)
a.Assemble()
print("finished assemble")
print(a.mat)
a.CalcInverse() #V.FreeDofs)
print("finished inverse")
# print(a.invmat)

# # uex = lambda x,y: 0.5 * x * (1-x) * y * (1-y)
# ff = lambda x,y: y*(1-y) + x * (1-x)
ff = lambda x,y: 1-x-y
f = Source(V, load=ff)
f.Assemble()
print(f.vec)
uh.vec = a.invmat @ f.vec
# #

# for i,e in enumerate(mesh.els):
#     print(V.GetDofs(i))

# uex = lambda x,y: 0.5 * x * (1-x) * y * (1-y)
# ud = lambda x,y: y * (1-y)
# f = np.zeros(V.ndof)
# uh.SetBND(ud, "right")
# f = -a.mat @ uh.vec
# uh.vec += a.invmat @ f



# Draw(uh)



xx=[]
yy=[]
zz=[]
dxzz = []
dyzz = []


sd = 1
for i, e in enumerate(mesh.els):
    for ei in mesh.GetSubDivision(i, sd):
        for v in ei:
            x = v[0]
            y = v[1]
            xx.append(x)
            yy.append(y)
            
            dxzz.append(uh.DEvaluate(x,y, i)[0])
            dyzz.append(uh.DEvaluate(x,y, i)[1])
            zz.append(uh.Evaluate(x,y, i))
            # zz.append(uh.vec[v])

tri_idx = []
for i in range(len(mesh.els)):
    for j in range(4):
        tri_idx.append((12 * i + 3 * j, 12 * i + 3 * j + 1, 12 * i + 3 * j + 2))

# tri_idx = [(3 * i, 3 * i + 1, 3 * i + 2) for i in range(len(mesh.els))]


import matplotlib.pyplot as plt
import matplotlib.tri as tri

fig, axs = plt.subplots(subplot_kw={"projection": "3d"})
axs.view_init(elev=45.)

axs.plot_trisurf(xx,yy, zz, triangles=tri_idx, cmap="inferno")
# print(dxzz)
# axs.plot_trisurf(xx,yy, dxzz, triangles=tri_idx, cmap="inferno")

triag = tri.Triangulation(x=[1,0,0],y=[0,1,0])
axs.triplot(triag) #, "ko-")
plt.show()

# fig, axs = plt.subplots()
# axs.tricontourf(triang, zz)



# axs.plot_trisurf(xx, yy, zz, linewidth=0.2)


# axs.plot_trisurf(triang, zz_ex, cmap="inferno")

# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
# basic geometry features (for the background mesh)
from netgen.geom2d import SplineGeometry
# features for different sub steps
import time


#################################### Solver for velocity field ############################################

def VelSolver():
    
    beta1  = CoefficientFunction(0)
    beta2  = CoefficientFunction(1)
    betavec = CoefficientFunction((beta1,beta2))
    betanorm = sqrt(beta1*beta1+beta2*beta2)
    #Draw(betanorm,mesh,"norm")
    #Draw(betavec,mesh,"vector")


    timestep = 0.1
    k = timestep/nots
    #lsetold = GridFunction(H1(mesh,order=1))
    lsetold = lsetp1
    

    for i in range(0,nots):
        
        u,v = Vh.TrialFunction(), Vh.TestFunction()
        a = BilinearForm(Vh,Symmetric = False)
        a += SymbolicBFI(1/k*u*v)
        a += SymbolicBFI(0.5*((beta1*grad(u)[0]+beta2*grad(u)[1])*v))
        #njmpu = grad(u)*n_outer-grad(u.Other())*n_outer
        #njmpv = grad(v)*n_outer-grad(v.Other())*n_outer
        #njmpl = grad(lsetold)*n_outer-grad(lsetold).Other()*n_outer
        #a += SymbolicBFI( 0.5*1*h*h*njmpu*njmpv, VOL, skeleton=True)
        
        f = LinearForm(Vh)
        f += SymbolicLFI(1/k*lsetold*v-0.5*(beta1*grad(lsetold)[0]+beta2*grad(lsetold)[1])*v)
        #f += SymbolicLFI( -0.5*1*h*h*njmpl*njmpv, VOL, skeleton=True)
        
        a.Assemble()
        f.Assemble()
        
        lsetnew = GridFunction(Vh)
        lsetnew.Set(0*x+1,BND)
        
        
        res = f.vec.CreateVector()
        res.data = f.vec - a.mat * lsetnew.vec
        lsetnew.vec.data += a.mat.Inverse(Vh.FreeDofs()) * res
        lsetold = lsetnew

    return lsetnew



######################## Geometry ##############################
square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bcs=["b","r","t","l"])

######################## Mesh ##################################
meshsize = 0.05
mesh = Mesh (square.GenerateMesh(maxh=meshsize, quad_dominated=False))

############# initial level set function (negative inside Omega) #######
levelset = 1-1*y
lsetp1 = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset,lsetp1)
dom_indicator = IfPos(lsetp1,0,1)  
n_levelset = 1.0/Norm(grad(lsetp1)) * grad(lsetp1)

############# FESpaces  (standard) ##########################
order = 1
Vh = H1(mesh, order = order, dgjumps = True, dirichlet = "b")  
Vh2 = H1(mesh, order = order, dgjumps = False, dirichlet = "b")
print("das sind die dofs")
print(Vh.ndof)
print(Vh2.ndof)
Xh = FESpace([Vh,Vh])

############ Definitions and Parameters ######################
n_outer = specialcf.normal(mesh.dim)
h = specialcf.mesh_size   
nots = 50

##################### Solve Problems #############################


Draw(lsetp1,mesh,"domains")
for i in range(0,1):
  
    lsetp1= VelSolver()
    dom_indicator = IfPos(lsetp1,0,1)
    Draw(lsetp1,mesh,"domains")
    
    
###### Testkommentar ###################









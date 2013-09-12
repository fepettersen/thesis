from dolfin import *

class Left(SubDomain):
    def inside(self,x,on_boundary):
        return near(x[0],0.0)

class Right(SubDomain):
    def inside(self,x,on_boundary):
        return near(x[0],1.0)

class Top(SubDomain):
    def inside(self,x,on_boundary):
        return near(x[1],1.0)

class Bottom(SubDomain):
    def inside(self,x,on_boundary):
        return near(x[1],0.0)

class Obstacle(SubDomain):
    def inside(self,x,on_boundary):
        return (between(x[1], (0.5, 0.7)) and between(x[0], (0.2, 1.0)))

left = Left()
top = Top()
right = Right()
bottom = Bottom()
obstacle = Obstacle()

mesh = UnitSquareMesh(64,64)

domains = CellFunction("size_t",mesh)
domains.set_all(0)
obstacle.mark(domains,1)

boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)

a0 = Constant(1.0)
a1 = Constant(0.01)
g_L = Expression("- 10*exp(- pow(x[1] - 0.5, 2))")
g_R = Constant("1.0")
f = Constant(1.0)

V = FunctionSpace(mesh, "CG", 2)
u = TrialFunction(V)
v = TestFunction(V)

bcs = [DirichletBC(V, 5.0, boundaries, 2), DirichletBC(V, 0.0, boundaries, 4)]

# Define new measures associated with the interior domains and
# exterior boundaries
dx = Measure("dx")[domains]
ds = Measure("ds")[boundaries]

# Define variational form
F = (inner(a0*grad(u), grad(v))*dx(0) + inner(a1*grad(u), grad(v))*dx(1)
     - g_L*v*ds(1) - g_R*v*ds(3)
     - f*v*dx(0) - f*v*dx(1))

# Separate left and right hand sides of equation
a, L = lhs(F), rhs(F)

# Solve problem
u = Function(V)
solve(a == L, u, bcs)

# Plot solution and gradient
plot(u, title="u")
plot(grad(u), title="Projected grad(u)")
interactive()
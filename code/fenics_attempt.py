# -*- coding: utf-8 
# Attempt at combining Fenics diffusion solver and random walk

from dolfin import *
from combine import MultiscaleSolver, np

class FenicsMultiscaleSolver(MultiscaleSolver):
	"""docstring for FenicsMultiscaleSolver"""
	# from dolfin import *

	def __init__(self, mesh,functiontype='Lagrange'):
		self.mesh = mesh
		self.dt = 0.05

		self.V = FunctionSpace(self.mesh,functiontype,1)
		self.U = TrialFunction(self.V)
		self.v = TestFunction(self.V)
		self.f = Constant(0)
		self.U0 = Constant(0)
		self.Up = interpolate(self.U0,self.V)
		print np.shape(self.Up.compute_vertex_values())
		self.L = (self.Up + dt*self.f)*self.v*dx
		self.b = None
		self.bc = DirichletBC(self.V,self.U0,u0_boundary)
		self.a = self.U*self.v*dx + self.dt*inner(grad(self.U),grad(self.v))*dx

		self.WalkSolvers = []
		self.Indeces = []
		self.IterationCounter = 0

	def MapAreaToIndex(self,area,eps=1e-14):
		xcoor = [area[0][0],area[1][0]]
		ycoor = [area[0][1],area[1][1]]
		x0boundary = compile_subdomains("near(x[0],%g) && on_boundary"%xcoor[0])
		x1boundary = compile_subdomains("near(x[0],%g) && on_boundary"%xcoor[1])
		y0boundary = compile_subdomains("near(x[1],%g) && on_boundary"%ycoor[0])
		y1boundary = compile_subdomains("near(x[1],%g) && on_boundary"%ycoor[1])
		return [x0boundary,x1boundary,y0boundary,y1boundary]

	def Solve(self):
		self.U = Function(self.V)
		A = assemble(self.a)
		self.b = assemble(self.L,tensor=self.b)
		self.bc.apply(A,self.b)
		solve(A,self.U.vector(),self.b)
		counter = 0
		for solver in self.WalkSolvers:
			x = self.Indeces[counter]
		# 	hole = self.U.vector().array()[x[0][0]:x[0][1]+1,x[1][0]:x[1][1]+1]
		# 	self.U.vector().array()[x[0][0]:x[0][1]+1,x[1][0]:x[1][1]+1] = solver.advance(hole)
		# 	counter += 1
		print self.U.compute_vertex_values()
		self.Up.assign(self.U)

alpha = 3.14
beta = 2.78
dt = 0.05
mesh = UnitSquareMesh(5,5)
V = FunctionSpace(mesh,'Lagrange',1)
f = Constant(beta -2 - 2*alpha)
u = TrialFunction(V)
v = TestFunction(V)

u0 = Expression('1+x[0]*x[0] + alpha*x[1]*x[1] + beta*t',alpha = alpha, beta = beta, t = 0)
u_1 = interpolate(u0,V)

def u0_boundary(x,on_boundary):
	return on_boundary

bc = DirichletBC(V,u0,u0_boundary)

a = u*v*dx + dt*inner(nabla_grad(u),nabla_grad(v))*dx
L = (u_1 +dt*f)*v*dx

A = assemble(a)
u = Function(V)
T = 2
t = dt
b = None

'''
print type(u_1), u_1.vector()
q =  u_1.compute_vertex_values()
numpy_U = np.zeros((6,6))
for i in xrange(len(numpy_U)):
	numpy_U[i,:] = q[6*i:6*(i+1)]
import matplotlib.pyplot as mpl
x,y = np.meshgrid(np.linspace(0,1,6),np.linspace(0,1,6))
# mpl.ion()
fig = mpl.figure()
ax = fig.add_subplot(111,projection='3d')
ax.plot_wireframe(x,y,numpy_U)
mpl.show()

while t<=T:
	b = assemble(L, tensor=b)
	u0.t = t
	bc.apply(A,b)
	solve(A,u.vector(),b)
	plot(u)
	#rescale=False
	#interactive()
	t+=dt
	u_1.assign(u)
	u_e = interpolate(u0, V)
	maxdiff = np.abs(u_e.vector().array()-u.vector().array())

'''
area = [[0.3,0.3],[0.5,0.5]]
test = FenicsMultiscaleSolver(mesh)
test.Up = interpolate(u0,test.V)
test.f = f
test.a = test.U*test.v*dx + test.dt*inner(grad(test.U),grad(test.v))*dx
test.L = (test.Up + test.dt*test.f)*test.v*dx

test.AddWalkArea(area)
test.Solve()
print test.U.set_vertex_values()
# 1 and 2D diffusion eq solver
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.axes3d import Axes3D
from walk import Walk
import time

class Diffusion:
	"""2D diffusion equation solver with random walk implementation"""
	def __init__(self,x=[0,1],y=[0,1],D=1,d=1,dt=-1,solver='FE'):
		self.d = d 				# spatial dimension
		self.x = None
		self.x0,self.x1 = x
		self.y0,self.y1 = y
		self.dt = dt
		self.D = D
		self.solver = solver
		self.t = 0


		# self.n = nx				# spatial resolution
		# self.x = np.linspace(x[0],x[1],nx+1)
		# self.dx = self.x[1]-self.x[0]
		# self.dt = self.dx**2/3.0 if dt==-1 else dt
		# self._D = D*self.dt/(self.dx**2)
		# self.solver = solver
		# self.t = 0
		# if self.d ==2:
		# 	self.y = np.linspace(y[0],y[1],nx+1)
		# 	self.dy = self.y[1]-self.y[0]
		# 	self.dt = self.dx**2/5.0 if dt==-1 else dt
		# 	self._D = D*self.dt/(self.dx**2)
		# 	self.X,self.Y = np.meshgrid(self.x,self.y)

	def advance(self,U,Up,plot=False):
		if self.d==1:
			if self.x == None:
				self.Initialize_1d()
			U = self.solve1D(U,Up,plot)
		elif self.d==2:
			if self.x == None:
				nx,ny = np.shape(U)
				self.Initialize_2d(nx,ny)
			U = self.solve2D(U,Up,plot)
		else:
			print 'error!'
			return None
		self.t += self.dt
		return U

	def solve1D(self,U,Up,plot):
		if self.solver == 'FE':
			U = self.ForwardEuler1D(U,Up)
			Up = U.copy()
			
		elif self.solver == 'BE':
			U = self.BackwardEuler1D(U,Up)
			Up = U.copy()
		return U

	def solve2D(self,U,Up,plot):
		if self.solver.upper() == 'FE':
			U = self.ForwardEuler2D(U,Up)
			Up = U.copy()
		elif self.solver.upper() == 'BE':
			U = self.BackwardEuler2D(U,Up)
			Up = U.copy()
		return U

	def ForwardEuler1D(self,U,Up):
		for i in xrange(1,self.n):
			U[i] = self._D*(Up[i+1]-2*Up[i]+Up[i-1]) +Up[i]
		U = self.boundary(U,Up)
		return U

	def ForwardEuler2D(self,U,Up):
		U[1:-1,1:-1] = self._D*(Up[2:,1:-1]-2*Up[1:-1,1:-1]+Up[:-2,1:-1]) \
		+ Up[1:-1,1:-1] + self._D*(Up[1:-1,2:]-2*Up[1:-1,1:-1]+Up[1:-1,:-2]) \
		+self.dt*self.f(self.X[1:-1,1:-1],self.Y[1:-1,1:-1],self.t)
		# U = self.boundary2D(U,Up)
		U = self.NeumannBC(U,Up)
		return U

	def boundary2D(self,U,Up):
		return U

	def NeumannBC(self,U,Up,boundary='all',a=0):
		"""Neumann boundaryconditions. du/dn = a"""
		if boundary == 'all' or boundary == 0:
			U[0,1:-1] = 2*self._D*(Up[1,1:-1]-Up[0,1:-1]) + Up[0,1:-1] +\
			 self._D*(Up[0,2:]-2*Up[0,1:-1]+Up[0,:-2])
			#The corner (0,0)
			U[0,0] = 2*self._D*(Up[1,0]-Up[0,0]) +Up[0,0] + \
			2*self._D*(Up[0,1]-Up[0,0])
		if boundary == 'all' or boundary == 1:
			U[-1,1:-1] = 2*self._D*(Up[-2,1:-1]-Up[-1,1:-1])+ Up[-1,1:-1] +\
			 self._D*(Up[-1,2:]-2*Up[-1,1:-1]+Up[-1,:-2])
			 #The corner (nx,0)
			U[-1,0] = 2*self._D*(Up[-2,0]-Up[-1,0]) +Up[-1,0] + \
			 2*self._D*(Up[-1,1]-Up[-1,0])
		if boundary == 'all' or boundary == 2:
			U[1:-1,0] = self._D*(Up[2:,0]-2*Up[1:-1,0]+Up[:-2,0]) + Up[1:-1,0] +\
			2*self._D*(Up[1:-1,1]-Up[1:-1,0])
			#The corner (0,ny)
			U[0,-1] = 2*self._D*(Up[1,-1]-Up[0,-1]) + Up[0,-1] + \
			2*self._D*(Up[0,-2]-Up[0,-1])
		if boundary == 'all' or boundary == 3:
			U[1:-1,-1] = self._D*(Up[2:,-1]-2*Up[1:-1,-1]+Up[:-2,-1]) + \
			Up[1:-1,-1] +2*self._D*(Up[1:-1,-2]-Up[1:-1,-1])
			#The corner (nx,ny)
			U[-1,-1] = 2*self._D*(Up[-2,-1]-Up[-1,-1]) +Up[-1,-1] +\
			 2*self._D*(Up[-1,-2]-Up[-1,-1])
		return U

	def Initialize_1d(self,nx):
		self.x = np.linspace(self.x0,self.x1,nx)
		self.dx = self.x[1] - self.x[0]
		if self.dt == -1:
			self.dt = self.dx**2/3.0
		self._D = self.D*self.dt/(self.dx**2)

	def Initialize_2d(self,nx,ny):
		self.x = np.linspace(self.x0,self.x1,nx)
		self.y = np.linspace(self.y0,self.y1,ny)
		self.dx = self.x[1] - self.x[0]
		self.dy = self.y[1] - self.y[0]
		if self.dt == -1:
			self.dt = self.dx**2/5.0
		self.X,self.Y = np.meshgrid(self.x,self.y)
		self._D = self.D*self.dt/(self.dx*self.dx)		# Insert for proper condotion

	def DirichletBC(self,U,Up,boundary='all'):
		"""Dirichlet boundaryconditions. TO DO!"""
		pass

	def SetInitialCondition(self,U0):
		self.U0 = U0

	def SetBoundary(self,boundary_type,affected_boundary='all'):
		pass

	def boundary(self,U,Up):
		"""Neumann boundary in 1D"""
		U[0] = 2*self._D*(Up[1]-Up[0])+Up[0]
		U[-1] = 2*self._D*(Up[-2]-Up[-1]) + Up[-1]
		return U

	def f(self,x,y,t):
		"""
		Method to be overwritten in each implementation. E.g:
		def MyOwnFunction(self,x,t):
			return x[:]**(-t)
		solver = Diffusion(parameters)
		solver.f = MyOwnFunction
		"""
		return 0



twoD = True
nx = 11
ny = 11
T = 120

def f(x,y,t):
	return x*y

if __name__ == '__main__':
	if twoD:
		U = np.zeros((nx,ny))
		Up = np.zeros(np.shape(U))/np.sqrt(2)
		# Up[nx/4:3*nx/4,nx/4:3*nx/4] = 1
		Up[:nx/2,:nx/2] = 1
		solver = Diffusion(d=2)
		# solver.f = np.vectorize(f)
		mpl.ion()
		fig = mpl.figure()
		ax = fig.add_subplot(111,projection='3d')
		ax.set_autoscaley_on(False)
		for i in xrange(T):	
			U = solver.advance(U,Up)
			wframe = ax.plot_wireframe(solver.X,solver.Y,Up)
			mpl.draw()
			ax.collections.remove(wframe)
			Up = U.copy()
	else:
		U = np.zeros(11)
		Up[:6] = 0
		Up *= np.pi
		x = np.linspace(0,1,11)
		solver = Diffusion(d=1)
		im = []
		fig = mpl.figure()
		for i in xrange(50):
			im.append(mpl.plot(x,Up,'b-'))
			U = solver.advance(U,Up)
			Up = U.copy()
			# wframe = ax.plot_wireframe(solver.X,solver.Y,Up)
		ani = animation.ArtistAnimation(fig,im)
		mpl.show()

	raw_input('Press Return..')
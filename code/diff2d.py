# 1 and 2D diffusion eq solver
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.axes3d import Axes3D
from walk import Walk
import time

class Diffusion:
	"""2D diffusion equation solver with random walk implementation"""
	def __init__(self,x=[0,1],y=[0,1],nx=10,D=1,d=1,dt=-1,solver='FE'):
		self.d = d 				# spatial dimension
		self.n = nx				# spatial resolution
		self.x = np.linspace(x[0],x[1],nx+1)
		self.dx = self.x[1]-self.x[0]
		self.dt = self.dx**2/3.0 if dt==-1 else dt
		self._D = D*self.dt/(self.dx**2)
		self.solver = solver
		self.t = 0
		if self.d ==2:
			self.y = np.linspace(y[0],y[1],nx+1)
			self.dy = self.y[1]-self.y[0]
			self.dt = self.dx**2/5.0 if dt==-1 else dt
			self._D = D*self.dt/(self.dx**2)
			self.X,self.Y = np.meshgrid(self.x,self.y)

	def advance(self,U,Up,plot=False):
		if self.d==1:
			U = self.solve1D(U,Up,plot)
		elif self.d==2:
			U = self.solve2D(U,Up,plot)
		else:
			print 'error!'
			return None
		self.t += self.dt
		return U

	def solve1D(self,U,Up,plot):
		# try:
		# 	Up = self.U0
		# except AttributeError:
		# 	Up = np.zeros(self.n+1)
		# 	Up[:(self.n+1)/2] = 1
		# if plot:
		# 	fig = mpl.figure()
		# 	ims = [mpl.plot(self.x,Up,'b-')]
		# U = np.zeros(self.n+1)
		if self.solver == 'FE':
			U = self.ForwardEuler1D(U,Up)
			Up = U.copy()
			
		elif self.solver == 'BE':
			U = self.BackwardEuler1D(U,Up)
			Up = U.copy()
		# if plot:
		# 	im = mpl.plot(self.x,U,'b-')
		# 	ims.append(im)
		# 	self.animate(fig,ims)
		return U

	def solve2D(self,U,Up,plot):
		wframe=None
		# try:
		# 	Up = self.U0
		# except AttributeError:
		# 	Up = np.zeros((self.n+1,self.n+1))
		# 	Up[:(self.n+1)/2,:(self.n+1)/2] = 1
		# if plot:
		# 	mpl.ion()
		# 	fig = mpl.figure()
		# 	ax = fig.add_subplot(111,projection='3d')
		# 	self.X,self.Y = np.meshgrid(self.x,self.y)
		# 	wframe = ax.plot_wireframe(self.X,self.Y,Up)
		# 	mpl.draw()
		# U = np.zeros((self.n+1,self.n+1))
		if self.solver.upper() == 'FE':
			U = self.ForwardEuler2D(U,Up)
			Up = U.copy()
		elif self.solver.upper() == 'BE':
			U = self.BackwardEuler2D(U,Up)
			Up = U.copy()
		# if plot:
		# 	wframe = ax.plot_wireframe(self.X,self.Y,U)
		# 	if oldframe is not None:
		# 		# ims.append(im)
		# 		ax.collections.remove(oldframe)
		# 	mpl.draw()
		# 	time.sleep(1)
		return U

	def ForwardEuler1D(self,U,Up):
		for i in xrange(1,self.n):
			U[i] = self._D*(Up[i+1]-2*Up[i]+Up[i-1]) +Up[i]
		U = self.boundary(U)
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
			U[0,0] = 2*self._D*(Up[1,0]-Up[0,0]) +Up[0,0] + \
			2*self._D*(Up[0,1]-Up[0,0])
		if boundary == 'all' or boundary == 1:
			U[-1,1:-1] = 2*self._D*(Up[-2,1:-1]-Up[-1,1:-1])+ Up[-1,1:-1] +\
			 self._D*(Up[-1,2:]-2*Up[-1,1:-1]+Up[-1,:-2])
			U[-1,0] = 2*self._D*(Up[-2,0]-Up[-1,0]) +Up[-1,0] + \
			 2*self._D*(Up[-1,1]-Up[-1,0])
		if boundary == 'all' or boundary == 2:
			U[1:-1,0] = self._D*(Up[2:,0]-2*Up[1:-1,0]+Up[:-2,0]) + Up[1:-1,0] +\
			2*self._D*(Up[1:-1,1]-Up[1:-1,0])
			U[0,-1] = 2*self._D*(Up[1,-1]-Up[0,-1]) + Up[0,-1] + \
			2*self._D*(Up[0,-2]-Up[0,-1])
		if boundary == 'all' or boundary == 3:
			U[1:-1,-1] = self._D*(Up[2:,-1]-2*Up[1:-1,-1]+Up[:-2,-1]) + \
			Up[1:-1,-1] +2*self._D*(Up[1:-1,-2]-Up[1:-1,-1])
			U[-1,-1] = 2*self._D*(Up[-2,-1]-Up[-1,-1]) +Up[-1,-1] +\
			 2*self._D*(Up[-1,-2]-Up[-1,-1])
		return U

	def DirichletBC(self,U,Up,boundary='all'):
		"""Dirichlet boundaryconditions. TO DO!"""
		pass

	def SetInitialCondition(self,U0):
		self.U0 = U0

	def SetBoundary(self,boundary_type,affected_boundary='all'):
		pass

	def boundary(self,U):
		U[0] = 1
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


U = np.zeros((11,11))
Up = np.zeros(np.shape(U))
Up[:11/2,:11/2] = 1

if __name__ == '__main__':
	solver = Diffusion(d=2)
	mpl.ion()
	fig = mpl.figure()
	ax = fig.add_subplot(111,projection='3d')
	ax.set_autoscaley_on(False)
	wframe = ax.plot_wireframe(solver.X,solver.Y,Up)
	for i in xrange(50):
		mpl.draw()
		ax.collections.remove(wframe)
		U = solver.advance(U,Up)
		Up = U.copy()
		wframe = ax.plot_wireframe(solver.X,solver.Y,Up)


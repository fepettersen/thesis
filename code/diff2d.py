# 1 and 2D diffusion eq solver
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.axes3d import Axes3D
from walk import Walk
import time

class Diffusion:
	"""2D diffusion equation solver with random walk implementation"""
	def __init__(self,x=[0,1],y=[0,1],nx=10,D=1,d=1,dt=-1):
		self.d = d 				# spatial dimension
		self.n = nx
		self.x = np.linspace(x[0],x[1],nx+1)
		self.dx = self.x[1]-self.x[0]
		self.dt = self.dx**2/3.0 if dt==-1 else dt
		self._D = D*self.dt/(self.dx**2)
		if self.d ==2:
			self.y = np.linspace(y[0],y[1],nx+1)
			self.dy = self.y[1]-self.y[0]
			self.dt = self.dx**2/5.0 if dt==-1 else dt
			self._D = D*self.dt/(self.dx**2)

	def advance(self,n_t=13,plot=True):
		if self.d==1:
			self.solve1D(n_t,plot)
		elif self.d==2:
			self.solve2D(n_t,plot)
		else:
			print 'error!'
			return None

	def solve1D(self,n_t,plot):
		try:
			Up = self.U0
		except AttributeError:
			Up = np.zeros(self.n+1)
			Up[:(self.n+1)/2] = 1
		if plot:
			fig = mpl.figure()
			ims = [mpl.plot(self.x,Up,'b-')]
		U = np.zeros(self.n+1)
		for timestep in xrange(n_t):
			U = self.ForwardEuler1D(U,Up)
			Up = U.copy()
			if plot:
				im = mpl.plot(self.x,U,'b-')
				ims.append(im)
		if plot:
			self.animate(fig,ims)

	def solve2D(self,n_t,plot):
		wframe=None
		try:
			Up = self.U0
		except AttributeError:
			Up = np.zeros((self.n+1,self.n+1))
			Up[:(self.n+1)/2,:(self.n+1)/2] = 1
		if plot:
			mpl.ion()
			fig = mpl.figure()
			ax = fig.add_subplot(111,projection='3d')
			self.X,self.Y = np.meshgrid(self.x,self.y)
			wframe = ax.plot_wireframe(self.X,self.Y,Up)
			mpl.draw()
		U = np.zeros((self.n+1,self.n+1))
		for timestep in xrange(n_t):
			oldframe=wframe
			U = self.ForwardEuler2D(U,Up)
			Up = U.copy()
			if plot:
				wframe = ax.plot_wireframe(self.X,self.Y,U)
				if oldframe is not None:
					# ims.append(im)
					ax.collections.remove(oldframe)
				mpl.draw()
				time.sleep(1)

	def ForwardEuler1D(self,U,Up):
		for i in xrange(1,self.n):
			U[i] = self._D*(Up[i+1]-2*Up[i]+Up[i-1]) +Up[i]
		U = self.boundary(U)
		return U

	def ForwardEuler2D(self,U,Up):
		for i in xrange(1,self.n):
			for j in xrange(1,self.n):
				U[i,j] = self._D*(Up[i+1,j]-2*Up[i,j]+Up[i-1,j]) +\
				Up[i,j] +self._D*(Up[i,j+1]-2*Up[i,j]+Up[i,j-1])
		U = self.boundary2D(U)
		return U		

	def boundary2D(self,U):
		return U


	def SetInitialCondition(self,U0):
		self.U0 = U0

	def boundary(self,U):
		U[0] = 1
		return U


if __name__ == '__main__':
	solver = Diffusion(d=2)
	solver.advance()

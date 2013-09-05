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
		self.will_walk = False

	def solve(self,n_t=13,plot=True):
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
			U = self.advance(U,Up)
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
			U = self.advance2D(U,Up)
			Up = U.copy()
			if plot:
				wframe = ax.plot_wireframe(self.X,self.Y,U)
				if oldframe is not None:
					# ims.append(im)
					ax.collections.remove(oldframe)
				mpl.draw()
				time.sleep(1)

	def advance(self,U,Up):
		for i in xrange(1,self.n):
			if self.InAreaOfWalk(self.x[i]):
				boundary = [[self.x0,U[i-1]],[self.x1,Up[i+1]]]
				U[i] = self.walk(boundary)
				# print 'walking'
			else:
				U[i] = self._D*(Up[i+1]-2*Up[i]+Up[i-1]) +Up[i]
		U = self.boundary(U)
		return U

	def advance2D(self,U,Up):
		if self.will_walk:
			boundary = [[Up[3,3]],[Up[3,4]],[Up[4,3]],[Up[4,4]]]
			boundary = self.walk.walk(boundary)
			U[3,3] = boundary[0,0]
			U[3,4] = boundary[1,0]
			U[4,3] = boundary[2,0]
			U[4,4] = boundary[3,0]
		area = [3,4]
		# rangex = [self.x0,self.x1]
		for i in xrange(1,self.n):
			for j in xrange(1,self.n):
				if i in area or j in area: #self.InAreaOfWalk(self.x[i]):
					#insert better test!
					pass
					# boundary = [[U[3,3]],[U[3,4]],[U[4,3]],[U[4,4]]]
					# U[i] = self.walk.walk(boundary)
				else:
					U[i,j] = self._D*(Up[i+1,j]-2*Up[i,j]+Up[i-1,j]) +\
					Up[i,j] +self._D*(Up[i,j+1]-2*Up[i,j]+Up[i,j-1])
		U = self.boundary2D(U)
		return U		

	def boundary2D(self,U):
		return U

	def mesh(self,x,y,M=20,eps=1e-14):
		"""Initialize area for random walk"""
		self.x0 = x[0]
		self.x1 = x[1]
		self.y0 = y[0]
		self.y1 = y[1]
		nx0 = np.where(self.x-self.x0<eps)
		nx1 = np.where(self.x-self.x1<eps)
		ny0 = np.where(self.y-self.y0<eps)
		ny1 = np.where(self.y-self.y1<eps)
		print nx0, ny1, self.x, self.y
		self.rangex = range(nx0,nx1+1)
		self.rangey = range(ny0,ny1+1)
		self.walk = Walk([[self.x0,self.y0],[self.x1,self.y1]])
		#failsafe for coordinates
		self.M = M
		self.will_walk = True
		self.walkers = []

	def InAreaOfWalk(self,x):
		if not self.will_walk:
			return False
		elif x>self.x0 and x<self.x1:
			return True
		else:
			return False

	def walk(self,boundary,steps=100):
		x0 = self.x0; x1 = self.x1; M = self.M
		Hc = 1
		H = []
		for i in xrange(len(boundary)):
			H.append(int(boundary[i][-1]*self.M/Hc))
		for k in range(len(H)):
			i=H[k]
			for j in xrange(i):
				self.walkers.append(boundary[k][0]+\
					10**-2*np.random.standard_normal(1))
		self.nwalkers = len(self.walkers)
		counter_left = 0
		counter_right = 0
		new_walkers = []
		for walker in xrange(self.nwalkers):
			r0 = self.walkers[walker]
			s = 10**-2*np.random.standard_normal([1,steps])
			r = np.zeros([1,steps+1])	
			r[0] = r0
			r0 += np.sum(s)
			# for i in xrange(steps):
			# 	r[0,i+1] = r[0,i] +s[0,i]
			if self.HasLeft(r0):
				counter_left += 1
			elif self.HasRight(r0):
				# del(self.walkers[walker])
				counter_right += 1
			else:
				new_walkers.append(r0)
		self.walkers = new_walkers
		# print 'counter_right: ',counter_right
		return counter_right/float(M)

	def SetInitialCondition(self,U0):
		self.U0 = U0

	def boundary(self,U):
		U[0] = 1
		return U

	def HasLeft(self,pos):
		if pos<self.x0:
			return True
		else:
			return False

	def HasRight(self,pos):
		if pos>self.x1:
			return True
		else:
			return False


if __name__ == '__main__':
	solver = Diffusion(d=2)
	solver.mesh([0.3,0.4],[0.3,0.4])
	solver.solve()

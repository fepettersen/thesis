# -*- coding: utf-8 
# Module to combine (diffusion) PDE solver and random walk module

from diff2d import Diffusion
from walk import Walk, np
import time
import matplotlib.pyplot as mpl, inspect
import matplotlib.animation as animation

class MultiscaleSolver:
	"""docstring for MultiscaleSolver
	Combine a normal diffusion PDE solver with random walk model for 
	diffusion in certain areas
	Only works in 2D -- beeing rewritten to support 1d as well!"""



	def __init__(self,mesh,PdeSolver=None):
		"""mesh is the x (and y) coordinates with correct spacing, say:
		mesh = [np.linspace(x0,x1,nx),np.linspace(y0,y1,ny)] or 
		mesh = [np.linspace(x0,x1,nx)] in 1d"""
		self.mesh = np.asanyarray(mesh)
		# print mesh
		self.d = np.shape(self.mesh)[0]
		if self.d >3:
			"This should be more robust"
			self.d = 1
		nx = len(mesh[0])
		ny = 1 if self.d == 1 else len(mesh[1])
		self.WalkSolvers = []
		self.Indeces = []
		self.PdeSolver = Diffusion(d=self.d) if PdeSolver is None else PdeSolver
		self.U = np.zeros((nx,ny)) if self.d == 2 else np.zeros(nx)
		# print self.U
		self.Up = np.zeros(np.shape(self.U))
		self.IterationCounter = 0

		
	def AddWalkArea(self,area):
		self.WalkSolvers.append(Walk(area))
		self.Indeces.append(self.MapAreaToIndex(area))

	def MapAreaToIndex(self,area,eps=1e-14):
		"""Needs some work to tackle 1D as well
			Only works if x = y"""
		indeces = [[],[]]
		if self.d ==1:
			xcoor = [area[0][0],area[0][1]]
			for i in xcoor:
				for j in xrange(len(self.mesh[0])):
					if i-self.mesh[0][j]<eps:
						break
				indeces[0].append(j)
			indeces[-1] = 0	
		elif self.d==2:
			xcoor = [area[0][0],area[1][0]]
			ycoor = [area[0][1],area[1][1]]
			for i in xcoor:
				for j in xrange(len(self.mesh[0])):
					if i-self.mesh[0][j]<eps:
						break
				indeces[0].append(j)
			for i in ycoor:
				for j in xrange(len(self.mesh[-1])):
					if i-self.mesh[-1][j]<eps:
						break
				indeces[-1].append(j)
		# print indeces
		return indeces

	def getBoundary(self,counter):
		"""Only works in 2D
		Not used (as far as I can see)"""
		print inspect.stack()[0][3]		#print name of current function

		x0 = self.Indeces[counter][0][0]
		x1 = self.Indeces[counter][0][1]
		y0 = self.Indeces[counter][1][0]
		y1 = self.Indeces[counter][1][1]
		tmp = [[self.U[x0,y0],self.U[x0,y1]],[self.U[x0,y0],self.U[x1,y0]],\
		[self.U[x0,y1],self.U[x1,y1]],[self.U[x1,y0],self.U[x1,y1]]]
		# print 'tmp: ',tmp
		return tmp

	def setInitialCondition(self,U0):
		if np.shape(self.U) != np.shape(U0):
			print "Wrong shape",np.shape(U0)," of initial condition! Must match shape of mesh",np.shape(self.mesh),". Exiting"
			# U0 = U0[:np.shape(self.U)[0],:np.shape(self.U)[-1]]
			import sys
			sys.exit(0)
		self.Up = U0

	def setBoundary(self,boundary,index):
		print inspect.stack()[0][3]		#print name of current function
		x0 = self.Indeces[index][0][0]
		x1 = self.Indeces[index][0][1]
		y0 = self.Indeces[index][1][0]
		y1 = self.Indeces[index][1][1]
		rangex = range(x0,x1)
		rangey = range(y0,y1)
		for i in rangex:
			# print i-x0
			self.U[i,y0] = boundary[1,i-x0]
			self.U[i,y1] = boundary[2,i-x0]
		for j in rangey:
			self.U[x0,j] = boundary[0,j-y0]
			self.U[x1,j] = boundary[-1,j-y0]
		# print 'boundary: ',boundary
		# self.U[]

	def SaveState(self,filename=None,format='npy'):
		t = time.gmtime()
		datetime = '%02d%02d%d_%d%d'%(t.tm_mday,t.tm_mon,t.tm_year,t.tm_hour,t.tm_min)
		if filename is None:
			filename = 'Results_%s_%04d'%(datetime,self.IterationCounter)
		if format == 'npy':
			np.save(filename,self.U)
		elif format == 'txt' or format =='.txt':
			np.savetxt(filename,self.U)
		else:
			raise TypeError
		self.IterationCounter += 1

	def Solve(self):
		"""Need to support 1d as well"""
		self.U = self.PdeSolver.advance(self.U,self.Up)
		counter = 0
		for solver in self.WalkSolvers:
			x = self.Indeces[counter]
			if self.d==2:
				"Average of the two solutions"
				hole = self.U[x[0][0]:x[0][1]+1,x[1][0]:x[1][1]+1].copy()
				self.U[x[0][0]:x[0][1]+1,x[1][0]:x[1][1]+1] = 0.5*(solver.advance(hole) + self.U[x[0][0]:x[0][1]+1,x[1][0]:x[1][1]+1])
				# print self.U[x[0][0]:x[0][1]+1,x[1][0]:x[1][1]+1]
			elif self.d==1:
				"Average of the two solutions"
				hole = self.U[x[0][0]:x[0][1]+1].copy()
				# print hole
				self.U[x[0][0]:x[0][1]+1] = 0.5*(solver.advance(hole) + self.U[x[0][0]:x[0][1]+1])
				# print self.U[x[0][0]:x[0][1]+1]
			counter += 1
		self.Up = self.U.copy()

t = 0
T = 20
def setup_plot():
	mpl.ion()
	fig  = mpl.figure()
	ax = fig.add_subplot(111,projection='3d')
	ax.set_autoscaley_on(False)
	return fig,ax

if __name__ == '__main__':
	if True:
		"2D"
		nx = 11; ny =11
		Up = np.zeros((nx,ny))
		Up[(nx)/2:,(nx)/2:] = 1
		area = [[0.3,0.3],[0.5,0.5]]
		fig,ax = setup_plot()
		X,Y = np.meshgrid(np.linspace(0,1,nx),np.linspace(0,1,ny))
		mesh = [np.linspace(0,1,nx),np.linspace(0,1,ny)]
		test = MultiscaleSolver(mesh)
		test.AddWalkArea(area)
		test.setInitialCondition(Up)
		wframe = ax.plot_wireframe(X,Y,test.Up)
		mpl.draw()
		while t<T:
			ax.collections.remove(wframe)
			test.Solve()
			wframe = ax.plot_wireframe(X,Y,test.Up)
			mpl.draw()
			# time.sleep(1)
			t+=1
			# test.SaveState()
	if Fasle:
		"1D"
		im = []; fig = mpl.figure()
		n = 11
		x = np.linspace(0,1,n)
		area = [[0.3,0.5]]
		mesh = np.linspace(0,1,n)
		Up = np.zeros(n)
		Up[:n/2] = 1
		test = MultiscaleSolver([mesh])
		test.AddWalkArea(area)
		test.setInitialCondition(Up)
		im.append(mpl.plot(x,Up,'b-'))
		while t<T:
			t += 1
			print '%d of %d'%(t,T),' sum(U) = %g'%np.sum(test.Up)
			im.append(mpl.plot(x,test.U,'b-'))
			test.Solve()
		ani = animation.ArtistAnimation(fig,im)
		mpl.show()
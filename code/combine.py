# -*- coding: utf-8 
# Module to combine (diffusion) PDE solver and random walk module

from diff2d import Diffusion
from walk import Walk, np
import time
import matplotlib.pyplot as mpl
import matplotlib.animation as animation

class MultiscaleSolver:
	"""docstring for MultiscaleSolver
	Combine a normal diffusion PDE solver with random walk model for 
	diffusion in certain areas
	Only works in 2D"""



	def __init__(self,mesh,PdeSolver=Diffusion(d=2)):
		# super(MultiscaleSolver, self).__init__()
		self.mesh = mesh
		self.WalkSolvers = []
		self.Indeces = []
		self.PdeSolver = PdeSolver
		print type(PdeSolver)
		# if type(PdeSolver) == type(module()):
		# 	self.dolfin = True
		self.U = np.zeros((len(mesh),len(mesh)))
		self.Up = np.zeros(np.shape(self.U))
		self.IterationCounter = 0

		
	def AddWalkArea(self,area):
		self.WalkSolvers.append(Walk(area))
		self.Indeces.append(self.MapAreaToIndex(area))

	def MapAreaToIndex(self,area,eps=1e-14):
		"""Needs some work to tackle 1D as well
			Only works if x = y"""
		xcoor = [area[0][0],area[1][0]]
		ycoor = [area[0][1],area[1][1]]
		indeces = [[],[]]
		for i in xcoor:
			for j in xrange(len(self.mesh)):
				if i-self.mesh[j]<eps:
					break
			indeces[0].append(j)
		for i in ycoor:
			for j in xrange(len(self.mesh)):
				if i-self.mesh[j]<eps:
					break
			indeces[-1].append(j)
		print indeces
		return indeces

	def getBoundary(self,counter):
		"""Only works in 2D"""
		x0 = self.Indeces[counter][0][0]
		x1 = self.Indeces[counter][0][1]
		y0 = self.Indeces[counter][1][0]
		y1 = self.Indeces[counter][1][1]
		tmp = [[self.U[x0,y0],self.U[x0,y1]],[self.U[x0,y0],self.U[x1,y0]],\
		[self.U[x0,y1],self.U[x1,y1]],[self.U[x1,y0],self.U[x1,y1]]]
		# print 'tmp: ',tmp
		return tmp

	def setInitialCondition(self,U0):
		self.Up = U0
		if np.shape(self.U) != np.shape(self.Up):
			self.U = np.zeros(np.shape(self.Up))

	def setBoundary(self,boundary,index):
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
		# print self.U
		self.U = self.PdeSolver.advance(self.U,self.Up)
		# print self.U
		counter = 0
		for solver in self.WalkSolvers:
			x = self.Indeces[counter]
			hole = self.U[x[0][0]:x[0][1]+1,x[1][0]:x[1][1]+1]
			# print hole
			self.U[x[0][0]:x[0][1]+1,x[1][0]:x[1][1]+1] = solver.advance(hole)
			# print self.U[x[0][0]:x[0][1]+1,x[1][0]:x[1][1]+1]
			counter += 1
		self.Up = self.U.copy()

area = [[0.3,0.3],[0.5,0.5]]
mesh = np.linspace(0,1,11)
Up = np.zeros((11,11))
Up[:(11)/2,:(11)/2] = 1
X,Y = np.meshgrid(np.linspace(0,1,11),np.linspace(0,1,11))
t = 0
T = 2
def setup_plot():
	mpl.ion()
	fig  = mpl.figure()
	ax = fig.add_subplot(111,projection='3d')
	ax.set_autoscaley_on(False)
	return fig,ax

if __name__ == '__main__':
	fig,ax = setup_plot()
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
		time.sleep(1)
		t+=1
		# test.SaveState()
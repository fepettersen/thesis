# -*- coding: utf-8 
# Module to combine (diffusion) PDE solver and random walk module

from diff2d import Diffusion
from walk import Walk, np

class MultiscaleSolver:
	"""docstring for MultiscaleSolver
	Combine a normal diffusion PDE solver with random walk model for 
	diffusion in certain areas"""



	def __init__(self,mesh,PdeSolver=Diffusion()):
		# super(MultiscaleSolver, self).__init__()
		self.mesh = mesh
		self.WalkSolvers = []
		self.PdeSolver = PdeSolver
		
	def AddWalkArea(self,area):
		self.WalkSolvers.append(Walk(area))
		indeces = self.MapAreaToIndex(area)

	def MapAreaToIndex(self,area,eps=1e-14):
		xcoor = [area[0][0],area[1][0]]
		ycoor = [area[0][1],area[1][1]]
		indeces = []
		for i in xcoor:
			for j in xrange(len(self.mesh)):
				if i-self.mesh[j]<eps:
					break
			indeces.append(j)
		print indeces
		return indeces


	def Solve(self):
		U = self.PdeSolver.advance(U,Up)
		# Find boundary for walk
		for solver in self.WalkSolvers:
			boundary = U[7]
			solver.advance(boundary)

area = [[0.3,0.3],[0.5,0.5]]
mesh = np.linspace(0,1,11)
if __name__ == '__main__':
	test = MultiscaleSolver(mesh)
	test.AddWalkArea(area)
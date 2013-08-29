# simple random walk program

import numpy as np
import matplotlib.pyplot as mpl


class Walk:
	"""
	Random walk class at the moment in 1D only
	"""
	def __init__(self,boundary,thresh,d=1,M=20):
		self.Hc = thresh			#conversion factor between walkers and concentration
		self.d = d 					#dimension of the walk
		self.boundary = boundary	#the boundary of the walk with initial conditions
		self.M = M
		self.walkers = []
		self.Initialize()

	def Initialize(self):
		"""Initialize area for walk, boundary conditions etc"""
		H = []			#temporary walker number to account for position on boundary
		for i in xrange(len(self.boundary)):
			H.append(int(self.boundary[i][-1]*M/self.Hc))
		self.nwalkers = sum(H) 			#calculate from initial conditions
		for i in H:
			for j in i:
				self.walkers.append(self.boundary[i][0:self.d]+\
					10**-2*np.random.standard_normal(self.d))
		
	def setup(self):
		"""Count previous walkers and add/remove new walkers from boundary"""
		self.nwalkers = 2+int(H*M)

		pass

	def HasLeftArea(self,pos):
		if pos>area or pos<area:
			return True
		else:
			return False

	def walk(self):
		H = concentration
		counter = 0
		for walker in xrange(self.nwalkers):	
			s = 10**-2*np.random.standard_normal([self.d,steps])
			r = np.zeros([self.d,steps+1])	
			r0 = self.walkers[walker]
			r[0] = r0
			for i in xrange(steps):
				r[0,i+1] = r[0,i] +s[0,i]
#			r0 += np.sum(s)				#Alternate version (faster?)
			if self.HasLeftArea(r[0,-1]):
				counter += 1
				del(self.walkers[walker])
		self.walkers = counter*M
		return self.walkers	

if __name__ == '__main__':
	seed = False
	if seed:
		np.random.seed(1)

	d = 2				#dimension of the walks
	d_ = [1 for i in range(d)]
	walkers = 32
	for walker in range(walkers):
		N = 100			#steps
		mu = 0; sigma = 1 	#unit variance, zero mean

		r0 = 1+10**-2*np.random.standard_normal(d_) 	
		steps = 10**-2*np.random.standard_normal([d,N])
		r = np.zeros([d,N+1])	
		r[:,0] = r0			
		for i in xrange(N):
			r[:,i+1] = r[:,i] +steps[:,i]
		mpl.plot(r[0,:],r[1,:])
	mpl.show()
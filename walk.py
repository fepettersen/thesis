# -*- coding: utf-8 
# simple random walk program

import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.animation as animation

class Walk:
	"""
	Random walk class at the moment in 1D and 2D only
	"""
	def __init__(self,area,thresh=1.0,M=20):
		self.Hc = thresh			#conversion factor between walkers and concentration
		self.d = 1 if len(area)==1 else 2 	#dimension of the walk
		self.area = area	#the coordinates of the area. [[x0,x1]] in 1D, [[x0,y0],[x1,y1]] in 2D
		# print "self.d = %d"%self.d
		self.M = M 					#conversion parameter
		self.walkers = []
		self.nwalkers = len(self.walkers)
		# self.Initialize()
		
		if len(area)==1:
			self.x0,self.x1 = area[0]
		elif len(area) == 2:
			self.x0,self.y0 = area[0]
			self.x1,self.y1 = area[1]

	def Initialize(self):
		"""Initialize area for walk, boundary conditions etc"""
		H = []			#temporary walker number to account for position on boundary
		for i in xrange(len(self.boundary)):
			H.append(int(self.boundary[i][-1]*self.M/self.Hc))
		self.nwalkers = sum(H) 			#calculate from initial conditions
		for i in xrange(len(H)):
			k = H[i]
			for j in xrange(k):
				self.walkers.append(self.boundary[i][0:self.d]+\
					10**-2*np.random.standard_normal(self.d))
		
	def setup(self):
		"""Count previous walkers and add/remove new walkers from boundary"""
		self.nwalkers = 2+int(H*M)

		pass

	def HasLeftArea(self,pos):
		"""pos = [x] in 1D; 
		pos = [x,y] in 2D"""
		if self.d==1:
			if pos<self.x0 or pos>self.x1:
				return True
		elif self.d==2:
			if pos[0]>self.x1 or pos[0]<self.x0:
				return True
			if pos[1]>self.y1 or pos[1]<self.y0:
				return True
		return False


	def walk(self,concentration):
		"""
		concentration is a 2d-by-N mattrix where each column is the 
		concentration at the boundary of the walk-area at the last 
		iteration.
		"""
		# H = concentration
		# convert boundary to walkers and so on
		self.InitializeTimestep(concentration)

		steps = 100
		# counter = 0
		new_walkers = []
		for walker in xrange(self.nwalkers):	
			s = 10**-2*np.random.standard_normal([self.d,steps])
			r = np.zeros([self.d,steps+1])	
			r0 = self.walkers[walker]
			r[0] = r0
			for i in xrange(steps):
				r[0,i+1] = r[0,i] +s[0,i]
#			r0 += np.sum(s)				#Alternate version (faster?)
			if self.HasLeftArea(r[0,-1]):
				# counter += 1
				new_walkers.append(self.walkers[walker])
		self.ReturnBoundary(new_walkers)
		# self.walkers = new_walkers 		# Veldig feil!!
		# self.walkers = counter*M
		return self.walkers 

	def InitializeTimestep(self,concentration):
		for i in xrange(len(concentration)):
			self.put_walkers(concentration[i],i)
		#check if there are walkers other places than on boundary 
		#(and if it is needed)
		#calculate "gradient"

	def put_walkers(self,col,i):
		"""Only works in 2D"""
		if self.d==1:
			print 'only 2D supported'
			raise SyntaxError
		n = int(sum(col)*self.M/self.Hc) if type(col)!=type(int()) else int(col*self.M/self.Hc)	#number of walkers along column
		if i==0:
			for walker in range(n):
				self.walkers.append([self.x0 +1e-2*(0.5-np.random.uniform()),\
					self.y0 +1e-2*(self.y1-self.y0)*(0.5-np.random.uniform())])
		elif i==1:
			for walker in range(n):
				self.walkers.append([self.x0 +1e-2*(self.x1-self.x0)*(0.5-np.random.uniform()),\
					self.y0 +1e-2*(0.5-np.random.uniform())])
		elif i==2:
			for walker in range(n):
				self.walkers.append([self.x0 +1e-2*(self.x1-self.x0)*(0.5-np.random.uniform()),\
					self.y1 +1e-2*(0.5-np.random.uniform())])
		elif i==3:
			for walker in range(n):
				self.walkers.append([self.x1 +1e-2*(0.5-np.random.uniform()),\
					self.y0 +1e-2*(self.y1-self.y0)*(0.5-np.random.uniform())])
		else:
			raise SyntaxError

	def ReturnBoundary(self,walkers,concentration=None):
		"""Convert the number of walkers back into the concentration 
		at the boundary and return the boundary so it can be returned 
		from self.walk()"""
		print walkers
		dy = self.y1-self.y0
		dx = self.x1-self.x0
		d0 = dx/len(concentration[0])
		d1 = dy/len(concentration[1])
		for walker in walkers:
			#find its positioon
			pass
		pass

area = [[0.3,0.3],[0.4,0.4]]
if __name__ == '__main__':
	walk = Walk(area,1.0)
	print walk.walk([1,1])
else:
	a = 0
	b = 1
	x = np.linspace(a,b,11)
	dx=x[1]-x[0]
	d = 1				#dimension of the walks
	d_ = [1 for i in range(d)]
	walkers = 50
	C = np.zeros(len(x))
	print 'len(C) = ',len(C)
	C[0:5] = 10
	t = 0
	im = []
	fig = mpl.figure()
	while t<10:
		for i in xrange(len(C)):
			for walker in xrange(int(C[i])):
				r0 = x[i] + 1e-2*np.random.uniform()
				r0 += sum(1e-2*(0.5-np.random.uniform(size=100)))
				if i==len(C)-1:
					if r0<x[i]:
						C[i] -=1
						C[i-1] +=1
					elif r0>x[i]+dx:
						C[i] -=1
						C[0] +=1
				else:
					if r0<x[i]:
						C[i] -=1
						C[i-1] += 1
					elif r0>x[i+1]:
						C[i] -= 1
						C[i+1] +=1							
			im.append(mpl.plot(x,C,'b-'))
		t+= 1
	ani = animation.ArtistAnimation(fig,im,interval=180,blit=True)
	mpl.show()
	# for walker in range(walkers):
	# 	N = 100			#steps
	# 	r0 = 1+10**-2*np.random.standard_normal(d_)
	# 	steps = 10**-2*np.random.standard_normal([d,N])
	# 	r = np.zeros([d,N+1])	
	# 	r[:,0] = r0			
	# 	for i in xrange(N):
	# 		r[:,i+1] = r[:,i] +steps[:,i]
	# 	mpl.plot(r[0,:],r[1,:])
	# mpl.show()
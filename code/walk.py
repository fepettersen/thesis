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
		# print 'HasLeftArea(self,pos)'
		if self.d==1:
			if pos[0]<self.x0 or pos[0]>self.x1:
				return True
		elif self.d==2:
			if pos[0]>self.x1 or pos[0]<self.x0:
				return True
			if pos[1]>self.y1 or pos[1]<self.y0:
				return True
		return False


	def advance(self,concentration):
		"""
		concentration is a 2d-by-N mattrix where each column is the 
		concentration at the boundary of the walk-area at the last 
		iteration.
		"""
		# H = concentration
		# convert boundary to walkers and so on
		self.InitializeTimestep(concentration)

		steps = 100
		counter = 0
		indices = []
		new_walkers = []		#store the walkers that have left walk-area
		for walker in xrange(self.nwalkers):	
			s = 10**-2*np.random.standard_normal([self.d,steps])
			# r = np.zeros([self.d,steps+1])	
			r0 = self.walkers[walker]
			for i in xrange(self.d):
				r0[i] += np.sum(s[i])+0.1*self.gradient[i]
			if self.HasLeftArea(r0):
				indices.append(counter)
				new_walkers.append(self.walkers[walker]) 
			counter += 1
		counter = 0
		boundary = self.ReturnBoundary(new_walkers,concentration)
		for i in sorted(indices):
			# Remove walkers that have left area
			del(self.walkers[i-counter])
			counter +=1
		boundary /= self.M
		return boundary

	def InitializeTimestep(self,C):
		"C = concentration"
		nwalkers_tmp = self.nwalkers
		self.gradient = [(C[0][0]-C[-1][0]),\
		(C[1][0]-C[2][0])]
		# print gradient
		for i in xrange(len(C)):
			self.put_walkers(C[i],i)
		self.nwalkers = len(self.walkers)
		#check if there are walkers other places than on boundary 
		#(and if it is needed)
		if nwalkers_tmp == 0:
			print "this will hopefully only print once"
			# pass


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
		# print concentration
		boundary = np.zeros(np.shape(concentration))
		dy = self.y1-self.y0
		dx = self.x1-self.x0
		d0 = dx/len(concentration[0])
		d1 = dy/len(concentration[1])
		for walker in walkers:
			#find its position on the boundary
			index = self.FindPosition(walker)
			if index != -1:
				boundary[index] += 1
			# print 'index = ',index
		return boundary

	def FindPosition(self,pos):
		"""return value of -1 implies that the walker is in 
		a corner and therefore to be ignored"""
		if self.d==1:
			return 0 if pos<self.x0 else 1
		elif self.d==2:
			# print pos, " X = (%.3f,%.3f) ; Y = (%.3f,%.3f)"%(self.x0,self.x1,self.y0,self.y1) 
			if pos[0]<self.x0 and pos[1]>self.y0:
				# print 'area 0'
				if pos[1]<self.y1:
					return 0
				else:
					return -1
			elif pos[0]<self.x0 and pos[1]<self.y0:
				return -1
			elif pos[1]<self.y0 and pos[0]>self.x0:
				# print 'area 1'
				if pos[0]<self.x1:
					return 1
				else:
					return -1
			elif pos[1]<self.y0 and pos[0]<self.x0:
				return -1
			elif pos[1]>self.y1 and pos[0]>self.x0:
				# print 'area 2'
				if pos[0]<self.x1:
					return 2
				else:
					return -1
			elif pos[1]<self.y1 and pos[0]<self.x0:
				return -1
			elif pos[1]<self.y1 and pos[0]>self.x1:
				# print 'area 3'
				if pos[1]>self.y0:
					return 3
				else:
					return -1
			else:
				return -1

area = [[0.3,0.3],[0.4,0.4]]
if __name__ == '__main__':
	walk = Walk(area,1.0)
	print walk.advance([[1],[1],[0],[0]])
if False:
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
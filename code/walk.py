# -*- coding: utf-8 
# simple random walk program

import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.animation as animation

class Walk:
	"""
	Random walk class at the moment in 1D and 2D only
	"""
	def __init__(self,area,n=3,thresh=1.0,M=20):
		self.Hc = thresh			#conversion parameter between walkers and concentration
		self.d = 1 if len(area)==1 else 2 	#dimension of the walk
		self.area = area	#the coordinates of the area. [[x0,x1]] in 1D, [[x0,y0],[x1,y1]] in 2D
		self.M = M 					#conversion parameter
		self.walkers = []
		self.nwalkers = len(self.walkers)
		self.first = True
		self.factor = 1e-2			#reduces the steplength 
		if len(area)==1:
			self.x0,self.x1 = area[0]
		elif len(area) == 2:
			self.x0,self.y0 = area[0]
			self.x1,self.y1 = area[1]
		self.X,self.Y = np.meshgrid(np.linspace(self.x0,self.x1,n),\
			np.linspace(self.y0,self.y1,n))
		print self.X

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
		"""Count previous walkers and add/remove new walkers from boundary
		Currently not used"""
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
		will become the new advance function. 
		Concentration is now a matrix containing the entire are of the
		walk. This makes the coding simpler.
		"""
		a = self.InitializeTimestep(concentration)

		steps = 100
		counter = 0
		indices = []
		walkers_leaving_area = []
		if not a:
			return np.zeros(np.shape(concentration))
		for walker in xrange(self.nwalkers):	
			s = self.factor*np.random.standard_normal([self.d,steps])
			# r = np.zeros([self.d,steps+1])	
			r0 = self.walkers[walker]
			for i in xrange(self.d):
				r0[i] += np.sum(s[i])+self.factor*self.gradient[i]
			if self.HasLeftArea(r0):
				indices.append(counter)
				walkers_leaving_area.append(self.walkers[walker]) 
			counter += 1
		print len(walkers_leaving_area)
		counter = 0
		# boundary = self.ReturnBoundary(walkers_leaving_area,concentration)
		boundary = self.ReturnBoundary(self.walkers,concentration)
		for i in sorted(indices):
			# Remove walkers that have left area
			del(self.walkers[i-counter])
			counter +=1
		boundary /= self.M
		return boundary		

	
	def InitializeTimestep(self,C):
		"C = concentration"
		nwalkers_tmp = self.nwalkers
		C_tot = np.sum(C)
		N_walkers_tot = self.M*self.Hc*C_tot - nwalkers_tmp 
		print 'C_tot = ',C_tot,' N_walkers_tot = ',N_walkers_tot
		if C_tot == 0:
			return None
		self.gradient = [(C[0][0]-C[-1][0]),\
		(C[1][0]-C[2][0])]
		# self.gradient = self.CalculateGradient(C)
		# print gradient
		for i in xrange(len(C)):
			for j in xrange(len(C[i])):
				n = int((C[i][j]/C_tot)*N_walkers_tot)
				print i,j
				self.put_walkers(n,i,j)
		self.nwalkers = len(self.walkers)
		return 1


	def put_walkers(self,N,i,j):
		"""Only works in 2D
		Must be redone!"""
		if self.d==1:
			print 'only 2D supported'
			raise SyntaxError
		x = self.X[i,j]
		y = self.Y[i,j]
		for k in range(N):
			self.walkers.append([x+self.factor*(0.5-np.random.uniform()),\
				y+self.factor*(0.5-np.random.uniform())])

	def ReturnBoundary(self,walkers,concentration=None):
		"""Convert the number of walkers back into the concentration 
		at the boundary and return the boundary so it can be returned 
		from self.adnvance()"""
		# print concentration
		boundary = np.zeros(np.shape(concentration))
		dx = self.X[0,1]-self.X[0,0]
		dy = self.Y[1,0]-self.Y[0,0]
		for walker in walkers:
			#find its position on the boundary
			index = self.FindPosition(walker,dx,dy)
			boundary[index[0],index[1]] += 1
		print boundary, np.sum(boundary)
		return boundary

	def FindPosition(self,pos,dx,dy):
		"""return value of -1 implies that the walker is in 
		a corner and therefore to be ignored --> this is wrong!
		If walkers in corners are ignored we lose energy conservation.
		must be redone!
		"""
		indx = [-1,-1]
		if self.d==1:
			return 0 if pos<self.x0 else 1 # This will not be correct!
		elif self.d==2:
			if pos[0]>self.X[0,-1]:
				pos[0] = np.fmod(pos[0],(self.X[0,-1]-self.X[0,0])+dx)+self.X[0,0]
			for i in xrange(len(self.X)):
				if pos[0]-self.X[i,i]<dx:
					indx[0] = i
					break
			if pos[1]>self.Y[-1,0]:
				pos[1] = np.fmod(pos[1],(self.Y[-1,0]-self.Y[0,0])+dy)+self.Y[0,0]
				print pos[1]
			for j in xrange(len(self.Y)):
				if pos[1]-self.Y[j,j]<dy:
					indx[1] = j
					break
			print indx
			return indx

	def CalculateGradient(self,C):
		"""calculate the concentration gradient in an smart way"""
		return 1


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
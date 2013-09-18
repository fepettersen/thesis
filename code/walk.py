# -*- coding: utf-8 
# simple random walk program

import numpy as np


class Walk:
	"""
	Random walk class at the moment in 1D and 2D only
	"""
	def __init__(self,area,n=3,thresh=1.0,M=50):
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
		# print self.X


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
		Concentration is now a matrix containing the entire area of the
		walk. This makes the coding simpler.
		"""
		print 'starting, %d walkers'%self.nwalkers
		if self.d !=1:
			nx,ny = np.shape(concentration)
			self.X,self.Y = np.meshgrid(np.linspace(self.x0,self.x1,nx),\
			np.linspace(self.y0,self.y1,ny)) 	#should be only if needed
		else:
			self.x = np.linspace(self.x0,self.x1,len(concentration))
		a = self.InitializeTimestep(concentration)

		steps = 100
		counter = 0
		indices = []
		walkers_leaving_area = []
		if not a:
			return np.zeros(np.shape(concentration))
		for walker in xrange(self.nwalkers):
			"loop over all walkers"	
			s = self.factor*np.random.standard_normal([self.d,steps])
			r0 = self.walkers[walker]
			for i in xrange(steps):
				# r0[:] += s[:,i]
				r0[:] = self.checkpos(r0,s[:,i])
			# for i in xrange(self.d):
			# 	r0[i] += np.sum(s[i])+self.factor*self.gradient[i]
			if self.HasLeftArea(r0):
				indices.append(counter)
				walkers_leaving_area.append(self.walkers[walker]) 
			counter += 1
		counter = 0
		# boundary = self.ReturnBoundary(walkers_leaving_area,concentration)
		boundary = self.ReturnBoundary(self.walkers,concentration)
		print 'done, %d walkers left'%self.nwalkers
		self.walkers = []
		self.nwalkers = len(self.walkers)
		boundary /= (self.M*self.Hc)
		return boundary		

	
	def InitializeTimestep(self,C):
		"""C = concentration
		should work in 1d as well now"""
		nwalkers_tmp = self.nwalkers
		C_tot = np.sum(C)
		N_walkers_tot = int(self.M*self.Hc*C_tot - nwalkers_tmp)

		if C_tot == 0 or N_walkers_tot <= 0:
			return None
		self.gradient = self.CalculateGradient(C)
		if self.d ==1:
			for i in xrange(len(C)):
				n = int((C[i]/C_tot)*N_walkers_tot)
				self.put_walkers(n,i,0)
		if self.d ==2:
			for i in xrange(len(C)):
				for j in xrange(len(C[i])):
					n = int((C[i][j]/C_tot)*N_walkers_tot)
					# print i,j
					self.put_walkers(n,i,j)
		self.nwalkers = len(self.walkers)
		return 1


	def put_walkers(self,N,i,j):
		"""Should work in 1d as well now"""
		if self.d==1:
			x = self.x[i]
			for k in xrange(N):
				self.walkers.append([x+self.factor*(0.5-np.random.uniform())])
		elif self.d==2:
			x = self.X[i,j]
			y = self.Y[i,j]
			for k in xrange(N):
				self.walkers.append([x+self.factor*(0.5-np.random.uniform()),\
					y+self.factor*(0.5-np.random.uniform())])

	def ReturnBoundary(self,walkers,concentration=None):
		"""Convert the number of walkers back into the concentration 
		at the boundary and return the boundary so it can be returned 
		from self.advance()"""
		# print concentration
		boundary = np.zeros(np.shape(concentration))
		print 'starting %d walkers'%self.nwalkers
		if self.d ==1:
			dx = self.x[1]-self.x[0]
			cont = 0
			for walker in walkers:
				cont += 1
				index = self.FindPosition(walker,dx,0)
				boundary[index] += 1
		print 'finished, have placed %d walkers'%np.sum(boundary)
		if self.d ==2:
			dx = self.X[0,1]-self.X[0,0]
			dy = self.Y[1,0]-self.Y[0,0]
			for walker in walkers:
				#find its position on the boundary
				index = self.FindPosition(walker,dx,dy)
				boundary[index[0],index[1]] += 1
		# print boundary, np.sum(boundary)
		return boundary

	def FindPosition(self,pos,dx,dy):
		"""Finds which index the walker belongs to. 
		Implements periodic boundary conditions on the walk-area
		Should work in 1d as well now"""
		indx = [-1,-1]
		if self.d==1:
			for i in xrange(len(self.x)+1):
				# print 'i = %d, x[i] = '%i,self.x[i]
				if abs(pos-self.x[i])<dx/2.0:
					return i
		elif self.d==2:
			if pos[0]>self.X[0,-1] or pos[0]<self.X[0,0]:
				# print 'this should not happen now'
				pos[0] = np.fmod(pos[0],(self.X[0,-1]-self.X[0,0])+dx)+self.X[0,0]
				pos[0] *= ((self.X[0,-1]-self.X[0,0])+dx) if pos[0]<0 else 1
			for i in xrange(len(self.X)):
				if pos[0]-self.X[i,i]<dx:
					indx[0] = i
					break
			if pos[1]>self.Y[-1,0] or pos[1]<self.Y[0,0]:
				# print 'this should not happen now'
				pos[1] = np.fmod(pos[1],(self.Y[-1,0]-self.Y[0,0])+dy)+self.Y[0,0]
				pos[1] *= ((self.Y[-1,0]-self.Y[0,0])+dx) if pos[1]<0 else 1
				# print pos[1]
			for j in xrange(len(self.Y)):
				if pos[1]-self.Y[j,j]<dy:
					indx[1] = j
					break
			# print indx
			return indx

	def CalculateGradient(self,C):
		"""calculate the concentration gradient in an smart way
		Not so smart yet"""
		grad = np.zeros(np.shape(C))
		if self.d==1:
			for i in xrange(len(C)):
				grad[i] = (C[i]-C[i-1])/(self.x[1]-self.x[0])
		elif self.d==2:
			for i in xrange(len(C)):
				for j in xrange(len(C[i])):
					grad[i,j] = (C[i][j]-C[i-1][j])/((self.X[-1,-1]-self.X[0,0])/len(self.X))
					grad[i,j] = (C[i][j]-C[i][j-1])/((self.Y[-1,-1]-self.Y[0,0])/len(self.Y))
		return grad

	def checkpos(self,r,s,eps = 1e-5):
		"""Implements reflecting boundaries -- Has a bug!"""
		tmp = r+s
		if not self.HasLeftArea(tmp):
			return tmp
		else:
			indx = 0
			b = 1 	#find the relevant boundary
			if self.d==1:
				b = self.x0 if tmp[indx]-self.x0<0 else self.x1
			elif self.d==2:
				if tmp[0]-self.x0<=0:
					indx = 0
					b = self.x0
				elif tmp[0]-self.x1>=0:
					indx = 0
					b = self.x1
				elif tmp[1]-self.y0<=0:
					indx = 1
					b = self.y0
				else:
					indx = 1
					b = self.y1
			f = 0.5
			it = 0
			while f>eps:
				if r[indx]+f*s[indx]-b<=eps:
					r += f*s
					r -= (1-f)*s
					# if it == 0: print 'tmp = ',tmp,' r = ',r,' s = ',s,' f = %g'%f
					return r
				elif r[indx]+f*s[indx]-b <0:
					f += f/2.0
				elif r[indx]+f*s[indx]-b >0:
					f -= f/2.0
				it += 1
			return r+f*s -(1-f)*s



import matplotlib.pyplot as mpl
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D 
# area = [[0.3,0.3],[0.4,0.4]]
t =0
def setup_plot():
	mpl.ion()
	fig  = mpl.figure()
	ax = fig.add_subplot(111,projection='3d')
	ax.set_autoscaley_on(False)
	return fig,ax

if __name__ == '__main__':
	if True:
		"1D"
		n = 11
		U = np.zeros(n)
		U[:n/2] = 1
		x = np.linspace(0,1,n)
		area = [[0,1]]
		im = []
		fig  = mpl.figure()
		walk = Walk(area,1.0)
		while t<10:
			print 't = %d'%t
			im.append(mpl.plot(x,U,'b-'))
			U = walk.advance(U)	#[[1,1],[0,0]]
			t+=1
		ani = animation.ArtistAnimation(fig,im,interval=180,blit=True)
		mpl.show()
	if False:
		"2D"
		X,Y = np.meshgrid(np.linspace(0,1,11),np.linspace(0,1,11))
		U = np.zeros((11,11))
		U[:11/2,:11/2] = 1
		area = [[0,0],[1,1]]
		fig,ax = setup_plot()
		walk = Walk(area,1.0)
		wframe = ax.plot_wireframe(X,Y,U)
		mpl.draw()
		while t<10:
			print 't = %d'%t
			ax.collections.remove(wframe)
			U = walk.advance(U)	#[[1,1],[0,0]]
			t+=1
			wframe = ax.plot_wireframe(X,Y,U)
			mpl.draw()
# 1D diffusion eq solver
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.animation as animation

class Diffusion:
	"""1D diffusion equation solver with random walk implementation"""
	def __init__(self, a=0,b=1,n=100,D=1,d=1,dt=-1):
		self.d = d 				# spatial dimension
		self.a = a
		self.b = b
		self.n = n
		self.x = np.linspace(a,b,n+1)
		self.dx = self.x[1]-self.x[0]
		self.dt = self.dx**2/3.0 if dt==-1 else dt
		self._D = D*self.dt/(self.dx**2)
		self.will_walk = False

	def solve(self,n_t=13,plot=True):
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


	def mesh(self,x0,x1,M=20,eps=1e-14):
		"""Initialize area for random walk"""
		self.x0 = x0
		self.x1 = x1

		# Divide solution array into subareas???
		# Dx = x1-x0
		# Dn = int(Dx/self.dx)
		# n0 = np.where(self.x-self.x0<eps)

		self.M = M
		self.will_walk = True
		self.walkers = []

	def animate(self,fig,ims,interv=180,blt=True):
		ani = animation.ArtistAnimation(fig,ims,interval=interv,blit=blt)
		mpl.show()

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

# def walk(concentration,x0,x1,M,steps=1000):
# 	H = concentration
# 	walkers = 2+int(H*M)
# 	counter = 0
# 	for walker in range(walkers):
# 		r0 = 1+10**-2*np.random.standard_normal(1) 	
# 		s = 10**-2*np.random.standard_normal([1,steps])
# 		r = np.zeros([1,steps+1])	
# 		r[0] = r0
# 		for i in xrange(steps):
# 			r[0,i+1] = r[0,i] +s[0,i]
# 		if r[0,-1]>=x1:
# 			counter += 1
# 	return counter*M


a = 0
b = 1
n = 10
x = np.linspace(a,b,n+1)
dx = 1.0/n
dt = dx*dx/3.0
D = 1.0
D_ = dt*D/dx**2
n_t = 20			#number of timesteps
M = 100

if __name__ == '__main__':
	solver = Diffusion()
	solver.mesh(0.3,0.4)
	solver.solve(n_t=100)
'''
C = np.zeros((n_t+1,n+1))
C[0,:(n+1)/2] = 1 #np.sin(4*np.pi*x)

# Forward Euler scheme
fig = mpl.figure()
ims = []
for i in xrange(n_t):
	C[i,0] = 1
	buff = C[i,-2]
	for j in xrange(1,n):
		C[i+1,j] = D_*(C[i,j+1]-2*C[i,j]+C[i,j-1]) +C[i,j]
		#print C[i,j-1], j
	if walk(buff,1,1.1,M)>0: print 'OMG!'
	im = mpl.plot(C[i,:],'b-')
	ims.append(im)
ani = animation.ArtistAnimation(fig,ims,interval=180,blit=True)

#mpl.show()
'''
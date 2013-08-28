# 1D diffusion eq solver
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.animation as animation

class Diffusion:
	"""1D diffusion equation solver with random walk implementation"""
	def __init__(self, a=0,b=1,n=100,D=1,dt=-1):
		self.a = a
		self.b = b
		self.n = n
		self.x = np.linspace(a,b,n+1)
		self.dx = self.x[1]-self.x[0]
		self.dt = self.dx**2/3.0 if dt==-1 else dt
		self._D = self.D*self.dt/(self.dx**2)

	def solve(self,n_t=10):
		for timestep in xrange(nt):
			U = self.advance(Up)
			Up = U.copy()

	def advance(self,Up):
		


	def mesh(self):
		pass

class Solver:
	def __init__(self):
		"""
		Need: mesh, problem, boundary conditions, initial conditions
		"""
		pass

	def solve(self):
		"call advance one timestep n times"
		pass

	def advance(self):
		#C[i+1,j] = D_*(C[i,j+1]-2*C[i,j]+C[i,j-1]) +C[i,j]
		pass

def walk(concentration,x0,x1,M,steps=1000):
	H = concentration
	walkers = 2+int(H*M)
	counter = 0
	for walker in range(walkers):
		r0 = 1+10**-2*np.random.standard_normal(1) 	
		s = 10**-2*np.random.standard_normal([1,steps])
		r = np.zeros([1,steps+1])	
		r[0] = r0
		for i in xrange(steps):
			r[0,i+1] = r[0,i] +s[0,i]
		if r[0,-1]>=x1:
			counter += 1
	return counter*M


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
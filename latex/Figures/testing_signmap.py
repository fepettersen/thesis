import numpy as np, matplotlib.pyplot as mpl
from math import sqrt



balle = 9
up = np.zeros(balle)
X = np.linspace(0,1,balle)

dx = X[1]-X[0]
dt = dx

def function(x,t):
	return np.exp(-np.pi*np.pi*t)*np.cos(np.pi*x)

def Assemble1d(D,a):
	n = len(D)
	N = n**2
	A = np.zeros((n,n))
	k=0
	A[k,k] = 1+2*a*(D[1]+D[0])
	A[k,k+1] = -2*a*(D[0]+D[1])
	k+=1
	for i in xrange(1,n-1):
		A[k,k] = 1+a*D[i+1]+2*a*D[i] +a*D[i-1]
		A[k,k+1] = -a*(D[i+1]+D[i])
		A[k,k-1] = -a*(D[i-1]+D[i])
		k+=1
	A[k,k] = 1+2*a*(D[n-2]+D[n-1])
	A[k,k-1] = -2*a*(D[n-1]+D[n-2])
	return A

def Precondition(A):
	n = int(sqrt(np.shape(A)[0]))
	N = n**2
	H = []
	D = []
	Aa = []
	i=0
	a = np.zeros((n,n))
	b = np.zeros((n,n))
	c = np.zeros((n,n))

	b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
	c[:,:] = A[i*n:(i+1)*n,(i+1)*n:(i+2)*n]

	D.append(np.linalg.inv(b))
	H.append(-1*np.dot(D[-1],c))
	for i in xrange(1,n-1):
		a[:,:] = A[i*n:(i+1)*n,(i-1)*n:i*n]
		b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
		c[:,:] = A[i*n:(i+1)*n,(i+1)*n:(i+2)*n]

		Aa.append(a.copy())
		D.append(np.linalg.inv(b+np.dot(a,H[-1])))
		H.append(-1*np.dot(D[-1],c))
		
	i=n-1

	a[:,:] = A[i*n:(i+1)*n,(i-1)*n:i*n]
	Aa.append(a.copy())
	b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
	D.append(np.linalg.inv(b+np.dot(a,H[-1])))

	return H,D,Aa
	
def modified_blockTridiag(H,D,a,Up):
	n = int(sqrt(np.shape(A)[0]))
	# print np.shape(A)
	N = n**2
	g = []
	k = np.zeros(n)
	x = np.zeros(N)

	i = 0
	k[:] = Up[i*n:(i+1)*n]
	g.append(np.dot(D[0],k))
	for i in xrange(1,n-1):
		gtmp = g[-1]
		k[:] = Up[i*n:(i+1)*n]
		g.append(np.dot(D[i],(k-np.dot(a[i-1],gtmp))))
	i = n-1
	gtmp = g[-1]
	k[:] = Up[i*n:(i+1)*n]
	g.append(np.dot(D[i],(k-np.dot(a[i-1],gtmp))))

	x[i*n:(i+1)*n] = g[i]
	for i in xrange(n-2,-1,-1):
		x[i*n:(i+1)*n] = g[i]+np.dot(H[i],x[(i+1)*n:(i+2)*n])

	return x

size = 16
mpl.rc('xtick', labelsize=size) 
mpl.rc('ytick', labelsize=size)
up = function(X,0)
signmap = np.where(up<0)

alpha = dt/dx*dx
Diff = np.ones(balle)

A = Assemble1d(Diff,alpha)


h,d,aabs = Precondition(A)

u = modified_blockTridiag(h,d,aabs,abs(up))
# print u
# exit()
tmp = np.ones(np.shape(u))
tmp[signmap] *= (-1)



leg = ['u(x,t)','abs(u(x,t))','0']
mpl.plot(X,up)
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
# mpl.rc('xtick', labelsize=14) 
# mpl.rc('ytick', labelsize=14) 
mpl.hold('on')
mpl.plot(X,abs(up),'r--o')
mpl.plot(X,np.zeros(np.shape(u)),'k--')
mpl.legend(leg,loc=0,fontsize=size)

mpl.figure()
leg = ['$|u(x,t+\Delta t)|$','$signmap\cdot|(u(t+\Delta t))|$','$u(x,t+\Delta t))$','0']
mpl.plot(X,u,'r--x')
u[signmap]*=(-1)
mpl.plot(X,u,'b-o')
u = modified_blockTridiag(h,d,aabs,up)
mpl.plot(X,u,'k-')
mpl.plot(X,np.zeros(np.shape(u)),'k--')
mpl.legend(leg,loc=0,fontsize=size)

mpl.show()
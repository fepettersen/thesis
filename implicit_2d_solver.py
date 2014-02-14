import numpy as np, matplotlib.pyplot as mpl, time
from mpl_toolkits.mplot3d.axes3d import Axes3D

def ADI_2D(Up,dt,dx,dy,D,v=0):
	'''
	Alternating Direction Implicit algorithm.
	'''
	q = 0
	def tridiag(f,a,b,c):
		N = (m-q)*(n-q)
		u = np.zeros(N)
		temp = np.zeros(N)
		btemp = b[0]
		u[0] = f[0]/btemp
		for i in xrange(1,N):
			"Forward substitution"
			temp[i] = c[i-1]/btemp;
			btemp = b[i] -a[i]*temp[i];
			u[i] = (f[i] -a[i]*u[i-1])/btemp;
		for i in xrange(N-2,0,-1):
			"Backward substitution"
			u[i] -= temp[i+1]*u[i+1];
		return u;
	dtdx2 = dt/(dx**2)
	dtdy2 = dt/(dy**2)
	m,n = np.shape(Up)
	Uptmp = np.zeros((m-q)*(n-q))
	b = np.zeros((m-q)*(n-q))
	a = np.zeros((m-q)*(n-q))
	c = np.zeros((m-q)*(n-q))
	k=0
	A = np.zeros((m*n,m*n))
	# print "alpha,beta = ", alpha,",",beta
	# for i in xrange(0,n):
	# 	if i==0:
	# 		Uptmp[k] = 2*dtdy2*(D[0,i+1]+D[0,i])*(Up[0,i+1]-Up[0,i])+Up[0,i]
	# 		a[k] = 0
	# 		c[k] = -2*dtdx2*(D[1,i]+D[0,i])
	# 	elif i==m-1:
	# 		Uptmp[k] = 2*dtdy2*(D[0,i-1]+D[0,i])*(Up[0,i-1]-Up[0,i])+Up[0,i]
	# 		a[k] = -2*dtdx2*(D[1,i]+D[0,i])
	# 		c[k] = 0
	# 	else:
	# 		a[k] = -dtdx2*(D[1,i]+D[0,i])
	# 		c[k] = -dtdx2*(D[1,i]+D[0,i])
	# 		Uptmp[k] = dtdy2*((D[0,i+1]+D[0,i])*(Up[0,i+1]-Up[0,i])-(D[0,i]+D[0,i-1])*(Up[0,i]-Up[0,i-1])) +Up[0,i]
	# 	b[k] = 1.+2*dtdx2*(D[1,i]+D[i,i])
	# 	k+=1
	# for i in xrange(1,m-1):
	# 	# Uptmp[k] =  2.*beta*Up[i,1] + (1.-2.*beta)*Up[i,0]
	# 	Uptmp[k] =  2*dtdy2*(D[i,1]+D[i,0])*(Up[i,1]-Up[i,0])+Up[i,0]
	# 	a[k] = 0
	# 	c[k] = -2*dtdx2*(D[i+1,0]+D[i,0])
	# 	b[k] = 1.+2*dtdx2*(D[i+1,0]+D[i,0])
	# 	k+=1
	# 	for j in xrange(1,n-1):
	# 		Uptmp[k] = dtdy2*((D[i,j+1]+D[i,j])*(Up[i,j+1]-Up[i,j])-(D[i,j]+D[i,j-1])*(Up[i,j]-Up[i,j-1])) +Up[i,j]
	# 		a[k] = -dtdx2*(D[i-1,j]+D[i,j])
	# 		c[k] = -dtdx2*(D[i+1,j]+D[i,j])
	# 		b[k] = 1.+dtdx2*(D[i+1,j]+D[i,j]+D[i-1,j])
	# 		k += 1
	# 	Uptmp[k] = 2*dtdy2*(D[i,n-2]+D[i,n-1])*(Up[i,n-2]-Up[i,n-1])+Up[i,n-1]
	# 	a[k] = -2*dtdx2*(D[i-1,0]+D[i,0])
	# 	c[k] = 0
	# 	b[k] = 1.+2*dtdx2*(D[i-1,0]+D[i,0])
	# 	k+=1
	# for i in xrange(0,n):
	# 	if i==0:
	# 		Uptmp[k] = 2*dtdy2*(D[m-1,i+1]+D[m-1,i])*(Up[m-1,i+1]-Up[m-1,i])+Up[m-1,i]
	# 		a[k] = 0
	# 		c[k] = -2*dtdx2*(D[m-1,i]+D[m-2,i])
	# 	elif i==m-1:
	# 		Uptmp[k] = 2*dtdy2*(D[m-1,n-2]+D[m-1,n-1])*(Up[m-1,n-2]-Up[m-1,n-1])+Up[m-1,n-1]
	# 		a[k] = -2*dtdx2*(D[m-1,i]+D[m-2,i])
	# 		c[k] = 0
	# 	else:
	# 		Uptmp[k] = dtdy2*((D[m-1,i+1]+D[m-1,j])*(Up[m-1,i+1]-Up[m-1,i])-(D[m-1,i]+D[m-1,i-1])*(Up[m-1,i]-Up[m-1,i-1])) +Up[m-1,i]
	# 		a[k] = -dtdx2*(D[m-1,i]+D[m-2,i])
	# 		c[k] = -dtdx2*(D[m-1,i]+D[m-2,i])
	# 	b[k] = 1.+2*dtdx2*(D[m-2,i]+D[m-1,i])
	# 	k+=1
	# # Everything below this comment is correct!!
	# np.fill_diagonal(A,b)
	# for i in xrange(m*n-1):
	# 	A[i,i+1] = c[i]
	# 	A[i+1,i] = a[i+1]
	# for i in xrange(m*n):
	# 	for j in xrange(m*n):
	# 		print A[i,j],",",
	# 	print ""

	#############################
	#Utmp = tridiag(Uptmp,a,b,c)##
	#############################
	# Utmp = np.linalg.solve(A,Uptmp)
	# Up = Utmp.reshape(m,n)
	k=0
	for i in xrange(0,n):
		Uptmp[k] = 2*dtdy2*(D[1,i]+D[0,i])*(Up[1,i]-Up[0,i])+Up[0,i]
		if i==0:
			a[k] = 0
			c[k] = -2*dtdx2*(D[0,1]+D[0,0])
			b[k] = 1+2*dtdy2*(D[0,1]+D[0,0])
		elif i==m-1:
			a[k] = -2*dtdx2*(D[0,i-1]+D[0,i])
			c[k] = 0
			b[k] = 1+2*dtdy2*(D[0,i-1]+D[0,i])
		else:
			a[k] = -2*dtdx2*(D[0,i-1]+D[0,i])
			c[k] = -2*dtdx2*(D[0,i+1]+D[0,i])
			b[k] = 1+dtdy2*(D[0,i+1]+2*D[0,i]+D[0,i-1])
		k+=1
	for i in xrange(1,m-1):
		for j in xrange(0,n):
			Uptmp[k] = dtdx2*((D[i+1,j]+D[i,j])*(Up[i+1,j]-Up[i,j])-(D[i,j]+D[i-1,j])*(Up[i,j]-Up[i-1,j]))+Up[i,j]
			if j==0:
				a[k] = 0
				c[k] = -2*dtdy2*(D[i,1]+D[i,0])
				b[k] = 1+2*dtdy2*(D[i,j+1]+D[i,j])
			elif j==n-1:
				a[k] = -2*dtdy2*(D[i,m-2]+D[i,m-1])
				c[k] = 0
				b[k] = 1+2*dtdy2*(D[i,j]+D[i,j-1])
			else:
				a[k] = -dtdy2*(D[i,j-1]+D[i,j])
				c[k] = -dtdy2*(D[i,j+1]+D[i,j])
				b[k] = 1+dtdy2*(D[i,j+1]+2*D[i,j]+D[i,j-1])
			k += 1

	for i in xrange(0,n):
		Uptmp[k] = 2*dtdy2*(D[m-1,i]+D[m-2,i])*(Up[m-2,i]-Up[m-1,i])+Up[m-1,i]
		if i==0:
			a[k] = 0
			c[k] = -2*dtdx2*(D[m-1,1]+D[m-1,0])
			b[k] = 1+2*dtdy2*(D[m-1,1]+D[m-1,0])
		elif i==m-1:
			a[k] = -2*dtdx2*(D[m-1,i-1]+D[m-1,i])
			c[k] = 0
			b[k] = 1+2*dtdy2*(D[m-1,i-1]+D[m-1,i])
		else:
			a[k] = -2*dtdx2*(D[m-1,i-1]+D[m-1,i])
			c[k] = -2*dtdx2*(D[m-1,i+1]+D[m-1,i])
			b[k] = 1+dtdy2*(D[m-1,i+1]+2*D[m-1,i]+D[m-1,i-1])
		k+=1
	np.fill_diagonal(A,b)
	for i in xrange(m*n-1):
		A[i,i+1] = c[i]
		A[i+1,i] = a[i+1]
	############################
	# Utmp = tridiag(Uptmp,a,b,c)#
	############################
	Utmp = np.linalg.solve(A,Uptmp)
	Up = Utmp.reshape(m,n)

	return Up

def BE_2D(Up,alpha,beta,gamma):
	'''
	Backward Euler algorithm.
	'''
	q = 0

	m,n = np.shape(Up)

	Uptmp = np.zeros((m-q)*(n-q))
	b = np.zeros((m-q)*(n-q))
	a = np.zeros((m-q)*(n-q))
	c = np.zeros((m-q)*(n-q))
	d = np.zeros((m-q)*(n-q))
	e = np.zeros((m-q)*(n-q))
	k=0
	A = np.zeros((m*n,m*n))
	for i in xrange(0,m):
		Uptmp[k] = Up[i,0]
		if i==0:
			a[k] = 0
			c[k] = -2*beta
		elif i==m-1:
			a[k] = -2*beta
			c[k] = 0
		else:
			a[k] = -beta
			c[k] = -beta
		b[k] = gamma
		k+=1
		for j in xrange(1,n-1,1):
			Uptmp[k] = Up[i,j]
			if j==0:
				a[k] = 0
				c[k] = -2.*beta
			elif j==n-1:
				a[k] = -2.*beta
				c[k] = 0
			else:
				a[k] = -beta
				c[k] = -beta
			b[k] = gamma
			k += 1
		Uptmp[k] = Up[i,n-1]
		if i==0:
			a[k] = 0
			c[k] = -2*beta
		elif i==m-1:
			a[k] = -2*beta
			c[k] = 0
		else:
			a[k] = -beta
			c[k] = -beta
		b[k] = gamma
		k+=1
	np.fill_diagonal(A,b)
	for i in xrange(m*n-1):
		A[i,i+1] = c[i]
		A[i+1,i] = a[i+1]
	for i in xrange(m*n-3):
		
		A[i,i+3] = -2*alpha
		A[i+3,i] = 0
	for i in xrange(m*n):
		for j in xrange(m*n):
			# print A[i,j]," ",
			print '%3.2f '%A[i,j],
		print " "
	
	Utmp = np.linalg.solve(A,Uptmp)
	Up = Utmp.reshape(m,n)
	return Up

def Initial(x,y):
	return np.cos(np.pi*y)*np.cos(np.pi*x)

N = 31
T = 3
U = np.zeros((N,N))
# Up = np.zeros((N,N))
x,y = np.meshgrid(np.linspace(0,1,N),np.linspace(0,1,N))
dx = 1./(N-1)
dt = dx**2/4.0
dt = 0.001
D = x+y
alpha = D*dt/(dx**2)
beta = alpha
gamma = (1+2*alpha+2*beta)

Up = Initial(x,y)

# print "dt,dx = ",dt,dx
# print alpha,beta
# for i in xrange(N):
# 	for j in xrange(N):
# 		print Up[i,j]," ",
# 	print " "
mpl.ion()
fig = mpl.figure()
ax = fig.add_subplot(111,projection='3d')
wframe = ax.plot_wireframe(x,y,Up)
ax.set_autoscaley_on(False)
mpl.draw()
def f(x,y,t):
	return np.exp(-t*np.pi**2)*np.cos(np.pi*x)*np.cos(np.pi*y)

for j in xrange(0,T):
	ax.collections.remove(wframe)
	time.sleep(1)
	U = BE_2D(Up,alpha,beta,gamma)
	# U = ADI_2D(Up,dt,dx,dx,D)
	Up = U.copy()
	wframe = ax.plot_wireframe(x,y,Up)
	# wframe = ax.plot_wireframe(x,y,f(x,y,(j+1)*dt))
	mpl.draw()

'''
A small script which generates some of the more conceptual 
plots in the thesis
'''

import numpy as np, matplotlib.pyplot as mpl

size = 16

# mpl.rc('text', usetex=True)
# mpl.rc('font', family='serif')
mpl.rc('xtick', labelsize=size) 
mpl.rc('ytick', labelsize=size) 

def initial(x,x0=0.5,sigma=0.12):
	return np.exp(-(x-x0)*(x-x0)/(2*sigma*sigma))

def f(t):
	return np.exp(-t)

def func(t,up,dt):
	return -dt*np.exp(-t) + up

N = 6
Hc = 1000
t = np.linspace(0,1,N)
dt = 1.0/(N-1)
'''
u = initial(t)
C = [int(Hc*u[i]) for i in range(N)]
mpl.plot(t,u,'b-')
mpl.ylabel('concentration, u',fontsize=size)
mpl.xlabel('x',fontsize=size)
mpl.legend('u(x)',fontsize=size)

mpl.figure()

mpl.bar(t,C,alpha=0.2,width=dt,)
mpl.ylabel('Number of walkers',fontsize=size)
mpl.xlabel('x',fontsize=size)
mpl.legend('C(x)',fontsize=size)
mpl.show()


'''
u = np.zeros(N)
u[0] = 1
for i in xrange(len(t)-1):
	u[i+1] = func(t[i],u[i],dt)

leg = ['$u^{n+1} = u^n - dt\cdot e^{-t}$','$u(t) = e^{-t}$']

mpl.plot(t,u,'b--o')
mpl.hold('on')
t = np.linspace(0,1,N**2)
mpl.plot(t,f(t),'r-')

mpl.xlabel('time,  [s]', fontsize=size)
mpl.ylabel('u(t)',fontsize=size)
# mpl.title('Error from numerical integration illustrated')
mpl.legend(leg,loc=0,fontsize=size)
mpl.show()

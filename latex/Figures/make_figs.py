'''
A small script which generates some of the more conceptual 
plots in the thesis
'''

import numpy as np, matplotlib.pyplot as mpl


def f(t):
	return np.exp(-t)

def func(t,up,dt):
	return -dt*np.exp(-t) + up

N = 6
t = np.linspace(0,1,N)
dt = 1.0/(N-1)

u = np.zeros(N)
u[0] = 1
for i in xrange(len(t)-1):
	u[i+1] = func(t[i],u[i],dt)

leg = ['$u^{n+1} = u^n - dt\cdot e^{-t}$','$u(t) = e^{-t}$']
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.rc('xtick', labelsize=14) 
mpl.rc('ytick', labelsize=14) 

mpl.plot(t,u,'b--o')
mpl.hold('on')
t = np.linspace(0,1,N**2)
mpl.plot(t,f(t),'r-')

mpl.xlabel('time,  [s]', fontsize=14)
mpl.ylabel('u(t)',fontsize=14)
mpl.title('Error from numerical integration illustrated')
mpl.legend(leg)
mpl.show()
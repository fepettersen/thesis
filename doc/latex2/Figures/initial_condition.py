import numpy as np, matplotlib.pyplot as mpl

def GaussianPulse(x,y,t=0,x0=2.5,sigma=1.0,A=2.5,Hc=15):
	return A*np.exp(-(x-x0)**2/(2*sigma**2))
	# return A*np.exp(-(x-x0)**2/(2*sigma**2)) + 1.0/((x+1.5)*Hc) +0.3*np.random.rand(len(x))

def Heaviside(x,a):
	H = np.zeros(np.shape(x))
	H[x<a] = 0
	H[x>a] = 1
	return H


N = 300
x = np.linspace(0,5,N)
T = 1000

a = 0.5
u0 = GaussianPulse(x,np.zeros(N))
u = np.zeros(N)
# u[:] = u0[:]
mpl.plot(x,u0,'b-')
mpl.fill_between(x[30:80],u0[30:80],0,facecolor='blue',alpha=0.5)
mpl.hold('on')
for t in xrange(T):
	# print 'balle'
	u[1:-1] = a*(u0[0:-2]-2*u0[1:-1] + u0[2:]) +u0[1:-1]
	u[0] =  2*a*(u0[1] - u0[0]) +u0[0]
	u[-1] =  2*a*(u0[-2] - u0[-1]) +u0[-1]
	u0 = u.copy()
mpl.plot(x,u,'r--')

mpl.rc('text',usetex=True)
mpl.xlabel('x ',fontsize=16)
mpl.ylabel(r'u',fontsize=16)
# mpl.xlabel('x   [$\mu$m]',fontsize=16)
# mpl.ylabel(r'u   [$\frac{nMol}{L}$]',fontsize=16)
# mpl.text(1, -0.6, r'Initial condition$', fontsize=16)
mpl.legend([r'$u(x,t=0\cdot\Delta t)$',r'$u(x,t=10^3\cdot\Delta t)$'],fontsize=16)
# x = np.linspace(0,1,101)
# step = Heaviside(x,0.5)
# fig = mpl.figure()
# ax = fig.add_subplot(111)
# ax.plot(x,step)
# ax.annotate('0.5', xy=(0.5, -0.1),  xycoords='data',
#                 xytext=(0.5, -0.1), textcoords='axes fraction',
#                 arrowprops=dict(facecolor='black', shrink=0.05),
#                 horizontalalignment='right', verticalalignment='bottom',
#                 )
# ax.legend(['H(x-0.5)'],loc=0)
# ax.set_ylim([-0.1,1.1])
# fig.xlabel('x')
# fig.ylabel('y')

mpl.show()
from walk import Walk
from diff2d import Diffusion
import sys, numpy as np

N = int(sys.argv[1]) if len(sys.argv)>1 else 20

Mean_U = np.zeros((11,11))
U = np.zeros(np.shape(Mean_U))
Up = np.zeros(np.shape(U))
Up[:11/2,:11/2] = 1
Upp = Up.copy()

walk = Walk([[0,0],[1,1]],1.0)

solver = Diffusion(d=2)
U = solver.advance(U,Up)

for i in xrange(N):
	print 'step %d of %d'%(i,N)
	Mean_U += 0.5*(walk.advance(Upp)+U)
	Upp = Up.copy()

Mean_U /= N
X,Y = np.meshgrid(np.linspace(0,1,11),np.linspace(0,1,11))
import matplotlib.pyplot as mpl 
fig = mpl.figure()
ax = fig.add_subplot(111,projection='3d')
ax.plot_wireframe(X,Y,Mean_U)
mpl.show()

print U-Mean_U

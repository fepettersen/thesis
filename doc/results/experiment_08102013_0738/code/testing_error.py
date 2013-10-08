from walk import Walk
from diff2d import Diffusion
import sys, numpy as np

N = int(sys.argv[1]) if len(sys.argv)>1 else 20

# Mean_U = np.zeros((11,11))
# U = np.zeros(np.shape(Mean_U))
# Up = np.zeros(np.shape(U))
# Up[:11/2,:11/2] = 1
# Upp = Up.copy()

# walk = Walk([[0,0],[1,1]],1.0)

# solver = Diffusion(d=2)
# U = solver.advance(U,Up)

# for i in xrange(N):
# 	print 'step %d of %d'%(i,N)
# 	Mean_U += 0.5*(walk.advance(Upp)+U)
# 	Upp = Up.copy()

# Mean_U /= N
# X,Y = np.meshgrid(np.linspace(0,1,11),np.linspace(0,1,11))
path = '/home/fredriep/Dropbox/uio/thesis/doc/results/experiment_07102013_0645/results/'
import glob

i = []
e = []
s = []

for excl in sorted(glob.glob(path+'Excl*')):
	e.append(np.load(excl))

for incl in sorted(glob.glob(path+'Incl*')):
	i.append(np.load(incl))

for j in xrange(len(i)):
	s.append(np.max(np.abs(i[j]-e[j])))

for k in xrange(len(s)):
	print s[k]

# import matplotlib.pyplot as mpl 
# fig = mpl.figure()
# ax = fig.add_subplot(111,projection='3d')
# ax.plot_wireframe(X,Y,Mean_U)
# mpl.show()

# print U-Mean_U

import os, sys, glob, numpy as np, matplotlib.pyplot as mpl

m = 51
n = 1
T = 50
Hc = 100

tofile = 1
x0 = 0.1
x1 = 0.2
y0 = 0.1
y1 = 0.2
result_path = 'blasdf'
filename = 'asfd'
dt = 0.001

os.system('./main_walk %d %f %f %f %f %d %d %d %s %s %f %g'%(tofile,x0,x1,y0,y1,m,n,T,result_path,filename,Hc,dt))
# exit()
mpl.ion()
for step in sorted(glob.glob('spine_*.txt')):
	tmp = np.loadtxt(step)
	if len(np.shape(tmp))==1:
		break
	for i in xrange(len(tmp)):
		mpl.plot(tmp[i,0],tmp[i,1],'r-x')
	mpl.draw()
a = raw_input('press return>')
import matplotlib.pyplot as mpl, numpy as np

infile = open('../doc/results/statistics/spine_stats.txt','r')
counter = 0
data = infile.readlines()
N = len(data)
bal = np.zeros((N,5))

for line in data:
	if counter != 0:
		line = line.split()
		for x in xrange(5):
			bal[counter,x] = float(line[x])
	counter +=1

mpl.plot(bal[:,2],bal[:,1],'rx')
mpl.xlabel('spine neck length [um]')
mpl.ylabel('diffusion time [s]')
mpl.show()
import matplotlib.pyplot as mpl, numpy as np, argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--plot_normal',default=False,action='store_true',help="--plot_normal gives individual plots rather than subplots and saves plots to latex/Figures")

args = parser.parse_args(sys.argv[1:])

infile = open('../doc/results/statistics/spine_stats.txt','r')
counter = 0
data = infile.readlines()
infile.close()

N = len(data)
stats = np.zeros((N,5))

for line in data:
	if counter != 0:
		line = line.split()
		if(len(line)>=5):
			for x in xrange(5):
				stats[counter,x] = float(line[x])
	counter +=1

ins = np.where(stats[:,1]<20)
# print stats[:,1][ins]
A = np.vstack([stats[:,1][ins], np.ones(len(stats[:,2][ins]))]).T
lstsq = np.linalg.lstsq(A,stats[:,2][ins])
a,c = lstsq[0]
x = np.linspace(min(stats[:,1][ins]),max(stats[:,1][ins]),len(stats[:,1][ins]))

ins = np.where(stats[:,4]<20)
A = np.vstack([stats[:,4][ins], np.ones(len(stats[:,2][ins]))]).T
a1,c1 = np.linalg.lstsq(A,stats[:,2][ins])[0]
x1 = np.linspace(min(stats[:,4][ins]),max(stats[:,4][ins]),len(stats[:,4][ins]))

A = np.vstack([stats[:,1], np.ones(len(stats[:,3]))]).T
a2,c2 = np.linalg.lstsq(A,stats[:,3])[0]
x2 = np.linspace(min(stats[:,1]),max(stats[:,1]),len(stats[:,3]))

long_necked_spines = np.where(stats[:,2]>=0.5)[0]
relevant_data = stats[long_necked_spines,1]

path = '../doc/latex2/Figures/'
size = 18

if args.plot_normal:
	# relative diffusion time vs neck length
	mpl.plot(stats[:,2][ins],stats[:,1][ins],'rx')
	mpl.xlabel('spine neck length [um]',fontsize=size)
	mpl.ylabel('diffusion time [s]',fontsize=size)
	mpl.hold('on')
	mpl.plot(a*x+c,x,'b-')
	mpl.legend(['observed','lstsq fit'],loc=0,fontsize=size)
	mpl.rc('xtick', labelsize=size)
	mpl.rc('ytick', labelsize=size)
	mpl.savefig(path+'spine_stats_reltime_nl.eps')
	print 'a,c = ',a,',',c
	print '---------------'

	mpl.figure()
	# total diffusion time vs neck length
	mpl.plot(stats[:,2][ins],stats[:,4][ins],'rx')
	mpl.xlabel('spine neck length [um]',fontsize=size)
	mpl.ylabel('total diffusion time [s]',fontsize=size)
	mpl.hold('on')
	mpl.plot(a1*x1+c1,x1,'b-')
	mpl.legend(['observed','lstsq fit'],loc=0,fontsize=size)
	mpl.rc('xtick', labelsize=size)
	mpl.rc('ytick', labelsize=size)
	mpl.savefig(path+'spine_stats_fulltime_nl.eps')
	print 'a,c = ',a1,',',c1
	print '----------------'

	mpl.figure()
	# relative diffusion time vs neck width
	mpl.plot(stats[:,3],stats[:,1],'rx')
	mpl.xlabel('spine neck width [um]',fontsize=size)
	mpl.ylabel('total diffusion time [s]',fontsize=size)
	mpl.hold('on')
	mpl.plot(a2*x2+c2,x2,'b-')
	mpl.legend(['observed','lstsq fit'],loc=0,fontsize=size)
	mpl.savefig(path+'spine_stats_reltime_nw.eps')

	mpl.figure()
	mpl.boxplot(relevant_data)
	mpl.title('diffusion times for spines with neck length > 0.5 um',fontsize=size)
	mpl.ylabel('diffusion time [s]',fontsize=size)
	mpl.savefig(path+'spine_stats_boxplot_reltime_longneck.eps')

else:
	mpl.subplot(2,2,1)
	# relative diffusion time vs neck length
	mpl.plot(stats[:,2][ins],stats[:,1][ins],'rx')
	mpl.xlabel('spine neck length [um]')
	mpl.ylabel('diffusion time [s]')
	mpl.hold('on')
	mpl.plot(a*x+c,x,'b-')
	mpl.legend(['observed','lstsq fit'],loc=0)
	print 'a,c = ',a,',',c
	print '---------------'

	mpl.subplot(2,2,2)
	# total diffusion time vs neck length
	mpl.plot(stats[:,2][ins],stats[:,4][ins],'rx')
	mpl.xlabel('spine neck length [um]')
	mpl.ylabel('total diffusion time [s]')
	mpl.hold('on')
	mpl.plot(a1*x1+c1,x1,'b-')
	mpl.legend(['observed','lstsq fit'],loc=0)
	print 'a,c = ',a1,',',c1
	print '----------------'

	mpl.subplot(2,2,3)
	# relative diffusion time vs neck width
	mpl.plot(stats[:,3],stats[:,1],'rx')
	mpl.xlabel('spine neck width [um]')
	mpl.ylabel('total diffusion time [s]')
	mpl.hold('on')
	mpl.plot(a2*x2+c2,x2,'b-')
	mpl.legend(['observed','lstsq fit'],loc=0)

	mpl.subplot(2,2,4)
	mpl.boxplot(relevant_data)
	mpl.title('diffusion times for spines with neck length > 0.5 um')
	mpl.ylabel('diffusion time [s]')

mpl.show()

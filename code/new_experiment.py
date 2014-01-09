# -*- coding: utf-8 
"""
 script that runs an experiment, creates a new subdirectory with a copy of 
 the relevant code and parameters, stores the results (plots) and updates 
 the thesis homepage on github with this information
"""
import os, sys, time, re, numpy as np, matplotlib.pyplot as mpl #,argparse??
import glob
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.axes3d import Axes3D

def right_split(s,delimiter):
	i = len(s)-1
	while i>0:
		if s[i] == delimiter:
			return s[:i]
		i -= 1

def setup_plot():
	mpl.ion()
	fig  = mpl.figure()
	ax = fig.add_subplot(111,projection='3d')
	return fig,ax

class Experiment:
	"""Usage:
	this_dir = right_split(os.getcwd(),'/')
	obj = Experiment(this_dir,True)
	Hc = [i for i in xrange(10)]
	obj.SetupRun(x0,x1,y0,y1,n,T)
	obj.RunDeterministic()
	for i in Hc:
		obj.RunSimulation(i)
	obj.CalculateError()
	obj.PlotError()
	obj.Finish()
	"""
	def __init__(self,path, DEBUG=False, save=True):
		self.debug = DEBUG
		this_dir = path

		t = time.gmtime()
		self.datetime = '%02d%02d%d_%02d%02d'%(t.tm_mday,t.tm_mon,t.tm_year,t.tm_hour,t.tm_min)
		self.url = 'https://raw.github.com/fepettersen/thesis/master/doc/results'+'/experiment_%s/results'%self.datetime
		self.parent_path = this_dir +'/doc/results/experiment_%s'%self.datetime
		# self.parent_path = this_dir + '/doc/results/experiment_04122013_1259_convergenceTest_combinedSimulation_2d'
		self.code_path = self.parent_path+'/code'
		self.parameter_path = self.parent_path+'/parameters'
		self.result_path = self.parent_path+ '/results'
		self.tofile = 1 if save else 0;

		os.system('mkdir %s'%self.parent_path)
		os.system('mkdir %s'%self.parameter_path)
		os.system('mkdir %s'%self.result_path)
		os.system('mkdir %s'%self.code_path)
		os.system('cp *.cpp *.h new_experiment.py %s'%self.code_path)
		
		self.walk = []
		self.no_walk = []
		self.T = 11
		self.n = 11
		self.x0 = self.x1 = self.y0 = self.y1 = 0
		self.runcounter =0
		self.samples = 0

	def compile(self):
		# os.system('g++ *.cpp -O2 -o -larmadillo -llapack -lblas main_walk')
		os.system('g++ *.cpp -o main_walk -O2 -larmadillo -llapack -lblas')

	def SetupRun(self,x0,x1,y0,y1,m,n,T,dt,filename="a"):
		self.x0 = x0; self.x1 = x1; self.y0 = y0; self.y1 = y1
		self.m = m; self.n = n; self.T = T; self.dt = dt
		self.filename = filename

	def RunDeterministic(self):
		# self.walk = [0]*self.T
		os.system('./main_walk %d %f %f %f %f %d %d %d %s %s %f %g'%(self.tofile,0,0,0,0,self.m,self.n,self.T,self.result_path,"Deterministic",0,self.dt))
		for step in sorted(glob.glob(self.result_path+'/Deterministic*.txt')):
			# self.walk.append(np.fromfile(step,sep=" "))
			self.walk.append(np.loadtxt(step))		

	def RunSimulation(self,Hc):
		tofile = self.tofile
		x0 = self.x0; x1 = self.x1; y0 = self.y0; y1 = self.y1
		m = self.m; n = self.n; T = self.T
		filename = self.filename
		result_path = self.result_path
		tmp = [0]*T
		os.system('./main_walk %d %f %f %f %f %d %d %d %s %s %f %g'%(tofile,x0,x1,y0,y1,m,n,T,result_path,filename,Hc,self.dt))
		self.runcounter += 1

	def CalculateError(self,Hc,exact=None):
		self.legends = Hc
		self.error = []
		tmp = np.zeros(self.T)
		if exact is None:
			for i in xrange(self.runcounter):
				for j in xrange(self.T):
					self.error[i][j] = np.max(np.absolute(self.walk[j]-self.no_walk[i][j]))
		else:
			if self.n <= 1:
				X = np.linspace(0,1,self.m)
				Y = np.zeros(self.m)
			else:
				X,Y = np.meshgrid(np.linspace(0,1,self.m),np.linspace(0,1,self.n))
			for a in range(len(Hc)):
				i = int(np.round(Hc[a]))
				print i
				for j in xrange(self.T): 
					step = np.loadtxt(self.result_path+'/results_FE_Hc%d_n%04d.txt'%(i,j))
					tmp[j] = np.max(np.linalg.norm(self.exact(X,Y,(j+1)*self.dt)-step))
				self.error.append(tmp.copy())



	def SaveError(self,header=None):
		fname = self.result_path+'/error.txt'
		f = np.asanyarray(self.error)
		if header is None:
			np.savetxt(fname,f)
		else:
			np.savetxt(fname,f)

	def PlotError(self,save=True):
		# mpl.hold('on')
		color = ['b-','r-','k-','c-','g-','m-','b-x','r-x','k-x','c-x','g-x','m-o','b-o','r-o','k-o','c-o','g-o','m-o']
		for i in range(len(self.error)):
			# mpl.plot(np.log(self.error[i]/self.dt))
			mpl.plot(self.error[i],color[i],label='Hc = %d'%(int(np.round(self.legends[i]))))
			# print self.error[i]
		mpl.legend(loc=0)
		mpl.xlabel('timestep no.')
		mpl.ylabel('max(abs(simulation-exact))')
		mpl.title('Error plot; dt = %g'%self.dt)
		if save:
			mpl.savefig(self.result_path+'/errorplot.png')
			mpl.savefig(self.result_path+'/errorplot.eps')
		mpl.show()

	def VerifyDeterministicError(self,save=True,DT=[]):
		if len(DT)==0:
			DT.append(self.dt)
		self.error = []
		leg = []
		tmp = np.zeros(self.T+1)
		color = ['b-','r-','k-','c-','g-','m-','b-x','r-x','k-x','c-x','g-x','m-x']
		if self.n>=1:
			X,Y = np.meshgrid(np.linspace(0,1,self.m),np.linspace(0,1,self.n))
		else:
			X = np.linspace(0,1,self.m)
			Y = np.zeros(self.m)
		for k in range(len(DT)):
			self.dt = DT[k]
			self.RunDeterministic()
			self.walk = []
			for i in range(self.T):
				infile = np.loadtxt(self.result_path+'/Deterministic_n%04d.txt'%i)
				tmp[i] = np.abs(np.linalg.norm(infile-self.exact(X,Y,(i+1)*self.dt)))
			self.error.append(tmp.copy())
			mpl.plot(self.error[k][:-1],color[k])
			leg.append('dt = %g'%DT[k])
		mpl.legend(leg)
		mpl.xlabel('timestep #')
		mpl.ylabel('errornorm')
		mpl.title('dt = %g'%self.dt)
		if save:
			mpl.savefig(self.result_path+'/deterministic_errorplot.png')
			mpl.savefig(self.result_path+'/deterministic_errorplot.eps')
		mpl.show()

	def UpdateWebpageDefault(self):
		datetime = self.datetime
		f = open(this_dir +'/doc/web/index.html','r')
		html = f.read()
		f.close()
		pattern = r'</table>' 	# There is an empty table at the top of the html-file
		match = re.search(pattern,html)
		part = match.span()[-1]+1

		html_code = """<h3><a name="%s" class="anchor" 
		href="#%s"><span class="octicon octicon-link">
		</span></a>New experiment %s.</h3><br>
		Ran new experiment using %d different values of 
		Hc. Results can be found 
		in doc/results/experiment_%s.
		"""%(datetime,datetime,time.ctime(),self.runcounter,self.datetime)
		
		sum_html = html[:part]+html_code+html[part:]
		f = open(this_dir +'/doc/web/index.html','w')
		f.write(sum_html)
		f.close()

	def UpdateSpecial(self):
		datetime = self.datetime
		f = open(this_dir +'/doc/web/index.html','r')
		html = f.read()
		f.close()
		pattern = r'</table>' 	# There is an empty table at the top of the html-file
		match = re.search(pattern,html)
		part = match.span()[-1]+1
		
		title = raw_input('Title on webpage (html is taken care of)>> ')
		html_code = """<h3><a name="%s" class="anchor" 
		href="#%s"><span class="octicon octicon-link">
		</span></a>%s</h3>
		"""%(datetime,datetime,title)
		explanaiton = raw_input('Add description (optional):  ')
		html_code += '<br>'+explanaiton+'\n'

		sum_html = html[:part]+html_code+html[part:]
		f = open(this_dir +'/doc/web/index.html','w')
		f.write(sum_html)
		f.close()


	def Visualize(self,path=None,filename=None,viz_type='difference'):
		if path is None:
			path = self.result_path
		if filename is None:
			filename = '/results'
		im = []
		if self.n<=1:
			fig = mpl.figure()
			x = np.linspace(0,1,self.m)
			counter = 1
			for step in sorted(glob.glob(path+filename+'*.txt')):
				tmp = np.loadtxt(step)
				if viz_type=='difference':
					im.append(mpl.plot(x,(self.exact(x,np.zeros(self.m),counter*self.dt)-tmp),'b-'))
				elif viz_type=='exact':
					im.append(mpl.plot(x,self.exact(x,np.zeros(self.m),counter*self.dt),'-b'))
				else:
					im.append(mpl.plot(x,tmp,'b-'))
				counter += 1
			ani = animation.ArtistAnimation(fig,im)
			mpl.show()
		else:
			X,Y = np.meshgrid(np.linspace(0,1,self.m),np.linspace(0,1,self.n))
			mpl.ion()
			fig = mpl.figure()
			ax = fig.add_subplot(111,projection='3d')
			counter = 1
			for step in sorted(glob.glob(path+filename+'*.txt')):
				# time.sleep(0.3)
				tmp = np.loadtxt(step)
				if viz_type=='difference':
					wframe = ax.plot_wireframe(X,Y,(self.exact(X,Y,(counter*self.dt))-tmp))
				elif viz_type=='exact':
					wframe = ax.plot_wireframe(X,Y,self.exact(X,Y,(counter*self.dt)))
				else:
					wframe = ax.plot_wireframe(X,Y,tmp)
				mpl.draw()
				if counter==1 and viz_type!='difference':
					# pass
					ax.set_autoscaley_on(False)
				ax.collections.remove(wframe)
				counter +=1
			# ax.plot_wireframe(X,Y,self.exact(X,Y,0))
			# mpl.draw()
			# time.sleep(1)

	def ConvergenceTest(self,Hc,save=True):
		E = []
		for i in xrange(len(Hc)):
			E.append(np.sqrt(Hc[i]*np.sum(self.error[i])))
		r = [0]*(len(E)-1)
		for j in xrange(len(E)-1):
			r[j] = np.log(E[j+1]/E[j])/np.log(float(Hc[j+1])/float(Hc[j]))
		print r
		print E
		mpl.plot((Hc[:-1]),np.abs(r),'b-x')
		mpl.xlabel('Hc (conversion rate)')
		mpl.ylabel('r')
		mpl.title('Convergence rate')
		if save:
			mpl.savefig(self.result_path+'/ConvergenceTest.eps')
			mpl.savefig(self.result_path+'/ConvergenceTest.png')
		mpl.show()

	def Compare(self,filename,func):
		counter = 1
		im = []
		error = []
		fig = mpl.figure()
		x = np.linspace(0,1,self.m)
		x,y = np.meshgrid(np.linspace(0,1,self.m),np.linspace(0,1,self.n))
		dx = x[1,1]-x[0,0]
		dy = y[1,1]-y[0,0]
		for i in sorted(glob.glob(self.result_path+filename)):
			infile = np.loadtxt(i)
			im.append(mpl.plot(func(counter,x,y,dx,dy,self.dt)-infile,'b-'))
			error.append(np.linalg.norm(func(counter,x,y,dx,dy,self.dt)-infile))
			counter += 1
		ani = animation.ArtistAnimation(fig,im)
		mpl.show()
		mpl.plot(error[1:])
		mpl.show()

	def SetInitialCondition(self,init):
		np.savetxt('InitialCondition.txt',init)

	def SetDiffusionTensor(self,tensor):
		np.savetxt('DiffusionTensor.txt',tensor)

	def exact(self,x,y,t):
		return 0


	def Finish(self):
		if self.debug:
			os.system('rm -rf %s'%self.parent_path)

# End of class Experiment



def f(x,y,t):
	# return np.exp(-t*np.pi**2)*np.cos(np.pi*x)
	# return np.ones(np.shape(x))*1.5
	# D = v = 1
	# tmp  = (1.0/np.sqrt(4*np.pi*D*t))*np.exp(-(x-v*t)**2/(4*D*t))
	# return tmp/np.sum(tmp)
	return np.exp(-t*np.pi**2)*np.cos(np.pi*x)*np.cos(np.pi*y)

def F(x,y,t):
	tmp = np.zeros(np.shape(x))
	j=0
	for i in xrange(1,500,2):
		tmp[:] += -2*(-1)**j/(i*np.pi)*np.exp(-(i*np.pi)*(i*np.pi)*t)*np.cos(i*np.pi*x)
		j+=1
	tmp[:] += 0.5
	return tmp

from scipy.special import binom

def numerical_exact(n,x,y,dx,dy,dt,D=1):
	u = np.zeros(np.shape(x))
	def u_xx(i,x,y):
		#return (-1)**i*(np.pi**(2*i))*np.cos(np.pi*x)
		if dx**(2*i)>0 and type(2**i)!=type(long()):
			return 2**(i-1)*np.cos(np.pi*y)*np.cos(np.pi*x)*((np.cos(np.pi*dx)-1)**i/(dx**(2*i)) +(np.cos(np.pi*dy)-1)**i/(dy**(2*i)))
		else:
			return 0
	
	for i in xrange(n+1):
		if dt**i>0:
			tmp = binom(n,i)*(D*dt)**i*u_xx(i,x,y)
		else:
			tmp = np.zeros(np.shape(u))
		u[:] += tmp
	return u

def D(x,y,t=0):
	return x
	# return np.ones(np.shape(x))*0.5

if __name__ == '__main__':
	DEBUG = True
	save_files = True
	mode = 'test'


	this_dir = right_split(os.getcwd(),'/')

	x0 = 0.4
	y0 = 0.6
	x1 = 0.6
	y1 = 0.7
	m = 23
	n = 23
	T = 25
	dx = 1.0/(m-1)
	dy = 1.0/(n-1) if n>1 else 0
	dt = dx*dy/4.0 if n>1 else dx**2/5.0
	dt = 0.001

	x,y = np.meshgrid(np.linspace(0,1,m),np.linspace(0,1,n))
	print 'Python: ',dt,' dx: ',dx
	Hc = [1600]
	# Hc = [1400,2000,3200,4400,5600,6800,8000,9200,10400,11600,13000]
	# Hc = [1000,2000,4000,8000,16000,32000,64000,128000,256000,512000,1024000,2048000]
	name = '/home/fredriep/Dropbox/uio/thesis/doc/results/experiment_18102013_1337/results/'

	run = Experiment(this_dir,DEBUG,save_files)
	run.exact = f
	run.SetInitialCondition(f(x,y,0))
	run.SetDiffusionTensor(D(x,y))
	run.compile()
	# run.SetInitialCondition(f(np.linspace(0,1,m),np.zeros(m),0))
	# dt = [dx*dy/5.0*10**(-i) for i in range(6)]
	# dt = [1e-4,1e-5,1e-6,1e-7,1e-8]
	run.SetupRun(x0,x1,y0,y1,m,n,T,dt)
	run.VerifyDeterministicError()

	# for i in Hc:
	# 	print "Hc = %g"%i
	# 	run.RunSimulation(i)
	# time.sleep(1)
	# run.CalculateError(Hc,exact=True)
	# run.PlotError()
	# h = [1./Hc[i] for i in range(len(Hc))]
	# run.ConvergenceTest(h)
	# run.ConvergenceTest(dt)
	# run.Compare('/Deterministic_n*',numerical_exact)

	# run.SaveError(header="max(abs(error)) for manifactured solution u(x,t) = exp(-t*pi**2*cos(pi*x) in 1D. Hc = %g"%Hc[0])
	# run.UpdateSpecial()
	# run.Visualize(viz_type=None)

	# run.Visualize(viz_type='difference')
	run.Visualize(filename='/Deterministic_n',viz_type=None)
	# run.Visualize(filename='/Deterministic_n',viz_type='exact')
	# run.Visualize(filename='/Deterministic_n',viz_type='difference')
	# a = raw_input('press return >>')

	run.Finish()
# g++ *.cpp -O2 -Wall -o main_walk
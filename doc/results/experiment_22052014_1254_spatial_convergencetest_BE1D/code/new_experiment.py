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
	def __init__(self,path, DEBUG=False, save=True, info = ''):
		self.debug = DEBUG
		this_dir = path

		t = time.gmtime()
		self.datetime = '%02d%02d%d_%02d%02d'%(t.tm_mday,t.tm_mon,t.tm_year,t.tm_hour,t.tm_min)
		self.url = 'https://raw.github.com/fepettersen/thesis/master/doc/results'+'/experiment_%s/results'%self.datetime
		if self.debug:
			self.parent_path = this_dir +'/doc/results/experiment_%s%s_DELETE_ME'%(self.datetime,info)
		else:
			self.parent_path = this_dir +'/doc/results/experiment_%s%s'%(self.datetime,info)

		# self.parent_path = this_dir + '/doc/results/experiment_04122013_1259_convergenceTest_combinedSimulation_2d'
		self.code_path = self.parent_path+'/code'
		self.parameter_path = self.parent_path+'/parameters'
		self.result_path = self.parent_path+ '/results'
		self.tofile = 1 if save else 0;
		self.error = []

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


	def SaveParameters(self):
		parameters = open('%s/parameters.txt'%self.parameter_path,'w')
		string = "dt = %g\n"%self.dt
		parameters.write(string)
		string = "length of dendrite = %g micrometers\n"%x_end
		parameters.write(string)
		string = "simulating for %g seconds"%(self.T*self.dt)
		parameters.close()


	def compile(self):
		# os.system('g++ *.cpp -O2 -o -larmadillo -llapack -lblas main_walk')
		# os.system('make')
		os.system('g++ *.cpp -o main_walk -O2 -larmadillo -llapack -lblas')

	def SetupRun(self,x0,x1,y0,y1,m,n,T,dt,filename="a"):
		self.x0 = x0; self.x1 = x1; self.y0 = y0; self.y1 = y1
		self.m = m; self.n = n; self.T = T; self.dt = dt
		self.filename = filename

	def RunDeterministic(self):
		# Run a simulation with the parameters specified in SetupRun 
		# disregarding random walks
		os.system('./main_walk %d %f %f %f %f %d %d %d %s %s %f %g'%(self.tofile,0,0,0,0,self.m,self.n,self.T,self.result_path,"Deterministic",0,self.dt))
		name = self.result_path+'/Deterministic*.txt'
		self.CalculateError(name)

	def RunSimulation(self,Hc):
		# Run a simulation with the parameters specified in SetupRun 
		tofile = self.tofile
		x0 = self.x0; x1 = self.x1; y0 = self.y0; y1 = self.y1
		m = self.m; n = self.n; T = self.T
		filename = self.filename
		result_path = self.result_path
		tmp = [0]*T
		os.system('./main_walk %d %f %f %f %f %d %d %d %s %s %f %g'%(tofile,x0,x1,y0,y1,m,n,T,result_path,filename,Hc,self.dt))
		self.runcounter += 1
		name = self.result_path+'/results_FE_Hc%d*.txt'%Hc if filename=='a' else result_path+filename
		self.CalculateError(name)
		self.SaveParameters()

	def CalculateError(self,filename):
		if self.n <= 1:
			X = np.linspace(0,1,self.m)
			Y = np.zeros(self.m)
		else:
			X,Y = np.meshgrid(np.linspace(0,1,self.m),np.linspace(0,1,self.n))
		exact = True if self.exact(X,Y,0)!=None else False
		if exact:
			tmp = np.zeros(self.T)
			j=0
			for step in sorted(glob.glob(filename)):
				# print step
				infile = np.loadtxt(step)
				tmp[j] = np.linalg.norm(self.exact(X,Y,(j+1)*self.dt)-infile)
				j+=1
			self.error.append(tmp.copy())


	def SaveError(self,header=None):
		fname = self.result_path+'/error.txt'
		f = np.asanyarray(self.error)
		if header is None:
			np.savetxt(fname,f)
		else:
			np.savetxt(fname,f)

	def PlotError(self,legend,save=True):
		# mpl.hold('on')
		color = ['b-','r--','k:','c-','g-','m-','b-x','r-x','k-x','c-x','g-x','m-o','b-o','r-o','k-o','c-o','g-o','m-o']
		if self.runcounter!=len(self.error):
			# print self.error,'\n \n'
			# print self.runcounter
			mpl.plot(self.error[0],color[0],label='Deterministic')
			for i in range(self.runcounter):
				# mpl.plot(np.log(self.error[i]/self.dt))
				mpl.plot(self.error[i+1],color[i+1],label=legend[i])
		else:
			for i in range(self.runcounter):
				# mpl.plot(np.log(self.error[i]/self.dt))
				mpl.plot(self.error[i],color[i],label=legend[i])
		mpl.legend(loc=0)
		mpl.xlabel('timestep no.')
		mpl.ylabel('errornorm (l2)')
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
		tmp = np.zeros(self.T)
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
			# # for i in range(self.T):
			# # 	infile = np.loadtxt(self.result_path+'/Deterministic_n%04d.txt'%i)
			# # 	tmp[i] = np.abs(np.linalg.norm(infile-self.exact(X,Y,(i+1)*self.dt)))
			# # self.error.append(tmp.copy())
			# mpl.plot(self.error[k][:],color[k])
			leg.append('dt = %g'%DT[k])
		mpl.legend(leg)
		mpl.xlabel('timestep #')
		mpl.ylabel('errornorm')
		mpl.title('dt = %g'%self.dt)
		if save:
			mpl.savefig(self.result_path+'/deterministic_errorplot.png')
			mpl.savefig(self.result_path+'/deterministic_errorplot.eps')
		# mpl.show()

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

	def UpdateWebpageSpecial(self):
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


	def Visualize(self,path=None,filename=None,viz_type='difference',save_video=False):
		if path is None:
			path = self.result_path
		if filename is None:
			filename = '/results'
		im = []
		if self.n<=1:
			# 1d
			fig = mpl.figure()
			x = np.linspace(x_start,x_end,self.m)
			counter = 1
			for step in sorted(glob.glob(path+filename+'*.txt')):
				tmp = np.loadtxt(step)
				mpl.title('t=%.3f/%.3f'%(counter*self.dt,self.T*self.dt))
				if viz_type=='difference':
					im.append(mpl.plot(x,(self.exact(x,np.zeros(self.m),counter*self.dt)-tmp),'b-'))
				elif viz_type=='exact':
					im.append(mpl.plot(x,self.exact(x,np.zeros(self.m),counter*self.dt),'-b'))
				else:
					im.append(mpl.plot(x,tmp,'b-'))
				if save_video:
					fig.savefig(step.split('.')[0]+'.png')
				counter += 1
			ani = animation.ArtistAnimation(fig,im)
			if save_video:
				self.MakeVideo(sorted(glob.glob(path+filename+'*.png')))
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
				if save_video:
					fig.savefig(step.split('.')[0]+'.png')
				ax.collections.remove(wframe)
				counter +=1
			if save_video:
				self.MakeVideo(sorted(glob.glob(path+filename+'*.png')))
					
	def MakeVideo(self,files):
		import scitools.std as sci
		sci.movie(files, encoder='mencoder', fps=25, output_file=self.result_path+'/movie.mpeg')
		for i in files:
			os.system('rm %s'%i)



	def ConvergenceTest(self,h,save=True):
		# h is the discretization parameter which the convergence is beeing tested for
		E = []
		if self.runcounter != len(self.error):
			for i in xrange(len(h)):
				E.append(np.sqrt(h[i]*np.sum(self.error[i+1]**2)))
		else:
			for i in xrange(len(h)):
				E.append(np.sqrt(h[i]*np.sum(self.error[i]**2)))
				# E.append(np.sqrt(h[i]*np.max(self.error[i])**2'))
		r = [0]*(len(E)-1)
		for j in xrange(len(E)-1):
			r[j] = np.log(E[j+1]/E[j])/np.log(float(h[j+1])/float(h[j]))
		print r
		print E
		# mpl.ylim(0.0,1.1)
		mpl.plot((h[:-1]),r,'b-x')
		mpl.xlabel('dx')
		mpl.ylabel('r')
		mpl.title('Convergence rate')
		if save:
			mpl.savefig(self.result_path+'/ConvergenceTest.eps')
			mpl.savefig(self.result_path+'/ConvergenceTest.png')
		mpl.show()

	def Compare(self,filename,func,mat=None,U0=None):
		counter = 1
		im = []
		error = []
		fig = mpl.figure()
		x = np.linspace(0,1,self.m)
		y = np.zeros(np.shape(x))
		if mat!=None and U0!=None:
			temp = mat
			U0 = U0.reshape(-1)
			print 'dt = ',self.dt
			for i in sorted(glob.glob(self.result_path+filename)):
				infile = np.loadtxt(i)
				temp = np.linalg.matrix_power(mat,counter)
				# temp = np.dot(mat,temp)
				infile = infile.reshape(-1)
				error.append(np.linalg.norm(infile-np.dot(temp,U0)))
				# im.append(mpl.plot(infile-np.dot(temp,U0),'b-'))
				counter += 1
		else:
			# x,y = np.meshgrid(np.linspace(0,1,self.m),np.linspace(0,1,self.n))
			# dx = x[1,1]-x[0,0]
			# dy = y[1,1]-y[0,0]
			dx = x[1]-x[0]
			dy = 1
			for i in sorted(glob.glob(self.result_path+filename)):
				infile = np.loadtxt(i)
				# im.append(mpl.plot(func(counter,x,y,dx,dy,self.dt),'b-'))
				error.append(np.linalg.norm(func(counter,x,y,dx,dy,self.dt)-infile))
				counter += 1
		# ani = animation.ArtistAnimation(fig,im)
		# mpl.show()

		mpl.plot(error[1:])
		mpl.title('Error compared to numerical exact')
		mpl.xlabel('timestep')
		mpl.ylabel('errornorm')
		mpl.show()

	def SetInitialCondition(self,init):
		# save the initial condition (init) as a numpy-array
		np.savetxt('InitialCondition.txt',init,fmt = '%.18f')

	def SetDiffusionTensor(self,tensor):
		# save the anisotropic diffusion "tensor" as a numpy-array
		np.savetxt('DiffusionTensor.txt',tensor)

	def exact(self,x,y,t):
		# This function can be set by user:
		# instance.exact = your_own_function
		return 0


	def Finish(self):
		if self.debug:
			print 'Test done! \n Removing folder \n "%s" and all contents.\n'%self.parent_path
			os.system('rm -rf %s'%self.parent_path)
		else:
			os.system('cp spine_info.txt %s'%self.parameter_path)
			print 'Done! all results and copy of code with parameters in folder \n \n %s \n'%self.parent_path

# End of class Experiment



def f(x,y,t):
	# return np.exp(-t*np.pi**2)*np.cos(np.pi*x)
	# return np.ones(np.shape(x))*1.5
	# D = v = 1
	# tmp  = (1.0/np.sqrt(4*np.pi*D*t))*np.exp(-(x-v*t)**2/(4*D*t))
	# return tmp/np.sum(tmp)
	return np.exp(-t*np.pi**2)*np.cos(np.pi*x)*np.cos(np.pi*y) +1
	# return np.pi*np.ones(np.shape(x))
	# return x+np.pi*t

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
		# return (-1)**i*(np.pi**(2*i))*np.cos(np.pi*x)
		if dx**(2*i)>0 and type(2**i)!=type(long()):
			# return 2**(i-1)*np.cos(np.pi*y)*np.cos(np.pi*x)*((np.cos(np.pi*dx)-1)**i/(dx**(2*i)) +(np.cos(np.pi*dy)-1)**i/(dy**(2*i)))
			return 2**(i)*(np.cos(np.pi*dx)-1)**i/(dx**(2*i))
		else:
			return 0
	
	for i in xrange(n+1):
		if dt**i>0:
			tmp = binom(n,i)*(D*dt)**i*u_xx(i,x,y)
		else:
			tmp = np.zeros(np.shape(u))
		u[:] += tmp
	return np.cos(np.pi*x)*u

def D(x,y,t=0):
	# return np.pi*x
	return np.ones(np.shape(x))

def GaussianPulse(x,y,t=0,x0=0,sigma=1.0,A=2.5,Hc=15):
	# return A*np.exp(-(x-x0)**2/(2*sigma**2)) + 1.0/((x+1.5)*Hc) +0.3*np.random.rand(len(x))
	return A*np.exp(-(x-x0)**2/(2*sigma**2))

	
if __name__ == '__main__':
	DEBUG = False
	save_files = True
	mode = 'test'


	this_dir = right_split(os.getcwd(),'/')

	x0 = 0.0
	y0 = 0.0
	x1 = 1.0
	y1 = 1.0
	m = 1800
	n = 1
	T = 9000		# no.of timesteps, [dt*T] = seconds

	x_start = 0
	x_end = 50.0 		#um

	dx = (x_end-x_start)/(m-1.0)
	dy = 1.0/(n-1) if n>1 else 0
	# dt = dx*dy/4.0 if n>1 else dx**2/5.0
	bla = 0.3*dx**2
	# dt = [bla]
	dt = [0.01]


	if n>1:
		x,y = np.meshgrid(np.linspace(0,1,m),np.linspace(0,1,n))
	else:
		x = np.linspace(x_start,x_end,m)
		y = np.zeros(m)
	# Hc = [1500]
	Hc = [20]#,2000,20000]
	# Hc = [5600, 10000, 50000]
	info='_Testrun_for_PKCg_diffusion'
	info='_spatial_convergencetest_BE1D'
	run = Experiment(this_dir,DEBUG,save_files,info)
	run.exact = f

	run.SetInitialCondition(run.exact(x,y,0))
	# run.SetDiffusionTensor(D(x,y))
	# run.SetInitialCondition(GaussianPulse(x,y,0))
	run.SetDiffusionTensor(D(x,y))

	run.compile()
	run.SetupRun(x0,x1,y0,y1,m,n,T,dt[0])
	# run.VerifyDeterministicError()
	# run.RunSimulation(Hc[0])
	# M = np.loadtxt("BE_matrix_inverse.txt")
	# u0 = run.exact(x,y,0)
	# run.Compare('/results_FE_Hc*.txt',numerical_exact)				# FE version
	# run.Compare('/Deterministic_n*.txt',numerical_exact,M,u0)				# BE version
	# run.PlotError('dt = %g'%dt[0])
	### --- Run for walkers --- ###

	# for i in Hc:
	# 	print "Hc = %g"%i
	# 	run.SetupRun(x0,x1,y0,y1,m,n,T,dt[0])
	# 	run.RunSimulation(i)
	# run.ConvergenceTest(Hc)
	# leg = ['Hc = %g'%i for i in Hc]
	# # os.system('python spine_statistics.py')

	## --- Run for time-step --- ###
	# dt = [0.05,0.02,0.01,0.005]

	# for i in dt:
	# 	print "dt = %g"%i
	# 	run.SetupRun(x0,x1,y0,y1,m,n,T,i)
	# 	run.RunSimulation(Hc[0])
	# run.ConvergenceTest(dt)
	# leg = ['dt = %g'%i for i in dt]
	# run.PlotError(leg)

	## --- Run for h --- ###
	dt = 6e-5
	h = [0.1,0.05,0.02]
	# dx = 0.025
	# h = [0.1,0.05,0.02,0.025]
	# h = [0.025]
	for i in h:
		# i = 1.0/j
		timestep = dt
		# dt.append(i*i)
		m = (1/i)+1
		# m = 21
		n = 1
		if n==1:
			x = np.linspace(0,1,m)
			y = np.zeros(m)
		elif n==m:
			x,y = np.meshgrid(np.linspace(0,1,m),np.linspace(0,1,m))
		run.SetupRun(x0,x1,y0,y1,m,n,T,timestep)
		run.SetInitialCondition(run.exact(x,y,0))
		run.SetDiffusionTensor(D(x,y))
		# hc = int(round(1.0/(timestep**2)))
		hc = 100
		print m," , ", timestep, " , ", hc
		run.RunSimulation(1*hc)
	leg = ['dx = %g'%i for i in h]
	run.ConvergenceTest(h)
	# leg = ['dt = %g'%i for i in dt]
	run.PlotError(leg)
	# # print 'dx =',dx
	# # run.SaveError(header="max(abs(error)) for manufactured solution u(x,t) = exp(-t*pi**2*cos(pi*x) in 1D. Hc = %g"%Hc[0])
	# # run.UpdateWebpageSpecial()
	# run.Visualize(viz_type=None)
	# run.Visualize(viz_type='exact')

	# run.Visualize(viz_type='difference',save_video=True)
	# run.Visualize(filename='/Deterministic_n',viz_type=None)
	# run.Visualize(filename='/Deterministic_n',viz_type='exact')
	# run.Visualize(filename='/Deterministic_n',viz_type='difference')
	# a = raw_input('press return >>')

	run.Finish()
